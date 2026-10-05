"""Tests for cemd.analysis.util and cemd.analysis.velocities.

Synthetic trajectories carry the answer in the way they are built; the
tests on ``test/traj_solution`` (untracked, half a gigabyte, skipped when
absent) compare each function with plain numpy on the same frames.

``velocity_profile`` cannot be run on ``traj_solution`` at all: a DCD file
stores positions only, so ``sel.velocities`` raises ``NoDataError``. Its
tests therefore build velocities.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

mda = pytest.importorskip("MDAnalysis")
from MDAnalysis.coordinates.memory import MemoryReader  # noqa: E402

from cemd.analysis.util import (  # noqa: E402
    com,
    mean_pos,
    minmax_position,
    require_box,
    shift2com,
    write_dcd,
)
from cemd.analysis.velocities import velocity_profile  # noqa: E402

SOLUTION = Path(__file__).parents[1] / "traj_solution"
requires_solution = pytest.mark.skipif(
    not (SOLUTION / "traj.dcd").exists(), reason="test/traj_solution not present"
)

BOX = 20.0


def memory_universe(
    positions, types=None, masses=None, velocities=None, box=BOX
) -> mda.Universe:
    """A universe from a (frames, atoms, 3) array, in a cubic box."""
    positions = np.asarray(positions, dtype=np.float32)
    n = positions.shape[1]
    universe = mda.Universe.empty(
        n,
        n_residues=n,
        atom_resindex=np.arange(n),
        residue_segindex=np.zeros(n, dtype=int),
        trajectory=True,
    )
    universe.add_TopologyAttr("type", types or ["A"] * n)
    universe.add_TopologyAttr("mass", masses or [1.0] * n)
    kwargs = {} if velocities is None else {"velocities": velocities}
    universe.load_new(positions, format=MemoryReader, **kwargs)
    for ts in universe.trajectory:
        ts.dimensions = [box, box, box, 90.0, 90.0, 90.0]
    return universe


def random_universe(n_atoms=50, n_frames=12, seed=0, types=None):
    rng = np.random.default_rng(seed)
    return memory_universe(
        rng.uniform(0, BOX, (n_frames, n_atoms, 3)),
        types=types,
        masses=list(rng.uniform(1, 20, n_atoms)),
    )


# ---------------------------------------------------------------------------
# require_box
# ---------------------------------------------------------------------------


class _Source:
    def __init__(self, dimensions):
        self.dimensions = dimensions


@pytest.mark.parametrize(
    "dimensions",
    [
        None,
        [0.0, 0.0, 0.0, 90.0, 90.0, 90.0],
        [10.0, -5.0, 10.0, 90.0, 90.0, 90.0],
        [10.0, np.nan, 10.0, 90.0, 90.0, 90.0],
        [10.0, np.inf, 10.0, 90.0, 90.0, 90.0],
    ],
)
def test_require_box_refuses_a_missing_or_degenerate_box(dimensions):
    with pytest.raises(ValueError, match="needs the periodic box"):
        require_box(_Source(dimensions), "somebody")


def test_require_box_names_the_caller_and_returns_the_box():
    box = [10.0, 20.0, 30.0, 90.0, 90.0, 90.0]
    assert require_box(_Source(box), "rdf") == pytest.approx(box)
    with pytest.raises(ValueError, match="^rdf needs"):
        require_box(_Source(None), "rdf")


# ---------------------------------------------------------------------------
# write_dcd
# ---------------------------------------------------------------------------


def test_write_dcd_round_trips_a_selection_and_a_frame_slice(tmp_path, capsys):
    universe = random_universe(types=["A", "B"] * 25)
    out = tmp_path / "sub.dcd"

    write_dcd(universe, str(out), selection="type B", start=2, end=10, step=3)

    selection = universe.select_atoms("type B")
    written = mda.Merge(selection)
    written.load_new(str(out))
    frames = list(range(2, 10, 3))
    assert len(written.trajectory) == len(frames)
    assert written.dimensions == pytest.approx([BOX, BOX, BOX, 90, 90, 90])
    for k, frame in enumerate(frames):
        universe.trajectory[frame]
        written.trajectory[k]
        assert written.atoms.positions == pytest.approx(selection.positions, abs=1e-4)


def test_write_dcd_without_an_end_writes_to_the_last_frame(tmp_path):
    universe = random_universe(n_frames=12)
    out = tmp_path / "tail.dcd"

    write_dcd(universe, str(out), start=9)

    written = mda.Merge(universe.atoms)
    written.load_new(str(out))
    assert len(written.trajectory) == 3


def test_write_dcd_refuses_a_wrong_extension_or_an_empty_selection(tmp_path):
    universe = random_universe()
    with pytest.raises(TypeError, match=r"\.dcd"):
        write_dcd(universe, str(tmp_path / "x.xtc"))
    with pytest.raises(ValueError, match="matched 0 atoms"):
        write_dcd(universe, str(tmp_path / "x.dcd"), selection="type Z")


# ---------------------------------------------------------------------------
# shift2com
# ---------------------------------------------------------------------------


def test_shift2com_puts_the_reference_centre_of_mass_at_the_origin(tmp_path):
    universe = random_universe(types=["A"] * 10 + ["B"] * 40)
    out = tmp_path / "recentred.dcd"

    shift2com(universe, ["A"], str(out))

    recentred = mda.Merge(universe.atoms)
    recentred.load_new(str(out))
    reference = recentred.select_atoms("type A")
    assert len(recentred.trajectory) == 12
    for _ in recentred.trajectory:
        assert reference.center_of_mass() == pytest.approx(np.zeros(3), abs=1e-4)


def test_shift2com_is_mass_weighted(tmp_path):
    # A 1-amu atom at x=0 and a 3-amu atom at x=4: COM x=3.
    universe = memory_universe([[[0.0, 10, 10], [4.0, 10, 10]]] * 2, masses=[1.0, 3.0])
    out = tmp_path / "recentred.dcd"

    shift2com(universe, ["A"], str(out))

    back = mda.Merge(universe.atoms)
    back.load_new(str(out))
    assert back.atoms.positions[:, 0] == pytest.approx([-3.0, 1.0], abs=1e-5)


def test_shift2com_says_so_when_the_reference_matches_nothing(tmp_path):
    with pytest.raises(ValueError, match="No atoms found"):
        shift2com(random_universe(), ["Z"], str(tmp_path / "x.dcd"))


# ---------------------------------------------------------------------------
# minmax_position
# ---------------------------------------------------------------------------


def test_minmax_position_averages_the_extremes_frame_by_frame(capsys):
    universe = random_universe(types=["A"] * 25 + ["B"] * 25)
    b = universe.select_atoms("type B")

    lows, highs = [], []
    for ts in universe.trajectory[0:8]:
        lows.append(b.positions[:, 1].min())
        highs.append(b.positions[:, 1].max())

    low, high = minmax_position(universe, ["B"], axis="y", start=0, end=8)
    assert low == pytest.approx(np.mean(lows))
    assert high == pytest.approx(np.mean(highs))


def test_minmax_position_takes_several_types_and_each_axis(capsys):
    universe = random_universe(types=["A"] * 10 + ["B"] * 10 + ["C"] * 30)
    ab = universe.select_atoms("type A B")

    for axis, column in (("x", 0), ("y", 1), ("z", 2)):
        expected = np.mean(
            [ab.positions[:, column].min() for _ in universe.trajectory[0:5]]
        )
        low, _ = minmax_position(universe, ["A", "B"], axis=axis, start=0, end=5)
        assert low == pytest.approx(expected)


def test_minmax_position_refuses_an_unknown_axis(capsys):
    with pytest.raises(ValueError):
        minmax_position(random_universe(), ["A"], axis="w", start=0, end=3)


# ---------------------------------------------------------------------------
# com
# ---------------------------------------------------------------------------


def test_com_is_the_mean_over_frames_of_the_mass_weighted_centre(capsys):
    universe = random_universe(types=["A"] * 20 + ["B"] * 30)
    a = universe.select_atoms("type A")
    expected = np.mean([a.center_of_mass() for _ in universe.trajectory], axis=0)

    assert com(universe, ["A"]) == pytest.approx(expected)


def test_com_unwrap_keeps_a_group_straddling_the_boundary_whole(capsys):
    # x = 0.5 and 19.5 in a 20 A box: 1 A apart, centred on x = 0 (== 20).
    # Without unwrap, the plain average puts it in the middle of the box.
    universe = memory_universe([[[0.5, 10, 10], [19.5, 10, 10]]] * 3)

    x = com(universe, ["A"], unwrap=True)[0]
    # Distance from the true centre, modulo the box.
    assert min(x % BOX, BOX - x % BOX) == pytest.approx(0.0, abs=1e-4)
    assert com(universe, ["A"])[0] == pytest.approx(BOX / 2)


def test_com_unwrap_follows_a_group_drifting_through_the_boundary(capsys):
    # A 1 A dimer drifting +0.5 A per frame from x = 18.5 to x = 22.5,
    # wrapped: its centre crosses the boundary halfway through. Its mean
    # centre is x = 20.5 (== 0.5), not an average of 19-ish and 1-ish.
    centres = 18.5 + 0.5 * np.arange(9)
    frames = [
        [[(c - 0.5) % BOX, 10, 10], [(c + 0.5) % BOX, 10, 10]] for c in centres
    ]
    x = com(memory_universe(frames), ["A"], unwrap=True)[0]
    assert x == pytest.approx(0.5, abs=1e-4)


def test_shift2com_unwrap_recentres_on_a_group_cut_by_the_boundary(tmp_path):
    universe = memory_universe([[[0.5, 10, 10], [19.5, 10, 10]]] * 2)
    out = tmp_path / "recentred.dcd"

    shift2com(universe, ["A"], str(out), unwrap=True)

    back = mda.Merge(universe.atoms)
    back.load_new(str(out))
    x = np.sort(back.atoms.positions[:, 0] % BOX)
    assert x == pytest.approx([0.5, 19.5], abs=1e-4)
    assert back.atoms.positions[:, 1] == pytest.approx([0.0, 0.0], abs=1e-4)


# ---------------------------------------------------------------------------
# mean_pos
# ---------------------------------------------------------------------------


def bonded_dimers(n_frames=30, seed=1):
    """Two-atom molecules, 1.2 A long, one of them across the boundary."""
    rng = np.random.default_rng(seed)
    centres = np.array([[10.0, 10.0, 10.0], [0.0, 5.0, 5.0]])
    positions = np.zeros((n_frames, 4, 3))
    positions[:, [0, 2]] = centres + [-0.6, 0, 0]
    positions[:, [1, 3]] = centres + [0.6, 0, 0]
    positions += rng.normal(0.0, 0.02, positions.shape)
    positions %= BOX  # wrap: the second dimer is split by the boundary

    universe = memory_universe(positions, types=["A", "B", "A", "B"])
    universe.add_TopologyAttr("bonds", [(0, 1), (2, 3)])
    return universe


def test_mean_pos_reports_one_row_per_atom_with_its_type():
    universe = bonded_dimers()
    frame = mean_pos(universe)

    assert list(frame.columns) == ["type", "x", "y", "z"]
    assert list(frame["type"]) == ["A", "B", "A", "B"]


def test_mean_pos_keeps_molecules_whole_across_the_boundary():
    # The second dimer is split by the boundary in every frame; unwrapped,
    # its two atoms stay 1.2 A apart instead of ~ a box length.
    frame = mean_pos(bonded_dimers())

    xyz = frame[["x", "y", "z"]].astype(float).to_numpy()
    assert np.linalg.norm(xyz[1] - xyz[0]) == pytest.approx(1.2, abs=0.05)
    assert np.linalg.norm(xyz[3] - xyz[2]) == pytest.approx(1.2, abs=0.05)


def test_mean_pos_coordinates_are_numeric_columns():
    frame = mean_pos(bonded_dimers())
    assert all(frame[c].dtype.kind == "f" for c in "xyz")


def test_mean_pos_can_be_called_twice_on_one_universe():
    # Regression test: the unwrap transformation stayed attached to the
    # universe, and the second call raised "Can't add transformations again".
    universe = bonded_dimers()
    first = mean_pos(universe)
    second = mean_pos(universe)
    assert second[["x", "y", "z"]].to_numpy() == pytest.approx(
        first[["x", "y", "z"]].to_numpy()
    )


# ---------------------------------------------------------------------------
# velocity_profile
# ---------------------------------------------------------------------------


def two_layer_velocities(n_atoms=400, n_frames=3, seed=2):
    """Atoms below z = 10 move at 1 A/ps per axis, those above at 3."""
    rng = np.random.default_rng(seed)
    positions = rng.uniform(0, BOX, (n_frames, n_atoms, 3))
    speed = np.where(positions[..., 2] < BOX / 2, 1.0, 3.0)
    velocities = np.repeat(speed[..., None], 3, axis=2)
    return memory_universe(
        positions, types=["Ow"] * n_atoms, velocities=velocities.astype(np.float32)
    )


def test_velocity_profile_needs_velocities_and_says_so():
    universe = memory_universe(np.random.default_rng(0).uniform(0, BOX, (2, 10, 3)))
    with pytest.raises(mda.exceptions.NoDataError):
        velocity_profile(universe, "A")


def test_velocity_profile_takes_a_list_of_types():
    profile = velocity_profile(two_layer_velocities(), ["Ow"], bin_size=2.0)
    assert len(profile) > 0


def test_velocity_profile_follows_the_speed_of_each_layer():
    # Lower half moves three times slower than the upper half, whichever
    # way the per-component speeds are averaged.
    profile = velocity_profile(two_layer_velocities(), "Ow", bin_size=2.0)

    lower = profile[profile.index < BOX / 2 - 1].mean()
    upper = profile[profile.index > BOX / 2 + 1].mean()
    assert lower / upper == pytest.approx(1.0 / 3.0, rel=0.05)


# ---------------------------------------------------------------------------
# On the reference trajectory
# ---------------------------------------------------------------------------

OW = "8"


@pytest.fixture(scope="module")
def solution_head(tmp_path_factory):
    """The first 30 frames of traj_solution, as a small DCD of its own."""
    source = mda.Universe(str(SOLUTION / "solution.data"), str(SOLUTION / "traj.dcd"))
    head = tmp_path_factory.mktemp("head") / "head.dcd"
    write_dcd(source, str(head), start=0, end=30)
    return str(SOLUTION / "solution.data"), str(head)


@requires_solution
def test_write_dcd_reproduces_the_reference_frames_exactly(tmp_path):
    topology = str(SOLUTION / "solution.data")
    source = mda.Universe(topology, str(SOLUTION / "traj.dcd"))
    out = tmp_path / "sub.dcd"

    write_dcd(source, str(out), selection=f"type {OW}", start=100, end=110, step=2)

    water = source.select_atoms(f"type {OW}")
    written = mda.Merge(water)
    written.load_new(str(out))
    assert written.dimensions == pytest.approx([35.0, 35.0, 35.0, 90, 90, 90])
    for k, frame in enumerate(range(100, 110, 2)):
        source.trajectory[frame]
        written.trajectory[k]
        assert written.atoms.positions == pytest.approx(water.positions, abs=1e-4)


@requires_solution
def test_minmax_position_on_the_reference_matches_numpy(solution_head, capsys):
    universe = mda.Universe(*solution_head)
    oxygens = universe.select_atoms(f"type {OW}")
    zs = [oxygens.positions[:, 2] for _ in universe.trajectory[0:20]]

    low, high = minmax_position(universe, [OW], axis="z", start=0, end=20)

    assert low == pytest.approx(np.mean([z.min() for z in zs]))
    assert high == pytest.approx(np.mean([z.max() for z in zs]))
    # Wrapped water fills the box, which runs from -0.51 to 34.49.
    assert low == pytest.approx(-0.5, abs=0.2)
    assert high == pytest.approx(34.5, abs=0.2)


@requires_solution
def test_com_of_the_reference_water_sits_at_the_box_centre(solution_head, capsys):
    universe = mda.Universe(*solution_head)
    centre = com(universe, [OW])

    # Uniform liquid in a box spanning [-0.5, 34.5]: the mean COM is the
    # middle of it, within the noise of 1433 molecules.
    assert centre == pytest.approx([17.0, 17.0, 17.0], abs=0.5)


@requires_solution
def test_shift2com_on_the_reference_recentres_every_frame(solution_head, tmp_path):
    topology, head = solution_head
    universe = mda.Universe(topology, head)
    out = tmp_path / "recentred.dcd"

    shift2com(universe, ["5", "9"], str(out))  # sodium and sulfur

    recentred = mda.Universe(topology, str(out))
    ions = recentred.select_atoms("type 5 9")
    assert len(recentred.trajectory) == 30
    for _ in recentred.trajectory:
        assert ions.center_of_mass() == pytest.approx(np.zeros(3), abs=1e-4)
