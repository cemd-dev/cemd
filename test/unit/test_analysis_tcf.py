"""Tests for cemd.analysis.tcf: bond survival and lifetime fits.

Two kinds of reference. On a synthetic trajectory the survival curve is
known exactly, because the moment each pair separates is chosen. On
``test/traj_solution`` (untracked, half a gigabyte, so those tests skip
themselves when it is absent) the reference is a brute-force
implementation written here from ``distance_array`` and nothing else, plus
one physical fact: SPC/E water is rigid, so an O-H bond never breaks.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

mda = pytest.importorskip("MDAnalysis")
from MDAnalysis.coordinates.memory import MemoryReader  # noqa: E402
from MDAnalysis.lib.distances import distance_array  # noqa: E402

from cemd.analysis.tcf import bondcorr, lifetime  # noqa: E402

SOLUTION = Path(__file__).parents[1] / "traj_solution"
requires_solution = pytest.mark.skipif(
    not (SOLUTION / "traj.dcd").exists(), reason="test/traj_solution not present"
)

BOX = 50.0


# bondcorr places blocks `gaplength` frames apart and needs at least one, so
# a single block starting at frame 0 means gaplength = n_frames - corrlength
# (40 frames throughout, unless a test says otherwise).


def pair_universe(break_frames, n_frames=40, box=BOX) -> mda.Universe:
    """Pairs of A and B atoms that separate at a frame of our choosing.

    Pair k sits 1 A apart until frame ``break_frames[k]``, then 3 A apart.
    Pairs are 10 A from one another, so the only bonds are the intended
    ones. A pair with ``break_frames[k] is None`` never separates.
    """
    n_pairs = len(break_frames)
    positions = np.zeros((n_frames, 2 * n_pairs, 3))
    for k, brk in enumerate(break_frames):
        base = np.array([5.0 + 10.0 * (k % 4), 5.0 + 10.0 * (k // 4), 25.0])
        for frame in range(n_frames):
            gap = 1.0 if brk is None or frame < brk else 3.0
            positions[frame, k] = base
            positions[frame, n_pairs + k] = base + [gap, 0.0, 0.0]

    n = 2 * n_pairs
    universe = mda.Universe.empty(
        n,
        n_residues=n,
        atom_resindex=np.arange(n),
        residue_segindex=np.zeros(n, dtype=int),
        trajectory=True,
    )
    universe.add_TopologyAttr("type", ["A"] * n_pairs + ["B"] * n_pairs)
    universe.add_TopologyAttr("mass", [1.0] * n)
    universe.load_new(positions.astype(np.float32), format=MemoryReader)
    for ts in universe.trajectory:
        ts.dimensions = [box, box, box, 90.0, 90.0, 90.0]
    return universe


# ---------------------------------------------------------------------------
# bondcorr on a trajectory whose answer is known
# ---------------------------------------------------------------------------


def test_bondcorr_survival_matches_the_chosen_break_times():
    # Ten pairs, one breaking at each of frames 1..10: after t frames,
    # exactly the pairs with break frame > t are left.
    breaks = list(range(1, 11))
    universe = pair_universe(breaks)

    tcf = bondcorr(
        universe, "A", "B", distance=2.0, dt=1000.0, corrlength=12, gaplength=28
    )

    expected = [sum(b > t for b in breaks) / len(breaks) for t in range(12)]
    assert tcf["TCF"].to_numpy() == pytest.approx(expected, abs=1e-6)


def test_bondcorr_time_axis_is_in_picoseconds():
    universe = pair_universe([None] * 4)
    tcf = bondcorr(
        universe, "A", "B", distance=2.0, dt=250.0, corrlength=10, gaplength=30
    )
    # 250 fs between frames is 0.25 ps.
    assert tcf["t [ps]"].to_numpy() == pytest.approx(np.arange(10) * 0.25)


def test_bondcorr_never_counts_a_bond_that_re_forms():
    # A pair that separates and comes back must stay dead: the TCF is a
    # survival probability, not a bond-present fraction.
    universe = pair_universe([5, None])
    positions = universe.trajectory.timeseries(order="fac").copy()
    positions[8:, 2, 0] = positions[8:, 0, 0] + 1.0  # pair 0's B atom returns
    universe.load_new(positions, format=MemoryReader)
    for ts in universe.trajectory:
        ts.dimensions = [BOX, BOX, BOX, 90.0, 90.0, 90.0]

    tcf = bondcorr(
        universe, "A", "B", distance=2.0, dt=1000.0, corrlength=12, gaplength=28
    )
    assert tcf["TCF"].iloc[4] == pytest.approx(1.0)
    assert tcf["TCF"].iloc[5:].to_numpy() == pytest.approx(0.5)


def test_bondcorr_of_bonds_that_never_break_is_one_throughout():
    universe = pair_universe([None] * 6)
    tcf = bondcorr(
        universe, "A", "B", distance=2.0, dt=1000.0, corrlength=20, gaplength=20
    )
    assert tcf["TCF"].to_numpy() == pytest.approx(1.0)


def test_bondcorr_of_like_atoms_ignores_self_pairs():
    # Every pair separates, so the TCF must end at zero. Two A atoms 1 A
    # apart, both then leaving.
    n_frames = 10
    positions = np.zeros((n_frames, 2, 3)) + 25.0
    positions[:, 1, 0] += 1.0
    positions[5:, 1, 0] += 5.0
    universe = mda.Universe.empty(
        2, n_residues=2, atom_resindex=[0, 1], residue_segindex=[0, 0], trajectory=True
    )
    universe.add_TopologyAttr("type", ["A", "A"])
    universe.add_TopologyAttr("mass", [1.0, 1.0])
    universe.load_new(positions.astype(np.float32), format=MemoryReader)
    for ts in universe.trajectory:
        ts.dimensions = [BOX, BOX, BOX, 90.0, 90.0, 90.0]

    tcf = bondcorr(
        universe, "A", "A", distance=2.0, dt=1000.0, corrlength=8, gaplength=2
    )
    assert tcf["TCF"].iloc[-1] == pytest.approx(0.0)


def test_bondcorr_honours_nblocks():
    # Pair 0 separates at frame 12, pair 1 at frame 25. With blocks at
    # frames 0, 10 and 20 the curve at t=20 averages 0.5, 0 and 0; the
    # first block alone gives 0.5.
    universe = pair_universe([12, 25], n_frames=60)
    tcf = bondcorr(
        universe,
        "A",
        "B",
        distance=2.0,
        dt=1000.0,
        nblocks=1,
        corrlength=30,
        gaplength=10,
    )
    assert tcf["TCF"].iloc[20] == pytest.approx(0.5)


def test_bondcorr_says_so_when_the_correlation_length_leaves_no_block():
    universe = pair_universe([None], n_frames=40)
    with pytest.raises(ValueError):
        bondcorr(universe, "A", "B", distance=2.0, corrlength=39, gaplength=10)


# ---------------------------------------------------------------------------
# lifetime
# ---------------------------------------------------------------------------


def double_exponential(a1, tau1, tau2, noise=0.0, seed=0) -> pd.DataFrame:
    t = np.arange(0.0, 100.0, 0.1)
    y = a1 * np.exp(-t / tau1) + (1 - a1) * np.exp(-t / tau2)
    y = y + np.random.default_rng(seed).normal(0.0, noise, len(t))
    return pd.DataFrame({"t [ps]": t, "TCF": y})


@pytest.mark.parametrize(
    ("a1", "tau1", "tau2", "corrtime"),
    [
        (0.3, 20.0, 2.0, 1.0),
        (0.6, 15.0, 1.0, 1.0),
        (0.5, 5.0, 0.5, 1.0),
        (0.3, 20.0, 2.0, 5.0),  # a poor first guess must not matter
        (0.5, 200.0, 3.0, 1.0),  # a tail longer than the fitting window
    ],
)
def test_lifetime_recovers_a_double_exponential(a1, tau1, tau2, corrtime):
    tau, p = lifetime(double_exponential(a1, tau1, tau2), corrtime)

    assert tau == pytest.approx(a1 * tau1 + (1 - a1) * tau2, rel=1e-3)
    assert p == pytest.approx([a1, tau1, tau2], rel=1e-3)


def test_lifetime_of_a_single_exponential_is_its_time_constant():
    # A single exponential is a degenerate double one (A1 and the two taus
    # are not separately identifiable), but the weighted lifetime is not.
    tau, _ = lifetime(double_exponential(0.5, 4.0, 4.0, noise=0.01), 1.0)
    assert tau == pytest.approx(4.0, rel=0.02)


def test_lifetime_survives_noise():
    tau, _ = lifetime(double_exponential(0.3, 20.0, 2.0, noise=0.01), 1.0)
    assert tau == pytest.approx(0.3 * 20.0 + 0.7 * 2.0, rel=0.05)


# ---------------------------------------------------------------------------
# On the reference trajectory
# ---------------------------------------------------------------------------

# Atom types in solution.data are the LAMMPS type numbers.
OW, HW, NA = "8", "3", "5"


@pytest.fixture(scope="module")
def solution():
    universe = mda.Universe(str(SOLUTION / "solution.data"), str(SOLUTION / "traj.dcd"))
    return universe


def brute_force_tcf(universe, a, b, cutoff, length, start=0, exclude_self=False):
    """One block of the survival function, straight from distance_array."""
    sel_a = universe.select_atoms(f"type {a}")
    sel_b = universe.select_atoms(f"type {b}")
    box = universe.dimensions

    universe.trajectory[start]
    alive = distance_array(sel_a.positions, sel_b.positions, box=box) < cutoff
    if exclude_self:
        np.fill_diagonal(alive, False)
    initial = alive.sum()

    curve = []
    for frame in range(start, start + length):
        universe.trajectory[frame]
        alive &= distance_array(sel_a.positions, sel_b.positions, box=box) < cutoff
        curve.append(alive.sum() / initial)
    return np.array(curve)


@requires_solution
def test_covalent_oh_bonds_of_rigid_water_never_break(solution):
    tcf = bondcorr(
        solution, OW, HW, distance=1.2, dt=100.0, corrlength=200, gaplength=500
    )
    assert tcf["TCF"].to_numpy() == pytest.approx(1.0)


@requires_solution
def test_sodium_hydration_shell_matches_a_brute_force_survival_curve(solution):
    # 3.1 A is the Na-O first minimum in SPC/E. A single block, so the
    # comparison is exact rather than statistical.
    length = 100
    tcf = bondcorr(
        solution,
        NA,
        OW,
        distance=3.1,
        dt=100.0,
        corrlength=length,
        gaplength=5000,
    )
    expected = brute_force_tcf(solution, NA, OW, 3.1, length)

    assert tcf["TCF"].to_numpy() == pytest.approx(expected, abs=1e-6)
    assert 0.0 < tcf["TCF"].iloc[-1] < 1.0  # it does decay, and not to zero


@requires_solution
def test_water_first_shell_decays_monotonically_and_matches_brute_force(solution):
    length = 100
    tcf = bondcorr(
        solution,
        OW,
        HW,
        distance=2.5,
        dt=100.0,
        corrlength=length,
        gaplength=5000,
    )
    expected = brute_force_tcf(solution, OW, HW, 2.5, length)

    curve = tcf["TCF"].to_numpy()
    assert curve == pytest.approx(expected, abs=1e-6)
    assert np.all(np.diff(curve) <= 1e-9)


@requires_solution
def test_water_oxygen_pairs_match_a_brute_force_curve_without_self_pairs(solution):
    length = 100
    tcf = bondcorr(
        solution,
        OW,
        OW,
        distance=3.5,
        dt=100.0,
        corrlength=length,
        gaplength=5000,
    )
    expected = brute_force_tcf(solution, OW, OW, 3.5, length, exclude_self=True)

    assert tcf["TCF"].to_numpy() == pytest.approx(expected, abs=1e-6)
