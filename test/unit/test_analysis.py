"""Tests for cemd.analysis.

The assertions are against values that are known analytically rather than
against a previous run: a uniformly random gas has g(r) = 1 everywhere by
definition, and a box of N atoms has a density of N / V whichever way you
project it. That is what catches a normalisation being wrong -- a density
map can be off by a factor of forty and still look perfectly sensible,
because a map is read for its contrast.
"""

from __future__ import annotations

import numpy as np
import pytest

mda = pytest.importorskip("MDAnalysis")
from MDAnalysis.coordinates.memory import MemoryReader  # noqa: E402

from cemd.analysis import compute_rdf, density_map, density_profile  # noqa: E402


BOX = 20.0
N_ATOMS = 8000
EXPECTED_DENSITY = N_ATOMS / BOX**3 * 1000  # atoms/nm^3


def uniform_universe(n_frames: int = 1, seed: int = 0) -> "mda.Universe":
    """A box of uniformly distributed atoms, in memory."""
    rng = np.random.default_rng(seed)
    universe = mda.Universe.empty(
        N_ATOMS,
        n_residues=N_ATOMS,
        atom_resindex=np.arange(N_ATOMS),
        residue_segindex=np.zeros(N_ATOMS, dtype=int),
        trajectory=True,
    )
    universe.add_TopologyAttr("type", ["Ow"] * N_ATOMS)
    universe.add_TopologyAttr("mass", [15.999] * N_ATOMS)
    universe.load_new(
        rng.uniform(0, BOX, (n_frames, N_ATOMS, 3)), format=MemoryReader
    )
    # Per frame: setting universe.dimensions only touches the current one,
    # and several analysis functions index the box without checking.
    for ts in universe.trajectory:
        ts.dimensions = [BOX, BOX, BOX, 90.0, 90.0, 90.0]
    return universe


# ---------------------------------------------------------------------------
# Density
# ---------------------------------------------------------------------------


def test_density_profile_recovers_a_uniform_density():
    universe = uniform_universe()
    profile = density_profile(universe, ["Ow"], axis="z", bin_size=0.5, end=1)

    assert profile["Ow"].mean() == pytest.approx(EXPECTED_DENSITY, rel=0.01)


def test_density_map_recovers_a_uniform_density():
    # Regression test: the map normalised by a 1D slab volume,
    # bin_size * box[a] * box[b], where a bin of a map is a column of
    # section bin_size^2. Densities came out low by
    # box[a] * box[b] / (bin_size * depth) -- 40x here, and a factor that
    # moves with the box and the binning rather than a constant offset.
    universe = uniform_universe()
    density = density_map(
        universe,
        "Ow",
        interface_coordinate=0.0,
        eps=BOX,  # the whole box
        axis="z",
        bin_size=0.5,
        end=1,
    )

    assert density.values.mean() == pytest.approx(EXPECTED_DENSITY, rel=0.02)


def test_density_profile_and_map_agree():
    # Two projections of the same atoms cannot disagree on the mean.
    universe = uniform_universe()
    profile = density_profile(universe, ["Ow"], axis="z", bin_size=0.5, end=1)
    density = density_map(
        universe, "Ow", interface_coordinate=0.0, eps=BOX, axis="z",
        bin_size=0.5, end=1,
    )

    assert density.values.mean() == pytest.approx(profile["Ow"].mean(), rel=0.02)


# ---------------------------------------------------------------------------
# RDF
# ---------------------------------------------------------------------------


def test_rdf_of_a_uniform_gas_is_flat_at_one():
    # No structure, so g(r) = 1 at every distance. This is what pins the
    # normalisation: a wrong shell volume or particle density tilts the
    # curve away from 1 without changing its shape.
    universe = uniform_universe()
    frame, _ = compute_rdf(universe, "Ow", "Ow", cutoff=8.0, dr=0.2)

    # The first bins hold very few pairs, so they are noisy; judge from 2 A.
    tail = frame.loc[frame.index > 2.0, "g_r"]
    assert tail.mean() == pytest.approx(1.0, abs=0.05)
    assert tail.std() < 0.05


def test_rdf_returns_the_pair_density_it_normalised_by():
    universe = uniform_universe()
    _, rho = compute_rdf(universe, "Ow", "Ow", cutoff=8.0, dr=0.2)

    assert rho == pytest.approx(N_ATOMS / BOX**3, rel=1e-6)


def test_rdf_coordination_number_matches_the_uniform_expectation():
    # For g(r) = 1, n(r) is just the number of atoms in a sphere:
    # (4/3) pi r^3 rho.
    universe = uniform_universe()
    frame, rho = compute_rdf(universe, "Ow", "Ow", cutoff=8.0, dr=0.2)

    r = 5.0
    expected = 4.0 / 3.0 * np.pi * r**3 * rho
    measured = frame.loc[frame.index <= r, "n_r"].iloc[-1]

    assert measured == pytest.approx(expected, rel=0.05)


# ---------------------------------------------------------------------------
# MSD and diffusion
# ---------------------------------------------------------------------------


def random_walk_universe(n_atoms=400, n_frames=600, step=0.30, seed=3):
    """Brownian walkers whose diffusion coefficient is known exactly.

    Each atom takes an independent Gaussian step of variance `step`^2 per
    axis per frame, so MSD_x(n) = n * step^2 and the Einstein slope is
    step^2 / dt. Nothing here is fitted to a previous run: the answer comes
    from the way the trajectory was built.
    """
    rng = np.random.default_rng(seed)
    box = 200.0  # far larger than any displacement, so no wrapping
    steps = rng.normal(0.0, step, (n_frames, n_atoms, 3))
    steps[0] = 0.0
    positions = np.cumsum(steps, axis=0) + box / 2

    universe = mda.Universe.empty(
        n_atoms,
        n_residues=n_atoms,
        atom_resindex=np.arange(n_atoms),
        residue_segindex=np.zeros(n_atoms, dtype=int),
        trajectory=True,
    )
    universe.add_TopologyAttr("type", ["Ow"] * n_atoms)
    universe.add_TopologyAttr("mass", [15.999] * n_atoms)
    universe.load_new(positions, format=MemoryReader)
    for ts in universe.trajectory:
        ts.dimensions = [box, box, box, 90.0, 90.0, 90.0]
    return universe, step


def test_msd_of_a_random_walk_grows_linearly():
    from cemd.analysis import msd

    universe, step = random_walk_universe()
    frame = msd(universe, "Ow", dt=1000.0)  # 1 ps between frames

    time = frame.index.to_numpy()
    isotropic = frame[["xx", "yy", "zz"]].mean(axis=1).to_numpy()

    # MSD_x(n frames) = n * step^2, and one frame is one picosecond here.
    expected = time * step**2
    half = len(time) // 2
    assert isotropic[1:half] == pytest.approx(expected[1:half], rel=0.10)


def test_diffusion_coefficient_matches_the_walk_it_was_given():
    from cemd.analysis import diffusion_coefficient, msd

    universe, step = random_walk_universe()
    frame = msd(universe, "Ow", dt=1000.0)
    result = diffusion_coefficient(frame, start=frame.index[len(frame) // 10])

    # D = step^2 / (2 dt), in A^2/ps, converted to m^2/s.
    expected = step**2 / 2.0 * 1e-8
    assert result.loc["DC (m2/s)", "3d"] == pytest.approx(expected, rel=0.10)


def test_diffusion_of_an_isotropic_walk_has_no_preferred_axis():
    from cemd.analysis import diffusion_coefficient, msd

    universe, _ = random_walk_universe()
    frame = msd(universe, "Ow", dt=1000.0)
    result = diffusion_coefficient(frame, start=frame.index[len(frame) // 10])

    diagonal = result.loc["DC (m2/s)", ["xx", "yy", "zz"]].to_numpy(dtype=float)
    assert diagonal.std() / diagonal.mean() < 0.10

    # Off-diagonal terms describe correlated motion between axes; there is
    # none here by construction.
    cross = result.loc["DC (m2/s)", ["xy", "xz", "yz"]].to_numpy(dtype=float)
    assert np.abs(cross).max() < 0.1 * diagonal.mean()
