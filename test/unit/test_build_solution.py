"""Tests for cemd.build.SolutionBuilder."""

from __future__ import annotations

import pytest

from cemd.build import SolutionBuilder

from conftest import requires_packmol

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_negative_density_raises():
    with pytest.raises(ValueError, match="Density must be positive"):
        SolutionBuilder(density=-1.0)


def test_negative_molarity_raises():
    with pytest.raises(ValueError, match="Molarity must be >= 0"):
        SolutionBuilder(molarities={"NaCl": -0.1})


def test_non_numeric_molarity_raises():
    with pytest.raises(TypeError, match="must be a number"):
        SolutionBuilder(molarities={"NaCl": "0.1"})


def test_negative_count_raises():
    with pytest.raises(ValueError, match="Count must be >= 0"):
        SolutionBuilder(counts={"Na": -5})


def test_non_int_count_raises():
    with pytest.raises(TypeError, match="must be an integer"):
        SolutionBuilder(counts={"Na": 5.0})


def test_from_water_is_empty_composition():
    blueprint = SolutionBuilder.from_water(density=1.0)
    assert blueprint.molarities == {}
    assert blueprint.counts == {}
    assert "pure water" in repr(blueprint)


# ---------------------------------------------------------------------------
# to_counts / get_solute_mass / get_water_count
# ---------------------------------------------------------------------------


def test_to_counts_converts_molarity_to_integer():
    blueprint = SolutionBuilder(density=1.0, molarities={"Na": 1.0})
    volume = 30.0**3  # Å^3
    counts = blueprint.to_counts(volume)
    # 1 M in a 30^3 A^3 box: n = M * N_A * V(L)
    assert counts["Na"] == pytest.approx(16, abs=1)


def test_to_counts_rejects_species_in_both_counts_and_molarities():
    blueprint = SolutionBuilder(
        density=1.0, molarities={"Na": 1.0}, counts={"Na": 10}
    )
    with pytest.raises(ValueError, match="both"):
        blueprint.to_counts(30.0**3)


def test_to_counts_combines_counts_and_molarities():
    blueprint = SolutionBuilder(
        density=1.0, molarities={"Na": 1.0}, counts={"Cl": 5}
    )
    counts = blueprint.to_counts(30.0**3)
    assert counts["Cl"] == 5
    assert "Na" in counts


def test_get_solute_mass_unknown_species_raises():
    blueprint = SolutionBuilder(density=1.0, counts={"NotAnElementXYZ": 5})
    with pytest.raises(ValueError, match="not found in mass database"):
        blueprint.get_solute_mass(30.0**3)


def test_get_solute_mass_knows_the_bundled_polyatomic_species():
    # "HO", "CO3" and "SO4" ship with CEMD and are resolved by name when
    # packing, but `get_solute_mass` only consulted the element table, so
    # they could be packed and not weighed -- and had to be handed back in
    # through `structures=` to be usable at all.
    blueprint = SolutionBuilder(density=1.0, counts={"HO": 4})
    # 4 x (O + H)
    assert blueprint.get_solute_mass(30.0**3) == pytest.approx(4 * 17.007, abs=0.01)


def test_get_water_count_positive_for_reasonable_density():
    blueprint = SolutionBuilder.from_water(density=1.0)
    n = blueprint.get_water_count(30.0**3)
    assert n > 0
    # Sanity check against the well-known ~33.3 waters/nm^3 for d=1 g/cm^3.
    expected = round(1.0 * 30.0**3 * 6.02214076e23 * 1e-24 / 18.015)
    assert n == pytest.approx(expected, rel=0.02)


def test_get_water_count_raises_when_density_too_low_for_solutes():
    # A huge count of a heavy solute in a tiny box can't fit under the
    # target density once solute mass alone exceeds the target total mass.
    blueprint = SolutionBuilder(density=0.01, counts={"Xe": 10000})
    with pytest.raises(ValueError, match="too low"):
        blueprint.get_water_count(30.0**3)


# ---------------------------------------------------------------------------
# build() / build_hemisphere() (require packmol)
# ---------------------------------------------------------------------------


@requires_packmol
def test_build_pure_water_box():
    blueprint = SolutionBuilder.from_water(density=1.0)
    system = blueprint.build(box=[18.0, 18.0, 18.0])

    assert system.num_atoms > 0
    assert set(system.atom_types) == {"Ow", "Hw"}
    np_water = system.get_count("Ow")
    assert np_water == pytest.approx(blueprint.get_water_count(18.0**3), abs=1)

    # Regression test: every packmol-built system used to come back with
    # every atom's charge silently overwritten to 1.0 (a PDB round-trip
    # bug -- see test_io.py::test_pdb_roundtrip_does_not_corrupt_charges).
    # Water is neutral and its H/O charges must not both be 1.0.
    assert system.total_charge == pytest.approx(0.0, abs=1e-6)
    assert system.charges["Ow"] != pytest.approx(1.0)
    assert system.charges["Hw"] != pytest.approx(1.0)


@requires_packmol
def test_build_with_ion_counts():
    blueprint = SolutionBuilder(density=1.0, counts={"Na": 4, "Cl": 4})
    system = blueprint.build(box=[18.0, 18.0, 18.0])

    assert system.get_count("Na") == 4
    assert system.get_count("Cl") == 4
    assert "Ow" in system.atom_types
    assert system.total_charge == pytest.approx(0.0, abs=1e-6)


@requires_packmol
def test_build_accepts_a_bundled_polyatomic_species_by_name():
    # Regression test: packmol only reads PDB, but `ho.sdf` was handed to
    # it unconverted ("Could not read any atom from file", exit 171), so a
    # blueprint naming "HO" failed even though the structure ships with
    # CEMD. Its O-H bonds must survive the build too.
    blueprint = SolutionBuilder(density=1.0, counts={"HO": 4})
    system = blueprint.build(box=[18.0, 18.0, 18.0])

    assert system.get_count("O") == 4
    assert system.get_count("H") == 4
    bond_counts = system.bonds.groupby("type").size().to_dict()
    assert bond_counts.get("H-O") == 4


@requires_packmol
def test_build_hemisphere_uses_consistent_water_naming():
    # Regression test: build_hemisphere() used to skip the H1/H2/O1 ->
    # Hw/Ow rename that build() applies, leaving water under the raw
    # template type names and breaking downstream topology-rule matching.
    blueprint = SolutionBuilder.from_water(density=1.0)
    system = blueprint.build_hemisphere(radius=8.0)

    assert system.num_atoms > 0
    assert set(system.atom_types) == {"Ow", "Hw"}


@requires_packmol
def test_build_keeps_bonds_of_templates_that_have_no_ff_keys():
    # Regression test: a custom structure read from a file that carries no
    # force-field assignment (an SDF, say) still carries real connectivity.
    # `_rebuild_topology_from_templates` used to drop every such bond --
    # hydroxide ions came out as two unbonded atoms.
    from cemd import AtomicSystem
    from cemd.build._structures import STRUCTURES_DIR

    hydroxide = AtomicSystem.from_file(STRUCTURES_DIR / "ho.sdf")
    assert hydroxide.num_bonds == 1

    blueprint = SolutionBuilder(
        density=1.0, counts={"Na": 4, "HO": 4}, structures={"HO": hydroxide}
    )
    system = blueprint.build(box=[18.0, 18.0, 18.0])

    hydroxide_bonds = sum(1 for t in system.bonds["type"] if t == "H-O")
    assert hydroxide_bonds == 4


@requires_packmol
def test_build_hemisphere_rejects_invalid_axis():
    blueprint = SolutionBuilder.from_water(density=1.0)
    with pytest.raises(ValueError, match="Invalid axis"):
        blueprint.build_hemisphere(radius=8.0, axis="q")


@requires_packmol
def test_build_hemisphere_rejects_non_positive_radius():
    blueprint = SolutionBuilder.from_water(density=1.0)
    with pytest.raises(ValueError, match="radius must be positive"):
        blueprint.build_hemisphere(radius=0.0)


# ---------------------------------------------------------------------------
# Cancellation
# ---------------------------------------------------------------------------


def test_waiting_on_packmol_stops_when_cancelled():
    import subprocess
    import threading
    import time

    from cemd.build import PackmolCancelledError
    from cemd.build._packmol import _wait_for_packmol

    process = subprocess.Popen(["sleep", "30"])
    cancel = threading.Event()
    threading.Timer(0.2, cancel.set).start()

    start = time.perf_counter()
    with pytest.raises(PackmolCancelledError):
        _wait_for_packmol(process, cancel)
    assert time.perf_counter() - start < 5
    # Killed and reaped, not left running.
    assert process.poll() is not None


@requires_packmol
def test_a_build_inside_cancel_packmol_on_can_be_cancelled():
    import threading
    import time

    from cemd.build import PackmolCancelledError, cancel_packmol_on

    # ~13 000 waters: several seconds of Packmol, far more than the delay.
    blueprint = SolutionBuilder.from_water(density=1.0)
    cancel = threading.Event()
    threading.Timer(0.5, cancel.set).start()

    start = time.perf_counter()
    with cancel_packmol_on(cancel), pytest.raises(PackmolCancelledError):
        blueprint.build(box=[75.0, 75.0, 75.0])
    assert time.perf_counter() - start < 10


def test_cancellation_is_not_a_runtime_error():
    # Code that retries a failed packing catches RuntimeError; it must not
    # take a cancellation for a failure and start Packmol again.
    from cemd.build import PackmolCancelledError

    assert not issubclass(PackmolCancelledError, RuntimeError)
