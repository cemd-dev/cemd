"""Tests for cemd.build.Splitter."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cemd import AtomicSystem
from cemd.build import Splitter, SolutionBuilder

from conftest import requires_packmol


def _make_chain_system() -> AtomicSystem:
    """Four atoms in a chain along z: bonds 1-2 and 3-4 stay on one side
    of coordinate=5, bond 2-3 crosses it."""
    atoms = pd.DataFrame(
        {
            "type": ["X", "X", "X", "X"],
            "charge": [0.0] * 4,
            "x": [0.0] * 4,
            "y": [0.0] * 4,
            "z": [1.0, 4.0, 6.0, 9.0],
        },
        index=range(1, 5),
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [10.0, 10.0, 10.0, 90.0, 90.0, 90.0],
            "masses": {"X": 1.0},
            "charges": {},
        }
    )
    system.add_bond([1, 2])
    system.add_bond([2, 3])
    system.add_bond([3, 4])
    return system


def test_split_opens_gap_and_grows_box():
    system = _make_chain_system()
    result = Splitter(system, coordinate=5.0, axis="z", gap_size=10.0).split()

    assert result is system
    assert result.box[2] == pytest.approx(20.0)

    # Atoms below the cut are untouched; atoms above are shifted by gap_size.
    assert result.atoms.loc[1, "z"] == pytest.approx(1.0)
    assert result.atoms.loc[2, "z"] == pytest.approx(4.0)
    assert result.atoms.loc[3, "z"] == pytest.approx(16.0)
    assert result.atoms.loc[4, "z"] == pytest.approx(19.0)


def test_split_removes_bond_crossing_the_gap():
    system = _make_chain_system()
    result = Splitter(system, coordinate=5.0, axis="z", gap_size=10.0).split()

    bond_pairs = {
        frozenset((int(row.atom_1), int(row.atom_2)))
        for row in result.bonds.itertuples()
    }
    assert frozenset((2, 3)) not in bond_pairs
    assert frozenset((1, 2)) in bond_pairs
    assert frozenset((3, 4)) in bond_pairs


def test_bond_through_an_unrelated_periodic_boundary_survives_the_cut():
    """A bond that wraps through the box's own periodic boundary -- common
    in any packed/replicated system -- has nothing to do with the cut
    plane elsewhere in the box. Measured as a raw Cartesian difference it
    looks artificially long even before the split, and moving only one of
    its two atoms then falsely reads as 'stretched by the cut'."""
    atoms = pd.DataFrame(
        {
            "type": ["X", "X"],
            "charge": [0.0, 0.0],
            "x": [0.0, 0.0],
            "y": [0.0, 0.0],
            "z": [19.6, 0.6],  # true (minimum-image) separation: 1.0 A
        },
        index=[1, 2],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
            "masses": {"X": 1.0},
            "charges": {},
        }
    )
    system.add_bond([1, 2])

    # Cut at z=10 moves atom 1 (z=19.6 >= 10) but not atom 2 (z=0.6 < 10);
    # neither atom is anywhere near the cut plane.
    result = Splitter(system, coordinate=10.0, axis="z", gap_size=10.0).split()

    assert result.bonds is not None
    assert len(result.bonds) == 1


def test_axis_accepts_string_or_int():
    a = Splitter(_make_chain_system(), coordinate=5.0, axis="z", gap_size=10.0)
    b = Splitter(_make_chain_system(), coordinate=5.0, axis=2, gap_size=10.0)
    assert a.axis == b.axis == 2


def test_with_coordinate_sets_cut_position():
    splitter = Splitter(_make_chain_system(), coordinate=0.0, axis="z", gap_size=10.0)
    splitter.with_coordinate(5.0)
    assert splitter.coordinate == 5.0


def test_build_is_an_alias_for_split():
    system = _make_chain_system()
    result = Splitter(system, coordinate=5.0, axis="z", gap_size=10.0).build()
    assert result.box[2] == pytest.approx(20.0)


# ---------------------------------------------------------------------------
# add_solution / split (requires packmol)
# ---------------------------------------------------------------------------


@requires_packmol
def test_split_with_solution_fills_the_gap():
    atoms = pd.DataFrame(
        {
            "type": ["Na", "Cl"],
            "charge": [1.0, -1.0],
            "x": [0.0, 6.0],
            "y": [0.0, 6.0],
            "z": [1.0, 1.0],
        },
        index=[1, 2],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [12.0, 12.0, 10.0, 90.0, 90.0, 90.0],
            "masses": {"Na": 22.989769282, "Cl": 35.4532},
            "charges": {},
        }
    )
    blueprint = SolutionBuilder.from_water(density=1.0)

    result = (
        Splitter(system, coordinate=5.0, axis="z", gap_size=20.0)
        .add_solution(blueprint, padding=2.0, vacuum=0.0)
        .split()
    )

    assert result.num_atoms > 2
    assert "Ow" in result.atom_types
    assert "Na" in result.atom_types and "Cl" in result.atom_types


@requires_packmol
def test_solution_is_placed_inside_the_pore_not_stacked_on_top():
    # Regression test: `_insert_solution` centred the liquid in the gap and
    # then handed it to `_merge_structure`, the surface-builder helper that
    # lays a film *on top* of a slab -- which translated it above the
    # topmost atom and stretched the cell a second time. The pore came out
    # empty, the fluid sat outside the solid, and the box was far taller
    # than `original + gap_size`.
    atoms = pd.DataFrame(
        {
            "type": ["Na", "Cl"],
            "charge": [1.0, -1.0],
            "x": [0.0, 6.0],
            "y": [0.0, 6.0],
            "z": [1.0, 1.0],
        },
        index=[1, 2],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [12.0, 12.0, 10.0, 90.0, 90.0, 90.0],
            "masses": {"Na": 22.989769282, "Cl": 35.4532},
            "charges": {},
        }
    )
    blueprint = SolutionBuilder.from_water(density=1.0)

    result = (
        Splitter(system, coordinate=5.0, axis="z", gap_size=20.0)
        .add_solution(blueprint, padding=2.0, vacuum=0.0)
        .split()
    )

    # The box grows by exactly one gap, not by the gap plus a liquid slab.
    assert result.box[2] == pytest.approx(30.0)

    # And the water sits within the gap, [5, 25] in the new coordinates.
    water_z = result.atoms.loc[result.atoms["type"] == "Ow", "z"]
    assert water_z.min() >= 5.0
    assert water_z.max() <= 25.0


@requires_packmol
def test_split_with_solution_keeps_the_angles_of_the_solution():
    # The merge used to copy only atoms/bonds back onto the system, so the
    # solution's angles (and dihedrals, ff keys and charges) were dropped.
    system = _make_chain_system()
    blueprint = SolutionBuilder.from_water(density=1.0)

    result = (
        Splitter(system, coordinate=5.0, axis="z", gap_size=20.0)
        .add_solution(blueprint, padding=2.0, vacuum=0.0)
        .split()
    )

    n_water = (result.atoms["type"] == "Ow").sum()
    assert n_water > 0
    assert result.angles is not None
    assert len(result.angles) == n_water


@requires_packmol
def test_split_with_solution_gap_too_small_raises():
    system = _make_chain_system()
    blueprint = SolutionBuilder.from_water(density=1.0)

    splitter = Splitter(
        system, coordinate=5.0, axis="z", gap_size=3.0
    ).add_solution(blueprint, padding=2.0, vacuum=0.0)

    with pytest.raises(ValueError, match="too small"):
        splitter.split()


# ---------------------------------------------------------------------------
# Severed-contact detection (distance-based, no explicit bonds needed)
# ---------------------------------------------------------------------------


def _make_layered_system() -> AtomicSystem:
    """Two silicate 'layers' (vertical Si-O pairs at z=0 and z=12) with a
    water interlayer at z=7. No explicit bonds anywhere -- the Si-O
    framework carries none under ClayFF/CSHFF conventions.
    """
    types: list[str] = []
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []

    for layer_z in (0.0, 12.0):
        for k in range(3):
            types += ["Si", "O"]
            xs += [k * 4.0, k * 4.0]
            ys += [0.0, 0.0]
            zs += [layer_z, layer_z + 1.6]

    for k in range(2):
        types += ["Ow", "Hw", "Hw"]
        xs += [k * 4.0, k * 4.0 + 0.96, k * 4.0 - 0.24]
        ys += [0.0, 0.0, 0.93]
        zs += [7.0, 7.0, 7.0]

    atoms = pd.DataFrame(
        {"type": types, "charge": [0.0] * len(types), "x": xs, "y": ys, "z": zs},
        index=range(1, len(types) + 1),
    )
    return AtomicSystem(
        {
            "atoms": atoms,
            "box": [20.0, 20.0, 24.0, 90.0, 90.0, 90.0],
            "masses": {"Si": 28.085, "O": 15.999, "Ow": 15.999, "Hw": 1.008},
            "charges": {},
        }
    )


def test_detects_contacts_that_carry_no_explicit_bond():
    system = _make_layered_system()
    assert system.bonds is None  # nothing explicit to rely on

    splitter = Splitter(system, coordinate=1.0, axis="z", gap_size=15.0)
    assert splitter.count_broken_bonds() == 3


def test_find_broken_bonds_returns_the_severed_pairs():
    system = _make_layered_system()
    splitter = Splitter(system, coordinate=1.0, axis="z", gap_size=15.0)

    broken = splitter.find_broken_bonds()
    assert len(broken) == 3
    for id_1, id_2 in broken:
        elements = {
            system.elements[str(system.atoms.loc[i, "type"])] for i in (id_1, id_2)
        }
        assert elements == {"Si", "O"}


def test_cut_through_the_interlayer_breaks_nothing():
    system = _make_layered_system()
    splitter = Splitter(system, coordinate=5.0, axis="z", gap_size=15.0)
    assert splitter.count_broken_bonds() == 0


def test_scan_broken_bonds_separates_layers_from_interlayer():
    system = _make_layered_system()
    splitter = Splitter(system, coordinate=0.0, axis="z", gap_size=15.0)

    scan = splitter.scan_broken_bonds(step=1.0)

    assert list(scan.columns) == ["coordinate", "n_broken", "n_ionic"]
    assert scan["n_broken"].max() == 3
    # The best cut sits in the interlayer, away from both silicate layers.
    best = scan.loc[scan["n_broken"] == 0, "coordinate"]
    assert ((best > 2.0) & (best < 12.0)).any()


def test_scan_broken_bonds_rejects_non_positive_step():
    splitter = Splitter(_make_layered_system(), coordinate=1.0, axis="z")
    with pytest.raises(ValueError, match="step must be positive"):
        splitter.scan_broken_bonds(step=0.0)


def test_type_keys_take_priority_over_element_keys():
    """A type-keyed dict can exclude contacts that element keys would catch."""
    atoms = pd.DataFrame(
        {
            "type": ["Si", "O", "Ow", "Hw"],
            "charge": [0.0] * 4,
            "x": [0.0, 0.0, 8.0, 8.0],
            "y": [0.0] * 4,
            "z": [0.4, 2.0, 0.6, 1.56],
        },
        index=[1, 2, 3, 4],
    )
    payload = {
        "atoms": atoms,
        "box": [20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
        "masses": {"Si": 28.085, "O": 15.999, "Ow": 15.999, "Hw": 1.008},
        "charges": {},
    }

    # Default table keys on elements, so the water O-H is caught too.
    default = Splitter(AtomicSystem(payload), coordinate=1.0, axis="z")
    assert sorted(default.find_broken_bonds()) == [(1, 2), (3, 4)]

    # Keying on types lets the water float free, as it should.
    framework_only = Splitter(
        AtomicSystem(payload),
        coordinate=1.0,
        axis="z",
        bonds_dict={("Si", "O"): 1.8},
    )
    assert framework_only.find_broken_bonds() == [(1, 2)]


def test_contacts_through_the_periodic_boundary_are_not_broken():
    # Si and O are 1.6 A apart *through* the z boundary; the split shifts
    # the box and the moving fragment together, so the contact survives.
    atoms = pd.DataFrame(
        {
            "type": ["O", "Si"],
            "charge": [0.0, 0.0],
            "x": [0.0, 0.0],
            "y": [0.0, 0.0],
            "z": [0.8, 19.2],
        },
        index=[1, 2],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
            "masses": {"Si": 28.085, "O": 15.999},
            "charges": {},
        }
    )
    splitter = Splitter(system, coordinate=10.0, axis="z", gap_size=10.0)
    assert splitter.count_broken_bonds() == 0


# ---------------------------------------------------------------------------
# Automatic repair
# ---------------------------------------------------------------------------


def _make_broken_pair_system() -> AtomicSystem:
    """A single vertical Si-O contact, severed by a cut at z=1.0."""
    atoms = pd.DataFrame(
        {
            "type": ["Si", "O"],
            "charge": [0.0, 0.0],
            "x": [0.0, 0.0],
            "y": [0.0, 0.0],
            "z": [0.4, 2.0],
        },
        index=[1, 2],
    )
    return AtomicSystem(
        {
            "atoms": atoms,
            "box": [20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
            "masses": {"Si": 28.085, "O": 15.999},
            "charges": {},
        }
    )


def test_split_without_repair_leaves_dangling_atoms():
    system = _make_broken_pair_system()
    Splitter(system, coordinate=1.0, axis="z", gap_size=10.0).split()
    assert system.num_atoms == 2


def test_repair_caps_both_sides_of_a_severed_contact():
    system = _make_broken_pair_system()
    result = Splitter(system, coordinate=1.0, axis="z", gap_size=10.0).split(
        repair=True
    )

    # Si keeps a hydroxyl (new O + H); the displaced O is capped with an H.
    assert result.num_atoms == 5
    assert result.get_count("H") == 2
    assert result.get_count("O") == 2
    assert result.get_count("Si") == 1


def test_repair_places_caps_along_the_broken_direction():
    system = _make_broken_pair_system()
    result = Splitter(system, coordinate=1.0, axis="z", gap_size=10.0).split(
        repair=True, oh_length=1.0
    )

    zs = sorted(result.atoms["z"].tolist())
    # Si stays at 0.4 with its restored O at 2.0 and H at 3.0; the original
    # O moved to 12.0 and is capped by an H at 11.0.
    assert zs == pytest.approx([0.4, 2.0, 3.0, 11.0, 12.0])


def test_repair_respects_a_restricted_bonds_dict():
    system = _make_broken_pair_system()
    # No cutoff covers Si-O here, so nothing is detected and nothing capped.
    result = Splitter(system, coordinate=1.0, axis="z", gap_size=10.0).split(
        repair=True, bonds_dict={("Ow", "Hw"): 1.1}
    )
    assert result.num_atoms == 2


def test_build_forwards_repair_to_split():
    system = _make_broken_pair_system()
    result = Splitter(system, coordinate=1.0, axis="z", gap_size=10.0).build(
        repair=True
    )
    assert result.num_atoms == 5


def test_repair_reports_what_it_did():
    system = _make_broken_pair_system()
    splitter = Splitter(system, coordinate=1.0, axis="z", gap_size=10.0)
    splitter.split(repair=True)

    assert splitter.repair_report == {"broken": 1, "capped": 3, "skipped": 0}


def test_repair_resets_a_specifically_typed_dangling_oxygen_to_generic():
    """An oxygen that already carries a specific force-field type (e.g.
    'Osi', a non-bridging Si-OH oxygen) loses that type's meaning once the
    cut removes its Si neighbor. `set_topology()` only ever revisits atoms
    matching a literal 'type O' selection, so leaving the stale specific
    type in place would make the atom permanently unreachable for
    re-typing."""
    atoms = pd.DataFrame(
        {
            "type": ["Si", "Osi"],
            "charge": [1.5, -0.8],
            "x": [0.0, 0.0],
            "y": [0.0, 0.0],
            "z": [0.4, 2.0],
        },
        index=[1, 2],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
            "masses": {"Si": 28.085, "Osi": 15.999},
            "charges": {},
        }
    )
    splitter = Splitter(
        system, coordinate=1.0, axis="z", gap_size=10.0,
        bonds_dict={("Si", "Osi"): 1.8},
    )
    result = splitter.split(repair=True)

    # The dangling former-Osi (atom 2) is reset to generic "O"; the Si-side
    # restored oxygen (atom 3) was already generic on creation.
    assert result.atoms.loc[2, "type"] == "O"
    assert result.atoms.loc[3, "type"] == "O"


def test_repair_does_not_stack_atoms_on_a_shared_bridging_oxygen():
    """A bridging O that loses both its cations must not have each of them
    restore a copy of it at the very same site."""
    atoms = pd.DataFrame(
        {
            "type": ["Si", "Si", "O"],
            "charge": [0.0] * 3,
            "x": [-1.0, 1.0, 0.0],
            "y": [0.0, 0.0, 0.0],
            "z": [0.6, 0.6, 1.8],
        },
        index=[1, 2, 3],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
            "masses": {"Si": 28.085, "O": 15.999},
            "charges": {},
        }
    )

    splitter = Splitter(system, coordinate=1.2, axis="z", gap_size=10.0)
    assert splitter.count_broken_bonds() == 2

    result = splitter.split(repair=True)

    assert splitter.repair_report["broken"] == 2
    assert splitter.repair_report["skipped"] >= 1

    # No two atoms may end up sitting on top of each other.
    positions = result.atoms[["x", "y", "z"]].to_numpy(dtype=float)
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            assert np.linalg.norm(positions[i] - positions[j]) > 0.5


# ---------------------------------------------------------------------------
# keep_molecules: cut along the molecules, not along a plane
# ---------------------------------------------------------------------------


def test_keep_molecules_breaks_no_bond_where_a_plane_would():
    system = _make_layered_system()
    plane = Splitter(_make_layered_system(), coordinate=1.0, axis="z", gap_size=15.0)
    assert plane.count_broken_bonds() == 3

    splitter = Splitter(
        system, coordinate=1.0, axis="z", gap_size=15.0, keep_molecules=True
    )
    assert splitter.count_broken_bonds() == 0

    splitter.split()
    assert splitter.split_report["bonds_broken"] == []
    # Each Si-O pair moved as a whole: still 1.6 A apart
    z = system.atoms["z"].to_numpy()
    np.testing.assert_allclose(z[1:12:2] - z[0:12:2], 1.6, atol=1e-5)


def test_keep_molecules_keeps_explicit_bonds_beyond_the_cutoffs():
    # O-H at 1.15 A: longer than the default H-O cutoff, but bonded
    atoms = pd.DataFrame(
        {
            "type": ["Ow", "Hw"],
            "charge": [0.0, 0.0],
            "x": [0.0, 0.0],
            "y": [0.0, 0.0],
            "z": [4.5, 5.65],
        },
        index=[1, 2],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [10.0, 10.0, 10.0, 90.0, 90.0, 90.0],
            "masses": {"Ow": 15.999, "Hw": 1.008},
            "charges": {},
        }
    )
    system.add_bond([1, 2])

    splitter = Splitter(system, coordinate=5.0, axis="z", gap_size=10.0)
    assert splitter.count_broken_bonds() == 1

    splitter = Splitter(
        system, coordinate=5.0, axis="z", gap_size=10.0, keep_molecules=True
    )
    splitter.split()
    assert system.num_bonds == 1
    assert splitter.split_report["bonds_broken"] == []
    assert splitter.split_report["explicit_bonds_removed"] == 0


def test_keep_molecules_keeps_a_molecule_wrapped_through_the_boundary():
    # Si at the top of the cell, its O through the boundary at the bottom,
    # and a cut right next to the O: the pair must still move together.
    atoms = pd.DataFrame(
        {
            "type": ["O", "Si"],
            "charge": [0.0, 0.0],
            "x": [0.0, 0.0],
            "y": [0.0, 0.0],
            "z": [0.8, 19.2],
        },
        index=[1, 2],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
            "masses": {"Si": 28.085, "O": 15.999},
            "charges": {},
        }
    )
    splitter = Splitter(
        system, coordinate=0.5, axis="z", gap_size=10.0, keep_molecules=True
    )
    assert splitter.count_broken_bonds() == 0
    splitter.split()
    assert splitter.split_report["bonds_broken"] == []


def test_keep_molecules_puts_a_cation_with_most_of_its_shell():
    # Ca just above the cut, three O below it and one above
    atoms = pd.DataFrame(
        {
            "type": ["Ca", "O", "O", "O", "O"],
            "charge": [0.0] * 5,
            "x": [5.0, 7.4, 2.6, 5.0, 5.0],
            "y": [5.0, 5.0, 5.0, 7.4, 5.0],
            "z": [5.2, 4.4, 4.4, 4.4, 7.6],
        },
        index=range(1, 6),
    )
    payload = {
        "atoms": atoms,
        "box": [20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
        "masses": {"Ca": 40.078, "O": 15.999},
        "charges": {},
    }

    plane = Splitter(AtomicSystem(payload), coordinate=5.0, axis="z")
    assert plane.count_broken_ionic() == 3

    splitter = Splitter(
        AtomicSystem(payload), coordinate=5.0, axis="z", keep_molecules=True
    )
    assert not splitter.side_mask()[0]  # Ca stays below, with its 3 O
    assert splitter.count_broken_ionic() == 1


def test_coordinate_none_cuts_through_an_interlayer():
    system = _make_layered_system()
    splitter = Splitter(system, axis="z", gap_size=15.0, keep_molecules=True)
    splitter.split()

    # Either interlayer will do (one runs through the z boundary), as long
    # as the plane itself crosses no layer.
    plane = Splitter(_make_layered_system(), coordinate=splitter.coordinate, axis="z")
    assert plane.count_broken_bonds() == 0
    assert splitter.split_report["bonds_broken"] == []


def test_periodic_covalent_network_cannot_be_kept_whole():
    # A Si-O chain running through the cell along z, closed on itself
    zs = [0.0, 1.6, 3.2, 4.8]
    atoms = pd.DataFrame(
        {
            "type": ["Si", "O", "Si", "O"],
            "charge": [0.0] * 4,
            "x": [0.0] * 4,
            "y": [0.0] * 4,
            "z": zs,
        },
        index=range(1, 5),
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [10.0, 10.0, 6.4, 90.0, 90.0, 90.0],
            "masses": {"Si": 28.085, "O": 15.999},
            "charges": {},
        }
    )
    splitter = Splitter(
        system, coordinate=2.0, axis="z", gap_size=10.0, keep_molecules=True
    )
    with pytest.raises(ValueError, match="periodic"):
        splitter.side_mask()


@requires_packmol
def test_solution_stays_clear_of_units_sticking_into_the_pore():
    system = _make_layered_system()
    # Cut between Si (z=12) and its O (z=13.6): the O stays with its Si,
    # 1.6 A above the cut plane
    splitter = Splitter(
        system, coordinate=13.0, axis="z", gap_size=20.0, keep_molecules=True
    ).add_solution(SolutionBuilder.from_water(density=1.0), padding=2.0)
    result = splitter.split()

    n_solid = 12 + 6
    z = result.atoms["z"].to_numpy()
    liquid = z[n_solid:]
    assert liquid.min() >= 13.6 + 2.0 - 1e-6
    assert splitter.split_report["bonds_broken"] == []


def test_stretched_si_o_bond_counts_as_a_bond():
    # Regression test: at 1.8 A, the default Si-O cutoff missed the
    # stretched bonds of ReaxFF-relaxed C-S-H (up to ~2.0 A)
    atoms = pd.DataFrame(
        {
            "type": ["Si", "O"],
            "charge": [0.0, 0.0],
            "x": [0.0, 0.0],
            "y": [0.0, 0.0],
            "z": [4.0, 5.9],
        },
        index=[1, 2],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
            "masses": {"Si": 28.085, "O": 15.999},
            "charges": {},
        }
    )
    assert Splitter(system, coordinate=5.0, axis="z").count_broken_bonds() == 1
