"""Tests for AtomicSystem construction, properties and basic queries."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cemd import AtomicSystem

from conftest import DATA_DIR, make_water_system


# ---------------------------------------------------------------------------
# Construction / copy / repr
# ---------------------------------------------------------------------------


def test_construction_requires_mandatory_keys():
    atoms = pd.DataFrame(
        {"type": ["O"], "charge": [0.0], "x": [0.0], "y": [0.0], "z": [0.0]},
        index=[1],
    )
    with pytest.raises(KeyError):
        AtomicSystem({"atoms": atoms, "masses": {"O": 16.0}, "charges": {}})


def test_repr_and_str_without_topology(numeric_type_system):
    text = repr(numeric_type_system)
    assert text == str(numeric_type_system)
    assert "3 atoms" in text
    assert "bonds" not in text


def test_repr_lists_active_interactions(water_system):
    text = repr(water_system)
    assert "6 atoms" in text
    assert "4 bonds" in text
    assert "2 angles" in text


def test_copy_is_independent(water_system):
    clone = water_system.copy()
    assert clone is not water_system
    assert clone.num_atoms == water_system.num_atoms

    clone.set_atom_position(1, [99.0, 99.0, 99.0])
    assert not np.allclose(
        water_system.atoms.loc[1, ["x", "y", "z"]].to_numpy(), [99.0, 99.0, 99.0]
    )

    clone.add_atom("Na", [1.0, 1.0, 1.0])
    assert clone.num_atoms == water_system.num_atoms + 1


# ---------------------------------------------------------------------------
# Properties: atoms/bonds/.../velocities getters and setters
# ---------------------------------------------------------------------------


def test_atoms_setter_normalizes_columns(numeric_type_system):
    df = numeric_type_system.atoms.copy()
    df["extra"] = 0
    numeric_type_system.atoms = df
    assert list(numeric_type_system.atoms.columns) == ["type", "charge", "x", "y", "z"]


def test_bonds_default_none(numeric_type_system):
    assert numeric_type_system.bonds is None
    assert numeric_type_system.num_bonds == 0
    assert numeric_type_system.bond_types == []


def test_bonds_present_after_add_bond(water_system):
    assert water_system.bonds is not None
    assert water_system.num_bonds == 4
    assert set(water_system.bond_types) == {"Hw-Ow"}


def test_velocities_setter(numeric_type_system):
    vel = pd.DataFrame(
        {"vx": [0.1, 0.2, 0.3], "vy": [0.0] * 3, "vz": [0.0] * 3},
        index=[1, 2, 3],
    )
    numeric_type_system.velocities = vel
    assert list(numeric_type_system.velocities.columns) == ["vx", "vy", "vz"]


# ---------------------------------------------------------------------------
# Box / volume
# ---------------------------------------------------------------------------


def test_box_roundtrip_lattice(numeric_type_system):
    np.testing.assert_allclose(
        numeric_type_system.box, [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
    )


def test_volume_orthogonal_box(numeric_type_system):
    assert numeric_type_system.volume == pytest.approx(1000.0)


def test_volume_triclinic_box():
    atoms = pd.DataFrame(
        {"type": ["X"], "charge": [0.0], "x": [0.0], "y": [0.0], "z": [0.0]},
        index=[1],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [10.0, 10.0, 10.0, 90.0, 90.0, 60.0],
            "masses": {"X": 1.0},
            "charges": {},
        }
    )
    # V = a*b*c*sqrt(1 - cos^2(alpha) - cos^2(beta) - cos^2(gamma) + 2*cos*cos*cos)
    expected = 10.0 * 10.0 * 10.0 * np.sin(np.radians(60.0))
    assert system.volume == pytest.approx(expected, rel=1e-6)


# ---------------------------------------------------------------------------
# masses / charges / elements
# ---------------------------------------------------------------------------


def test_masses_and_charges_are_read_only_mappings(numeric_type_system):
    masses = numeric_type_system.masses
    with pytest.raises(TypeError):
        masses[1] = 999.0

    charges = numeric_type_system.charges
    with pytest.raises(TypeError):
        charges[1] = 999.0


def test_masses_fallback_to_masses_dict():
    atoms = pd.DataFrame(
        {"type": ["Ca"], "charge": [0.0], "x": [0.0], "y": [0.0], "z": [0.0]},
        index=[1],
    )
    system = AtomicSystem(
        {"atoms": atoms, "box": [10, 10, 10, 90, 90, 90], "masses": {}, "charges": {}}
    )
    assert system.masses["Ca"] == pytest.approx(40.0784, rel=1e-3)


def test_elements_inferred_from_mass(water_system):
    elements = water_system.elements
    assert elements["Ow"] == "O"
    assert elements["Hw"] == "H"


# ---------------------------------------------------------------------------
# Type lists and counts
# ---------------------------------------------------------------------------


def test_atom_types_sorted_as_strings(numeric_type_system):
    assert numeric_type_system.atom_types == ["1", "2"]


def test_num_and_type_counts(water_system):
    assert water_system.num_atoms == 6
    assert water_system.num_atom_types == 2
    assert water_system.num_bonds == 4
    assert water_system.num_bond_types == 1
    assert water_system.num_angles == 2
    assert water_system.num_angle_types == 1
    assert water_system.num_dihedrals == 0
    assert water_system.num_dihedral_types == 0
    assert water_system.num_impropers == 0
    assert water_system.num_improper_types == 0


def test_total_charge_and_mass(water_system):
    assert water_system.total_charge == pytest.approx(0.0, abs=1e-9)
    expected_mass = 2 * 15.9994 + 4 * 1.007947
    assert water_system.total_mass == pytest.approx(expected_mass)


def test_density(numeric_type_system):
    # 2 atoms of mass 12.011 + 1 atom of mass 15.999, box volume 1000 A^3
    expected = (2 * 12.011 + 15.999) / 6.02214076e23 / 1000.0 / 1e-24
    assert numeric_type_system.density == pytest.approx(expected, rel=1e-5)


def test_get_count(water_system):
    assert water_system.get_count("Ow") == 2
    assert water_system.get_count("Hw") == 4
    assert water_system.get_count("Na") == 0


# ---------------------------------------------------------------------------
# summary() / get_center_of_mass()
# ---------------------------------------------------------------------------


def test_summary_runs_and_prints(capsys, water_system):
    water_system.summary()
    captured = capsys.readouterr()
    assert "Box" in captured.out
    assert "Atoms" in captured.out
    assert "Bonds" in captured.out
    assert "Total charge" in captured.out


def test_get_center_of_mass_symmetric_pair():
    atoms = pd.DataFrame(
        {
            "type": ["X", "X"],
            "charge": [0.0, 0.0],
            "x": [0.0, 10.0],
            "y": [0.0, 0.0],
            "z": [0.0, 0.0],
        },
        index=[1, 2],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [20, 20, 20, 90, 90, 90],
            "masses": {"X": 1.0},
            "charges": {},
        }
    )
    com = system.get_center_of_mass()
    np.testing.assert_allclose(com, [5.0, 0.0, 0.0])


def test_get_center_of_mass_empty_system(empty_atoms_system):
    com = empty_atoms_system.get_center_of_mass()
    np.testing.assert_allclose(com, [0.0, 0.0, 0.0])


# ---------------------------------------------------------------------------
# ff_keys / ff_params active filtering
# ---------------------------------------------------------------------------


def test_ff_keys_and_params_are_filtered_to_active_types(water_system):
    water_system.set_bond_ff_keys({"Hw-Ow": "flexible_spc.oh"})
    # Remove all bonds -> the bond type is no longer "active".
    water_system.remove_all_connections()

    assert water_system.ff_keys.bond == {}
    # The underlying (unfiltered) storage still has it.
    assert "Hw-Ow" in water_system._ff_keys.bond


def test_make_water_system_helper_is_independent():
    a = make_water_system()
    b = make_water_system()
    a.add_atom("Na", [1.0, 1.0, 1.0])
    assert a.num_atoms == b.num_atoms + 1


def test_replace_internals_updates_every_box_representation():
    # Regression test: `_replace_internals` assigned `_box_lmp` directly and
    # left `_box` / `_box_vectors` untouched, so `.box` and `.volume` kept
    # reporting the previous cell after add_layer(), add_droplet() or
    # set_topology().
    system = AtomicSystem.from_file(DATA_DIR / "calcite.cif")
    other = system.copy()
    other.set_box([50.0, 50.0, 50.0, 90.0, 90.0, 90.0])

    system._replace_internals(other)

    np.testing.assert_allclose(system.box, [50.0, 50.0, 50.0, 90.0, 90.0, 90.0])
    assert system.volume == pytest.approx(50.0**3)


def test_total_charge_sums_per_atom_charges():
    # Regression test: total_charge rebuilt the sum from the *per-type*
    # `charges` mapping (count x representative charge). A force field that
    # assigns a charge per atom rather than per type -- ReaxFF's EEM
    # charges, ATB/GROMOS partial charges -- has no single per-type value,
    # so a strictly neutral system was reported as heavily charged.
    atoms = pd.DataFrame(
        {
            "type": ["O", "O", "O", "O"],
            "charge": [-0.7, -0.5, 0.5, 0.7],
            "x": [0.0, 1.0, 2.0, 3.0],
            "y": [0.0, 0.0, 0.0, 0.0],
            "z": [0.0, 0.0, 0.0, 0.0],
        },
        index=range(1, 5),
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [10.0, 10.0, 10.0, 90.0, 90.0, 90.0],
            "masses": {"O": 15.999},
            "charges": {},
        }
    )

    assert system.total_charge == pytest.approx(0.0)


def test_total_charge_still_follows_set_charges():
    # Per-type force fields must keep working: set_charges writes the
    # per-atom column too, so the sum tracks it.
    system = make_water_system()
    system.set_charges({"Ow": -0.8, "Hw": 0.4})
    n_water = (system.atoms["type"] == "Ow").sum()

    assert system.total_charge == pytest.approx(0.0, abs=1e-9)
    assert system.charges["Ow"] == pytest.approx(-0.8)
    assert n_water > 0


def test_elements_skips_a_mass_matching_nothing():
    # Regression test: the element was taken as the nearest mass in the
    # table, however far away. The table stops at barium, so platinum
    # (195.08) came back as barium (137.33) -- a 58 amu error, reported
    # with no warning at all.
    atoms = pd.DataFrame(
        {
            "type": ["O", "Pt"],
            "charge": [0.0, 0.0],
            "x": [0.0, 3.0],
            "y": [0.0, 0.0],
            "z": [0.0, 0.0],
        },
        index=[1, 2],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
            "masses": {"O": 15.9994, "Pt": 195.084},
            "charges": {},
        }
    )

    with pytest.warns(UserWarning, match="No element matches the mass"):
        elements = dict(system.elements)

    assert elements == {"O": "O"}
    assert "Pt" not in elements


def test_elements_still_resolves_ordinary_types():
    # The guard must not cost anything on a normal system: every element in
    # the bundled force fields matches its mass to better than 0.05 amu.
    system = make_water_system()
    assert dict(system.elements) == {"Ow": "O", "Hw": "H"}
