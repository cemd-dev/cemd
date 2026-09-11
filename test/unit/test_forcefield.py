"""Tests for AtomicSystem.ForceFieldMixin (masses/charges/keys/params)."""

from __future__ import annotations

import pandas as pd
import pytest

from cemd.forcefield.models import LJParams

from cemd import AtomicSystem

from conftest import DATA_DIR, make_water_system


# ---------------------------------------------------------------------------
# set_masses / set_charges
# ---------------------------------------------------------------------------


def test_set_masses_updates_value(water_system):
    water_system.set_masses({"Ow": 18.0})
    assert water_system.masses["Ow"] == pytest.approx(18.0)


def test_set_masses_rejects_non_positive(water_system):
    with pytest.raises(ValueError, match="strictly positive"):
        water_system.set_masses({"Ow": 0.0})
    with pytest.raises(ValueError, match="strictly positive"):
        water_system.set_masses({"Ow": -1.0})


def test_set_masses_rejects_non_dict(water_system):
    with pytest.raises(TypeError):
        water_system.set_masses([1.0, 2.0])


def test_set_masses_warns_on_unknown_type(water_system):
    with pytest.warns(UserWarning, match="missing from current system"):
        water_system.set_masses({"Zz": 1.0})


def test_set_masses_numeric_key_matches_numeric_system(numeric_type_system):
    # Regression test: dict keys must be normalized the same way atom_types
    # are (string), or an int key silently fails to match.
    numeric_type_system.set_masses({1: 99.0})
    assert numeric_type_system.masses["1"] == pytest.approx(99.0)


def test_set_charges_updates_atoms_dataframe_and_property(water_system):
    water_system.set_charges({"Ow": -1.0})
    assert water_system.charges["Ow"] == pytest.approx(-1.0)
    assert (water_system.atoms.loc[water_system.atoms["type"] == "Ow", "charge"] == -1.0).all()


def test_total_charge_reflects_charge_updates(water_system):
    water_system.set_charges({"Ow": -1.0, "Hw": 0.5})
    expected = 2 * (-1.0) + 4 * 0.5
    assert water_system.total_charge == pytest.approx(expected)


# ---------------------------------------------------------------------------
# set_pair_params / apply_pair_mixing_rules
# ---------------------------------------------------------------------------


def test_set_pair_params_self_interaction(water_system):
    water_system.set_pair_params("Ow", params=LJParams(epsilon=0.15, sigma=3.16))
    key = water_system._normalize_binary_key("Ow", "Ow")
    assert key in water_system._ff_params.pair


def test_set_pair_params_unknown_type_raises(water_system):
    with pytest.raises(ValueError, match="does not exist"):
        water_system.set_pair_params("Zz", params=LJParams(epsilon=0.1, sigma=3.0))


def test_apply_pair_mixing_rules_arithmetic(water_system):
    water_system.set_pair_params("Ow", params=LJParams(epsilon=0.15535, sigma=3.166))
    water_system.set_pair_params("Hw", params=LJParams(epsilon=0.0, sigma=0.0))
    water_system.apply_pair_mixing_rules(rule="arithmetic")

    cross_key = water_system._normalize_binary_key("Ow", "Hw")
    cross = water_system._ff_params.pair[cross_key]
    assert cross.epsilon == pytest.approx((0.15535 * 0.0) ** 0.5)
    assert cross.sigma == pytest.approx((3.166 + 0.0) / 2.0)


def test_apply_pair_mixing_rules_geometric(water_system):
    water_system.set_pair_params("Ow", params=LJParams(epsilon=0.2, sigma=3.0))
    water_system.set_pair_params("Hw", params=LJParams(epsilon=0.05, sigma=1.0))
    water_system.apply_pair_mixing_rules(rule="geometric")

    cross_key = water_system._normalize_binary_key("Ow", "Hw")
    cross = water_system._ff_params.pair[cross_key]
    assert cross.sigma == pytest.approx((3.0 * 1.0) ** 0.5)


def test_apply_pair_mixing_rules_missing_self_params_raises(water_system):
    water_system.set_pair_params("Ow", params=LJParams(epsilon=0.1, sigma=3.0))
    with pytest.raises(ValueError, match="Missing parameters"):
        water_system.apply_pair_mixing_rules()


def test_apply_pair_mixing_rules_unknown_rule_raises(water_system):
    water_system.set_pair_params("Ow", params=LJParams(epsilon=0.1, sigma=3.0))
    water_system.set_pair_params("Hw", params=LJParams(epsilon=0.1, sigma=1.0))
    with pytest.raises(ValueError, match="Unknown mixing rule"):
        water_system.apply_pair_mixing_rules(rule="bogus")


def test_apply_pair_mixing_rules_does_not_overwrite_by_default(water_system):
    water_system.set_pair_params("Ow", params=LJParams(epsilon=0.1, sigma=3.0))
    water_system.set_pair_params("Hw", params=LJParams(epsilon=0.1, sigma=1.0))
    sentinel = LJParams(epsilon=999.0, sigma=999.0)
    water_system.set_pair_params("Ow", "Hw", params=sentinel)

    water_system.apply_pair_mixing_rules(overwrite=False)
    cross_key = water_system._normalize_binary_key("Ow", "Hw")
    assert water_system._ff_params.pair[cross_key] is sentinel


# ---------------------------------------------------------------------------
# set_bond/angle/dihedral/improper_params
#
# Regression coverage: these four public methods all delegated to a helper,
# `_apply_topology_param`, that did not exist anywhere in the codebase --
# every call raised AttributeError.
# ---------------------------------------------------------------------------


def test_set_bond_params_stores_canonicalized(water_system):
    water_system.set_bond_params("Ow-Hw", [450.0, 1.0])  # reverse of canonical "Hw-Ow"
    assert water_system._ff_params.bond["Hw-Ow"] == [450.0, 1.0]


def test_set_angle_params_stores_canonicalized(water_system):
    water_system.set_angle_params("Hw-Ow-Hw", [55.0, 104.5])
    assert water_system._ff_params.angle["Hw-Ow-Hw"] == [55.0, 104.5]


def test_set_dihedral_params_keeps_chain_order():
    system = make_water_system()
    system.add_dihedral([2, 1, 4, 5])
    system.set_dihedral_params("Hw-Ow-Ow-Hw", [1.0, 1, 0])
    assert system._ff_params.dihedral["Hw-Ow-Ow-Hw"] == [1.0, 1, 0]


def test_set_improper_params_warns_when_type_absent(water_system):
    with pytest.warns(UserWarning, match="not present"):
        water_system.set_improper_params("A-B-C-D", [1.0])
    assert water_system._ff_params.improper["A-B-C-D"] == [1.0]


# ---------------------------------------------------------------------------
# set_ff_keys / set_*_ff_keys
# ---------------------------------------------------------------------------


def test_set_atom_ff_keys_dict(water_system):
    water_system.set_atom_ff_keys({"Ow": "spc.ospc", "Hw": "spc.hspc"})
    assert water_system.ff_keys.atom == {"Ow": "spc.ospc", "Hw": "spc.hspc"}


def test_set_atom_ff_keys_sequence_must_match_length(water_system):
    with pytest.raises(ValueError, match="does not match"):
        water_system.set_atom_ff_keys(["only_one"])


def test_set_bond_ff_keys_canonicalizes(water_system):
    water_system.set_bond_ff_keys({"Ow-Hw": "spc.hspc-ospc"})
    assert water_system.ff_keys.bond == {"Hw-Ow": "spc.hspc-ospc"}


def test_set_ff_keys_overwrite_false_keeps_existing(water_system):
    water_system.set_ff_keys(atom={"Ow": "first"})
    water_system.set_ff_keys(atom={"Ow": "second"}, overwrite=False)
    assert water_system.ff_keys.atom["Ow"] == "first"


# ---------------------------------------------------------------------------
# set_ff_from_database (real database lookup, real SPC water model)
# ---------------------------------------------------------------------------


def test_set_ff_from_database_assigns_masses_charges_and_pair_params(water_system):
    water_system.set_ff_from_database(
        atom_assignments={"Ow": "spc.ospc", "Hw": "spc.hspc"},
        bond_assignments={"Hw-Ow": "spc.hspc-ospc"},
    )

    assert water_system.masses["Ow"] == pytest.approx(15.9994, rel=1e-2)
    assert water_system.charges["Ow"] == pytest.approx(-0.82)
    assert water_system.charges["Hw"] == pytest.approx(0.41)

    pair_key = water_system._normalize_binary_key("Ow", "Ow")
    assert pair_key in water_system._ff_params.pair
    assert water_system._ff_params.pair[pair_key].sigma == pytest.approx(3.166)

    assert "Hw-Ow" in water_system._ff_params.bond


def test_set_ff_from_database_warns_on_unresolvable_atom_abbreviation(water_system):
    with pytest.warns(UserWarning, match="not found in database"):
        water_system.set_ff_from_database(
            atom_assignments={"Ow": "not_a_real_forcefield_abbreviation"}
        )


def test_set_ff_from_database_keeps_charges_a_force_field_does_not_define():
    # Regression test: `_update_masses_and_charges` applied the database's
    # atom charge unconditionally, while the mass was guarded. Force fields
    # that carry their charges per atom rather than per type (GROMOS, the
    # CHARMM Interface parameter files) report no per-type charge, and that
    # was written back as 0.0 -- so a molecule solvated from an ATB topology
    # lost all its partial charges the moment its parameters were resolved.
    caffeine = AtomicSystem.from_file(DATA_DIR / "caffeine.lt")
    before = dict(caffeine.charges)
    assert any(abs(q) > 1e-9 for q in before.values())

    caffeine.set_ff_from_database()

    for atom_type, charge in before.items():
        assert caffeine.charges[atom_type] == pytest.approx(charge)


def test_set_ff_from_database_applies_the_charges_a_force_field_defines():
    # The other half of the contract checked just above: when the force
    # field does define a per-type charge, it wins over whatever the system
    # was carrying. Only a force field that defines none (GROMOS, the
    # CHARMM Interface files) leaves the existing charges alone.
    #
    # +/-9.99 is a sentinel, not a physical charge: no force field would
    # ever assign it, so the assertions below fail loudly if ClayFF's
    # values are not written over it.
    atoms = pd.DataFrame(
        {
            "type": ["ao", "ob"],
            "charge": [9.99, -9.99],
            "x": [0.0, 2.0],
            "y": [0.0, 0.0],
            "z": [0.0, 0.0],
        },
        index=[1, 2],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [10.0, 10.0, 10.0, 90.0, 90.0, 90.0],
            "masses": {"ao": 26.98, "ob": 16.0},
            "charges": {"ao": 9.99, "ob": -9.99},
        }
    )

    system.set_ff_keys(atom={"ao": "clayff.ao", "ob": "clayff.ob"})
    system.set_ff_from_database()

    assert system.charges["ao"] == pytest.approx(1.575)
    assert system.charges["ob"] == pytest.approx(-1.05)
    # The per-atom column follows, so the total is the force field's.
    assert system.total_charge == pytest.approx(1.575 - 1.05)
