"""Tests for AtomicSystem.TopologyMixin (types, connections, rule engine)."""

from __future__ import annotations

import pytest

from conftest import make_water_system


# ---------------------------------------------------------------------------
# set_types / set_types_from_elements / set_type2atoms
# ---------------------------------------------------------------------------


def test_set_types_renames_atoms_and_bonds():
    system = make_water_system()
    system.set_types({"Ow": "O", "Hw": "H"})
    assert set(system.atom_types) == {"H", "O"}
    assert set(system.bond_types) == {"H-O"}
    assert set(system.angle_types) == {"H-O-H"}


def test_set_types_renames_every_connection_kind():
    """Dihedrals keep their chain order; the other kinds are made canonical."""
    from cemd import AtomicSystem
    from cemd.core._format import canonical_ff_type

    from conftest import DATA_DIR

    system = AtomicSystem.from_file(DATA_DIR / "caffeine.lt")
    # Prefixes that reverse the alphabetical order of the types, so that
    # the canonical orientation of many connections flips.
    types = sorted(str(t) for t in system.atom_types)
    system.set_types({t: f"{chr(ord('z') - i)}{t}" for i, t in enumerate(types)})

    atom_type = system.atoms["type"].astype(str)
    for kind in ("bond", "angle", "dihedral", "improper"):
        table = getattr(system, kind + "s")
        assert table is not None and len(table), kind
        for _, row in table.iterrows():
            names = tuple(atom_type[row[c]] for c in table.columns if c != "type")
            expected = (
                "-".join(names)
                if kind == "dihedral"
                else canonical_ff_type(names, kind)
            )
            assert row["type"] == expected, kind


def test_set_types_warns_on_unknown_type(water_system):
    with pytest.warns(UserWarning, match="not currently in the system"):
        water_system.set_types({"Zz": "Q"})


def test_set_types_preserves_masses_and_charges():
    system = make_water_system()
    system.set_types({"Ow": "O"})
    assert system.masses["O"] == pytest.approx(15.9994)
    assert system.charges["O"] == pytest.approx(-0.8)


def test_set_types_from_elements_uses_inferred_symbols(water_system):
    water_system.set_types_from_elements()
    assert set(water_system.atom_types) == {"O", "H"}


def test_set_types_from_elements_respects_prevent(water_system):
    water_system.set_types_from_elements(prevent={"Ow": "OW_KEEP"})
    assert "OW_KEEP" in water_system.atom_types
    assert "H" in water_system.atom_types


def test_set_type2atoms_changes_subset(water_system):
    water_system.set_type2atoms([2], "Hnew")
    types = water_system.atoms["type"].tolist()
    assert types[0] == "Ow"
    assert types[1] == "Hnew"
    assert "Hnew" in water_system.masses


def test_set_type2atoms_new_element_type_keeps_its_mass_through_renames(
    water_system,
):
    # Regression test: `set_type2atoms` tested the new type against the
    # `masses` property, which already lists it (falling back on the
    # atomic mass), so nothing was stored; the next `set_types` then gave
    # it 1.008. A C-S-H turned C-A-S-H was written with Al at 1.008.
    water_system.set_type2atoms([1], "Al")
    water_system.set_types({"Hw": "H"})
    assert water_system.masses["Al"] == pytest.approx(26.9815, abs=1e-3)


def test_set_type2atoms_new_type_takes_the_charge_of_the_old_one(water_system):
    water_system.set_type2atoms([2], "Hnew")
    assert water_system.charges["Hnew"] == pytest.approx(0.4)


def test_set_type2atoms_type_mismatch_raises(water_system):
    with pytest.raises(TypeError):
        water_system.set_type2atoms([1], 42)


def test_set_type2atoms_numeric_ok(numeric_type_system):
    numeric_type_system.set_type2atoms([3], 9)
    assert 9 in [int(t) for t in numeric_type_system.atom_types]


# ---------------------------------------------------------------------------
# add_bond / add_angle / add_dihedral / add_improper
# ---------------------------------------------------------------------------


def test_add_bond_creates_canonical_type():
    system = make_water_system()
    system.remove_all_connections()
    system.add_bond([1, 2])
    assert system.num_bonds == 1
    assert system.bonds.iloc[0]["type"] == "Hw-Ow"


def test_add_bond_duplicate_is_noop(capsys):
    system = make_water_system()
    system.remove_all_connections()
    system.add_bond([1, 2])
    system.add_bond([1, 2])
    assert system.num_bonds == 1
    assert "already exists" in capsys.readouterr().out


def test_add_bond_numeric_types_raises(numeric_type_system):
    with pytest.raises(ValueError, match="numeric atom types"):
        numeric_type_system.add_bond([1, 2])


def test_add_angle_and_dihedral_and_improper():
    system = make_water_system()
    system.remove_all_connections()
    system.add_angle([2, 1, 3])
    assert system.angles.iloc[0]["type"] == "Hw-Ow-Hw"

    system.add_dihedral([2, 1, 4, 5])
    assert system.dihedrals.iloc[0]["type"] == "Hw-Ow-Ow-Hw"

    system.add_improper([1, 2, 3, 4])
    assert system.impropers.iloc[0]["type"] == "Ow-Hw-Hw-Ow"


def test_dihedral_type_keeps_chain_order_not_sorted():
    system = make_water_system()
    system.remove_all_connections()
    # atom 2 = Hw, atom1 = Ow, atom4 = Ow, atom6 = Hw
    system.add_dihedral([2, 1, 4, 6])
    assert system.dihedrals.iloc[0]["type"] == "Hw-Ow-Ow-Hw"


# ---------------------------------------------------------------------------
# remove_connection_types / keep_connection_types / remove_all_connections
# ---------------------------------------------------------------------------


def test_remove_connection_types(water_system):
    water_system.remove_connection_types(bond_types=["Hw-Ow"])
    assert water_system.bonds is None
    assert water_system.num_angles == 2  # angles untouched


def test_keep_connection_types_drops_unlisted(water_system):
    water_system.add_angle([2, 1, 3]) if False else None  # no-op, angle already exists
    water_system.keep_connection_types(bond_types=["Hw-Ow"], angle_types=[])
    assert water_system.num_bonds == 4
    assert water_system.angles is None


def test_remove_all_connections(water_system):
    result = water_system.remove_all_connections()
    assert result is water_system
    assert water_system.bonds is None
    assert water_system.angles is None
    assert water_system.dihedrals is None
    assert water_system.impropers is None


# ---------------------------------------------------------------------------
# guess_connections / guess_angles / guess_dihedrals / guess_impropers
# ---------------------------------------------------------------------------


def test_guess_connections_no_bonds_prints_warning(numeric_type_system, capsys):
    result = numeric_type_system.guess_connections()
    assert result is numeric_type_system
    assert "No bonds found" in capsys.readouterr().out


def test_guess_angles_from_bonds():
    system = make_water_system()
    system.remove_connection_types(
        angle_types=list(system.angle_types)
    )  # keep bonds only
    system.guess_angles()
    assert system.num_angles == 2
    assert set(system.angle_types) == {"Hw-Ow-Hw"}


def test_guess_dihedrals_and_impropers_on_chain():
    # Build a 4-atom chain A-B-C-D to exercise dihedral guessing.
    from cemd import AtomicSystem
    import pandas as pd

    atoms = pd.DataFrame(
        {
            "type": ["A", "B", "C", "D"],
            "charge": [0.0] * 4,
            "x": [0.0, 1.0, 2.0, 3.0],
            "y": [0.0, 0.0, 0.0, 0.0],
            "z": [0.0, 0.0, 0.0, 0.0],
        },
        index=[1, 2, 3, 4],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [20, 20, 20, 90, 90, 90],
            "masses": {"A": 1.0, "B": 1.0, "C": 1.0, "D": 1.0},
            "charges": {},
        }
    )
    system.add_bond([1, 2])
    system.add_bond([2, 3])
    system.add_bond([3, 4])

    system.guess_dihedrals()
    assert system.num_dihedrals == 1
    assert system.dihedrals.iloc[0]["type"] == "A-B-C-D"


def test_guess_topo_alias_matches_guess_connections():
    a = make_water_system()
    b = make_water_system()
    a.remove_connection_types(angle_types=list(a.angle_types))
    b.remove_connection_types(angle_types=list(b.angle_types))
    a.guess_connections()
    b.guess_topo()
    assert a.num_angles == b.num_angles


# ---------------------------------------------------------------------------
# set_topology (predefined string style + custom TopologyRule/DihedralRule)
# ---------------------------------------------------------------------------


def test_set_topology_with_custom_rule():
    from cemd.topology import NeighborCriterion, TopologyRule

    system = make_water_system()
    rule = TopologyRule(
        center="type Ow",
        neighbors=[NeighborCriterion("type Hw", 1.2, 2, "H_tagged")],
        new_type="O_tagged",
        bonds=True,
        angles=True,
    )
    system.set_topo(rule)
    assert "O_tagged" in system.atom_types
    assert "H_tagged" in system.atom_types


def test_set_topology_predefined_style_string():
    system = make_water_system()
    system.set_topology("clayff")
    assert set(system.atom_types) == {"Ow", "Hw"}
    assert system.num_bonds == 4
    assert system.num_angles == 2


def test_set_topology_unknown_style_raises():
    system = make_water_system()
    with pytest.raises(ValueError, match="Unknown predefined topology style"):
        system.set_topology("not_a_real_style")


def test_set_topology_wrong_type_raises(water_system):
    with pytest.raises(TypeError):
        water_system.set_topology(42)


def test_set_topology_list_of_rules():
    from cemd.topology import NeighborCriterion, TopologyRule

    system = make_water_system()
    rule = TopologyRule(
        center="type Ow",
        neighbors=[NeighborCriterion("type Hw", 1.2, 2, "H_tagged")],
        new_type="O_tagged",
    )
    system.set_topology([rule])
    assert "O_tagged" in system.atom_types


def test_set_topology_list_rejects_bad_items(water_system):
    with pytest.raises(TypeError):
        water_system.set_topology([42])
