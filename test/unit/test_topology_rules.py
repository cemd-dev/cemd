"""Tests for cemd.topology: predefined rule sets, DihedralRule, and the
low-level rule-application engine in cemd.topology._apply.
"""

from __future__ import annotations

import pandas as pd
import pytest

from cemd import AtomicSystem
from cemd.topology import CLAYFF_RULES, CSHFF_RULES, DihedralRule
from cemd.topology.rules import NeighborCriterion, TopologyRule

from conftest import make_water_system


def _make_silica_system() -> AtomicSystem:
    """Si1-O-Si2 bridge, Si1-O-H silanol, plus one free water molecule."""
    atoms = pd.DataFrame(
        {
            "type": ["Si", "Si", "O", "O", "H", "O", "H", "H"],
            "charge": [0.0] * 8,
            "x": [0.0, 3.2, 1.6, -1.6, -2.5, 10.0, 10.96, 9.76],
            "y": [0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.93],
            "z": [0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0],
        },
        index=range(1, 9),
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
            "masses": {"Si": 28.085, "O": 15.999, "H": 1.008},
            "charges": {},
        }
    )
    system.add_bond([1, 3])
    system.add_bond([2, 3])
    system.add_bond([1, 4])
    system.add_bond([4, 5])
    system.add_bond([6, 7])
    system.add_bond([6, 8])
    return system


# ---------------------------------------------------------------------------
# CLAYFF_RULES / CSHFF_RULES applied through AtomicSystem.set_topology
# ---------------------------------------------------------------------------


def test_clayff_rules_classify_bridging_and_silanol_oxygens():
    system = _make_silica_system()
    system.set_topology(CLAYFF_RULES)

    types = set(system.atom_types)
    assert "Ob" in types  # bridging Si-O-Si
    assert "Osih" in types  # Si-OH silanol oxygen
    assert "Hsi" in types  # silanol hydrogen
    assert "Ow" in types  # free water oxygen
    assert "Hw" in types  # free water hydrogen


def test_cshff_rules_add_calcium_water_coordination():
    atoms = pd.DataFrame(
        {
            "type": ["Ca", "O", "H", "H"],
            "charge": [0.0] * 4,
            "x": [0.0, 2.5, 3.46, 2.26],
            "y": [0.0, 0.0, 0.0, 0.93],
            "z": [0.0, 0.0, 0.0, 0.0],
        },
        index=range(1, 5),
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
            "masses": {"Ca": 40.078, "O": 15.999, "H": 1.008},
            "charges": {},
        }
    )
    system.add_bond([2, 3])
    system.add_bond([2, 4])

    system.set_topology(CSHFF_RULES)

    assert "Cw" in system.atom_types  # Ca within range of a water oxygen
    assert set(system.atom_types) == {"Cw", "Ow", "Hw"}


def test_clayff_string_and_list_forms_are_equivalent():
    a = _make_silica_system()
    b = _make_silica_system()

    a.set_topology("clayff")
    b.set_topology(CLAYFF_RULES)

    assert sorted(a.atom_types) == sorted(b.atom_types)


# ---------------------------------------------------------------------------
# DihedralRule
#
# Regression coverage: `apply_single_dihedral_rule_to_universe` used to
# iterate `for rule in dihedral_rules` even though every caller (both the
# single-rule and the list-of-rules branch in `set_topology`) passed a lone
# DihedralRule, not a list -- `TypeError: 'DihedralRule' object is not
# iterable` on every use.
# ---------------------------------------------------------------------------


def test_dihedral_rule_single_instance():
    system = make_water_system()
    rule = DihedralRule(
        i="type Hw", j="type Ow", k="type Ow", l_="type Hw", cutoffs=[3.0, 15.0, 3.0]
    )
    result = system.set_topology(rule)

    assert result is system
    assert system.num_dihedrals > 0


def test_dihedral_rule_inside_a_list():
    system = make_water_system()
    rule = DihedralRule(
        i="type Hw", j="type Ow", k="type Ow", l_="type Hw", cutoffs=[3.0, 15.0, 3.0]
    )
    system.set_topology([rule])
    assert system.num_dihedrals > 0


def test_dihedral_rule_mixed_with_topology_rule():
    system = make_water_system()
    topo_rule = TopologyRule(
        center="type Ow",
        neighbors=[NeighborCriterion("type Hw", 1.2, 2, "H_tagged")],
        new_type="O_tagged",
    )
    dih_rule = DihedralRule(
        i="type Hw", j="type Ow", k="type Ow", l_="type Hw", cutoffs=[3.0, 15.0, 3.0]
    )

    system.set_topology([topo_rule, dih_rule])

    # The TopologyRule renamed the water oxygens before the DihedralRule ran,
    # so the dihedral selections ("type Hw"/"type Ow") no longer match
    # anything -- this exercises the no-match early-return path rather than
    # asserting a particular dihedral count.
    assert "O_tagged" in system.atom_types


def test_dihedral_rule_no_matches_is_a_noop():
    system = make_water_system()
    rule = DihedralRule(i="type Xx", j="type Ow", k="type Ow", l_="type Hw")
    result = system.set_topology(rule)

    assert result is system
    assert system.num_dihedrals == 0


# ---------------------------------------------------------------------------
# Neighbor lists vs one selection per center
# ---------------------------------------------------------------------------


def _random_mix(n_atoms: int, seed: int) -> AtomicSystem:
    """A dense random mix, where centers share neighbors and renames chain."""
    import numpy as np

    rng = np.random.default_rng(seed)
    types = rng.choice(["Si", "Al", "O", "O", "O", "H", "H", "C", "S", "Ca"], n_atoms)
    xyz = rng.uniform(0, 14.0, (n_atoms, 3))
    masses = {
        "Si": 28.09,
        "Al": 26.98,
        "O": 16.0,
        "H": 1.008,
        "C": 12.01,
        "S": 32.06,
        "Ca": 40.08,
    }
    return AtomicSystem(
        {
            "atoms": pd.DataFrame(
                {
                    "type": types,
                    "charge": 0.0,
                    "x": xyz[:, 0],
                    "y": xyz[:, 1],
                    "z": xyz[:, 2],
                },
                index=range(1, n_atoms + 1),
            ),
            "box": [14.0, 14.0, 14.0, 90.0, 90.0, 90.0],
            "masses": masses,
            "charges": dict.fromkeys(masses, 0.0),
        }
    )


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_neighbor_lists_match_one_selection_per_center(seed):
    """The fast path must type and connect exactly as the per-center one.

    The rules rename atoms as they go and later centers see the new names,
    so this compares the whole sequence, not each rule on a fresh system.
    The extra rules cover inexact counts, impropers, a selection that is
    not a plain type list, and a neighbor renamed by a rule.
    """
    from cemd.topology import _apply

    extra = [
        TopologyRule(
            "type Ca",
            NeighborCriterion("type O Ob Osi Oh Ow", 3.0, 2, exact_match=False),
            new_type="Cx",
            bonds=True,
            impropers=True,
        ),
        TopologyRule(
            "type C",
            NeighborCriterion("name * and (type O or type Oc)", 1.6, 3),
            angles=True,
            impropers=True,
        ),
        TopologyRule(
            "type S",
            NeighborCriterion("type O*", 1.8, 1, "Oq", exact_match=False),
            bonds=True,
        ),
    ]

    system = _random_mix(800, seed)
    reference, fast = system.to_mda(), system.to_mda()
    for rule in CSHFF_RULES + extra:
        _apply._apply_rule_per_center(reference, rule)
        _apply._apply_rule_with_neighbor_lists(fast, rule)

    def state(universe):
        def indices(attr):
            if not hasattr(universe, attr):
                return []
            return sorted(map(tuple, getattr(universe, attr).indices.tolist()))

        return (
            list(universe.atoms.types),
            indices("bonds"),
            indices("angles"),
            indices("impropers"),
        )

    expected = state(reference)
    assert expected[1], "the mix should produce bonds"
    assert state(fast) == expected
