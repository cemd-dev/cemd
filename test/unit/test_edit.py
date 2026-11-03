"""Tests for AtomicSystem.EditMixin (add/remove atoms, box, geometry ops)."""

from __future__ import annotations

import numpy as np
import pytest

from conftest import make_water_system

# ---------------------------------------------------------------------------
# add_atom(s)
# ---------------------------------------------------------------------------


def test_add_atom_appends_row(numeric_type_system):
    numeric_type_system.add_atom("Na", [5.0, 5.0, 5.0], charge=1.0)
    assert numeric_type_system.num_atoms == 4
    last = numeric_type_system.atoms.iloc[-1]
    assert last["type"] == "Na"
    assert last["charge"] == pytest.approx(1.0)
    assert numeric_type_system.masses["Na"] == pytest.approx(22.989769282, rel=1e-3)


def test_add_atoms_batch_with_explicit_mass(numeric_type_system):
    numeric_type_system.add_atoms(
        ["Xx", "Xx"], [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], masses=[5.0, 5.0]
    )
    assert numeric_type_system.num_atoms == 5
    assert numeric_type_system.masses["Xx"] == pytest.approx(5.0)


def test_add_atoms_default_charge_is_zero(numeric_type_system):
    numeric_type_system.add_atoms(["Xx"], [[0.0, 0.0, 0.0]])
    assert numeric_type_system.atoms.iloc[-1]["charge"] == 0.0


def test_add_atom_ids_are_contiguous(empty_atoms_system):
    empty_atoms_system.add_atom("X", [0.0, 0.0, 0.0])
    empty_atoms_system.add_atom("X", [1.0, 0.0, 0.0])
    assert list(empty_atoms_system.atoms.index) == [1, 2]


# ---------------------------------------------------------------------------
# remove_atom(s)
# ---------------------------------------------------------------------------


def test_remove_atom_reindexes_and_updates_topology():
    system = make_water_system()
    system.remove_atom(1)  # the first Ow: drops both of its bonds/angle
    assert system.num_atoms == 5
    assert list(system.atoms.index) == [1, 2, 3, 4, 5]
    assert system.num_bonds == 2
    assert system.num_angles == 1


def test_remove_atoms_drops_unused_types():
    system = make_water_system()
    system.remove_atom(2)
    system.remove_atom(3)  # remove both H of the first water -> "Hw" type gone? no, second water still has Hw
    assert "Ow" in system.atom_types


def test_remove_atoms_clears_unused_type_from_masses():
    system = make_water_system()
    # Remove every Hw atom (indices 2, 3, 5, 6).
    system.remove_atoms([2, 3, 5, 6])
    assert "Hw" not in system.atom_types
    assert "Hw" not in system.masses


def test_remove_all_atoms_sets_topology_to_none():
    system = make_water_system()
    system.remove_atoms([1, 2, 3, 4, 5, 6])
    assert system.num_atoms == 0
    assert system.bonds is None
    assert system.angles is None


def test_remove_atoms_warns_on_missing_index(numeric_type_system):
    with pytest.warns(UserWarning, match="were found"):
        numeric_type_system.remove_atoms([999])
    assert numeric_type_system.num_atoms == 3


def test_remove_atoms_warns_on_empty_system(empty_atoms_system):
    with pytest.warns(UserWarning, match="no atoms"):
        empty_atoms_system.remove_atoms([1])


# ---------------------------------------------------------------------------
# set_box / set_atom_position
# ---------------------------------------------------------------------------


def test_set_box_updates_volume(numeric_type_system):
    numeric_type_system.set_box([20.0, 20.0, 20.0, 90.0, 90.0, 90.0])
    assert numeric_type_system.volume == pytest.approx(8000.0)


def test_set_atom_position(numeric_type_system):
    numeric_type_system.set_atom_position(1, [1.0, 2.0, 3.0])
    row = numeric_type_system.atoms.loc[1]
    np.testing.assert_allclose(row[["x", "y", "z"]].to_numpy(dtype=float), [1.0, 2.0, 3.0])


def test_set_atom_position_wrong_length_raises(numeric_type_system):
    with pytest.raises(ValueError):
        numeric_type_system.set_atom_position(1, [1.0, 2.0])


def test_set_atom_position_missing_index_warns(numeric_type_system, capsys):
    result = numeric_type_system.set_atom_position(999, [1.0, 2.0, 3.0])
    assert result is numeric_type_system
    assert "not found" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# replicate
# ---------------------------------------------------------------------------


def test_replicate_identity_is_noop(numeric_type_system):
    result = numeric_type_system.replicate([1, 1, 1])
    assert result is numeric_type_system
    assert numeric_type_system.num_atoms == 3


def test_replicate_multiplies_atoms_and_box(numeric_type_system):
    n0 = numeric_type_system.num_atoms
    numeric_type_system.replicate([2, 1, 1])
    assert numeric_type_system.num_atoms == 2 * n0
    np.testing.assert_allclose(
        numeric_type_system.box, [20.0, 10.0, 10.0, 90.0, 90.0, 90.0]
    )


def test_replicate_shifts_topology_indices():
    system = make_water_system()
    n0 = system.num_atoms
    b0 = system.num_bonds
    system.replicate([2, 1, 1])
    assert system.num_atoms == 2 * n0
    assert system.num_bonds == 2 * b0
    # Bond references must stay within the new atom index range.
    assert system.bonds[["atom_1", "atom_2"]].to_numpy().max() <= system.num_atoms


@pytest.mark.parametrize("factors", [(2, 1, 1), (3, 2, 1)])
def test_replicate_links_boundary_bonds_to_neighbour_copy(factors):
    # A water molecule split across the x and y faces: every bond and angle must stay
    # short (minimum image) in the replicated box, and so link to the neighbouring copy.
    system = make_water_system()
    system.set_atom_position(1, [29.8, 29.8, 5.0])
    system.set_atom_position(2, [0.75, 29.8, 5.0])
    system.set_atom_position(3, [29.5, 0.7, 5.0])
    system.replicate(list(factors))
    box = np.asarray(system.box[:3])
    pos = system.atoms[["x", "y", "z"]]
    for name, cols in [("bonds", ["atom_1", "atom_2"]), ("angles", ["atom_1", "atom_3"])]:
        df = getattr(system, name)
        d = pos.loc[df[cols[0]]].to_numpy() - pos.loc[df[cols[1]]].to_numpy()
        d -= box * np.rint(d / box)
        assert np.linalg.norm(d, axis=1).max() < 2.0
        assert len(set(map(tuple, df[cols].to_numpy()))) == len(df)


def test_replicate_drops_velocities():
    system = make_water_system()
    import pandas as pd

    system.velocities = pd.DataFrame(
        {"vx": [0.0] * 6, "vy": [0.0] * 6, "vz": [0.0] * 6}, index=range(1, 7)
    )
    system.replicate([2, 1, 1])
    assert system.velocities is None


# ---------------------------------------------------------------------------
# wrap
# ---------------------------------------------------------------------------


def test_wrap_brings_atoms_into_box(numeric_type_system):
    numeric_type_system.set_atom_position(1, [-1.0, 15.0, 5.0])
    numeric_type_system.wrap()
    row = numeric_type_system.atoms.loc[1]
    assert 0.0 <= row["x"] < 10.0
    assert 0.0 <= row["y"] < 10.0
    np.testing.assert_allclose([row["x"], row["y"]], [9.0, 5.0], atol=1e-8)


def test_wrap_empty_system_is_noop(empty_atoms_system):
    empty_atoms_system.wrap()
    assert empty_atoms_system.num_atoms == 0


# ---------------------------------------------------------------------------
# orthogonalize / unskew
# ---------------------------------------------------------------------------


def test_unskew_reduces_large_tilt():
    import pandas as pd

    from cemd import AtomicSystem

    atoms = pd.DataFrame(
        {"type": ["X"], "charge": [0.0], "x": [1.0], "y": [1.0], "z": [1.0]},
        index=[1],
    )
    # A very skewed but valid triclinic box (large xy tilt via a big gamma deviation
    # is awkward to express directly; instead build one via lattice2vectors semantics
    # using an oblique cell with a small acute angle).
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [10.0, 10.0, 10.0, 90.0, 90.0, 20.0],
            "masses": {"X": 1.0},
            "charges": {},
        }
    )
    box_before = system.box.copy()
    system.unskew()
    # unskew should not change the volume of the cell.
    assert system.volume == pytest.approx(
        10.0 * 10.0 * 10.0 * np.sin(np.radians(20.0)), rel=1e-6
    )
    assert not np.allclose(system.box, box_before)


def test_orthogonalize_orthogonal_box_is_noop(numeric_type_system):
    n0 = numeric_type_system.num_atoms
    numeric_type_system.orthogonalize()
    assert numeric_type_system.num_atoms == n0
    np.testing.assert_allclose(
        numeric_type_system.box, [10.0, 10.0, 10.0, 90.0, 90.0, 90.0], atol=1e-6
    )


# ---------------------------------------------------------------------------
# center_on_com (requires MDAnalysis via to_mda())
# ---------------------------------------------------------------------------


def test_center_on_com_moves_com_to_box_center(compact_molecule):
    compact_molecule.center_on_com()
    com = compact_molecule.get_center_of_mass()
    box_center = compact_molecule.box[:3] / 2.0
    np.testing.assert_allclose(com, box_center, atol=1e-6)


def test_center_on_com_with_atom_type_filter(compact_molecule):
    compact_molecule.center_on_com(atom_types=["Ow"])
    # The single oxygen should now sit exactly at the box center.
    ow_pos = compact_molecule.atoms.loc[
        compact_molecule.atoms["type"] == "Ow", ["x", "y", "z"]
    ].to_numpy(dtype=float)[0]
    np.testing.assert_allclose(ow_pos, compact_molecule.box[:3] / 2.0, atol=1e-6)


# ---------------------------------------------------------------------------
# protonate_atom(s) (requires MDAnalysis via to_mda())
# ---------------------------------------------------------------------------


def test_protonate_atom_adds_hydrogen():
    system = make_water_system()
    n0 = system.num_atoms
    # `atom_index` is a real atom id (DataFrame label), like every other
    # index parameter in this class -- atom 1 is the first Ow.
    system.protonate_atom(1, bond_length=1.0)
    assert system.num_atoms == n0 + 1
    last = system.atoms.iloc[-1]
    assert last["type"] == "H"
    assert last["charge"] == pytest.approx(1.0)

    ow_pos = system.atoms.loc[1, ["x", "y", "z"]].to_numpy(dtype=float)
    h_pos = last[["x", "y", "z"]].to_numpy(dtype=float)
    distance = np.linalg.norm(h_pos - ow_pos)
    assert distance == pytest.approx(1.0, rel=1e-6)


def test_protonate_atom_last_atom_id_does_not_crash():
    # Regression test: protonate_atoms used to index the atoms DataFrame
    # positionally (`.iloc[atom_index]`) while every other method in this
    # class takes a real atom id (`.loc`-style); passing the id of the
    # last atom used to raise IndexError instead of protonating it.
    system = make_water_system()
    last_id = int(system.atoms.index[-1])
    system.protonate_atom(last_id, bond_length=1.0)

    new_h = system.atoms.iloc[-1]
    target_pos = system.atoms.loc[last_id, ["x", "y", "z"]].to_numpy(dtype=float)
    distance = np.linalg.norm(
        new_h[["x", "y", "z"]].to_numpy(dtype=float) - target_pos
    )
    assert distance == pytest.approx(1.0, rel=1e-6)


def test_protonate_atoms_batch():
    system = make_water_system()
    n0 = system.num_atoms
    system.protonate_atoms([1, 4])  # both Ow atoms, by real atom id
    assert system.num_atoms == n0 + 2


# ---------------------------------------------------------------------------
# substitute_atoms
# ---------------------------------------------------------------------------


def make_tetrahedral_sheet(n: int = 6, spacing: float = 3.2):
    """A periodic square net of Si, with an O halfway along every Si-Si edge.

    Each Si shares an oxygen with its four neighbours, so Loewenstein's rule
    allows at most half of them to become Al (a checkerboard), and the
    periodic edges are as much bridges as the interior ones.
    """
    import pandas as pd

    from cemd import AtomicSystem

    rows = []
    for i in range(n):
        for j in range(n):
            x, y = i * spacing, j * spacing
            rows.append(("Si", x, y))
            rows.append(("O", x + spacing / 2, y))
            rows.append(("O", x, y + spacing / 2))
    atoms = pd.DataFrame(
        {
            "type": [r[0] for r in rows],
            "charge": [0.0] * len(rows),
            "x": [r[1] for r in rows],
            "y": [r[2] for r in rows],
            "z": [5.0] * len(rows),
        },
        index=range(1, len(rows) + 1),
    )
    box = [n * spacing, n * spacing, 10.0, 90.0, 90.0, 90.0]
    return AtomicSystem(
        {"atoms": atoms, "box": box, "masses": {"Si": 28.0855, "O": 15.9994}, "charges": {}}
    )


def _shared_oxygen_pairs(system, ids, cutoff=2.0):
    """Number of pairs among `ids` bonded to a common oxygen."""
    from MDAnalysis.lib.distances import capped_distance

    oxygens = system.atoms[system.atoms["type"] == "O"][["x", "y", "z"]].to_numpy()
    pos = system.atoms.loc[ids, ["x", "y", "z"]].to_numpy()
    pairs = capped_distance(pos, oxygens, cutoff, box=system.box, return_distances=False)
    per_oxygen = np.bincount(pairs[:, 1], minlength=len(oxygens))
    return int((per_oxygen > 1).sum())


def test_substitute_atoms_obeys_loewenstein_up_to_the_checkerboard():
    system = make_tetrahedral_sheet()
    ids = system.substitute_atoms("Si", "Al", fraction=0.5, seed=0)
    assert len(ids) == 18
    assert (system.atoms.loc[ids, "type"] == "Al").all()
    assert _shared_oxygen_pairs(system, ids) == 0


def test_substitute_atoms_refuses_more_than_loewenstein_allows():
    system = make_tetrahedral_sheet()
    with pytest.raises(ValueError, match="without two sharing"):
        system.substitute_atoms("Si", "Al", count=19, seed=0, max_attempts=20)
    assert (system.atoms["type"] != "Al").all()


def test_substitute_atoms_without_the_rule_can_pair_neighbours():
    system = make_tetrahedral_sheet()
    ids = system.substitute_atoms("Si", "Al", count=30, avoid_adjacent=False, seed=0)
    assert _shared_oxygen_pairs(system, ids) > 0


def test_substitute_atoms_is_reproducible_and_sets_charge_and_mass():
    first = make_tetrahedral_sheet().substitute_atoms("Si", "Al", count=6, seed=3)
    system = make_tetrahedral_sheet()
    second = system.substitute_atoms("Si", "Al", count=6, seed=3, charge=1.575)
    np.testing.assert_array_equal(first, second)
    assert system.atoms.loc[second, "charge"].tolist() == [1.575] * 6
    assert system.masses["Al"] == pytest.approx(26.98, rel=1e-3)


def test_substitute_atoms_stays_within_the_candidates():
    system = make_tetrahedral_sheet()
    silicons = system.atoms.index[system.atoms["type"] == "Si"].to_numpy()
    ids = system.substitute_atoms("Si", "Al", count=4, candidates=silicons[:10], seed=1)
    assert set(ids) <= set(silicons[:10])


def test_substitute_atoms_rejects_candidates_of_another_type():
    system = make_tetrahedral_sheet()
    oxygen = int(system.atoms.index[system.atoms["type"] == "O"][0])
    with pytest.raises(ValueError, match="not of type"):
        system.substitute_atoms("Si", "Al", count=1, candidates=[oxygen])


def test_substitute_atoms_needs_exactly_one_of_count_and_fraction():
    system = make_tetrahedral_sheet()
    with pytest.raises(ValueError, match="exactly one"):
        system.substitute_atoms("Si", "Al")
    with pytest.raises(ValueError, match="exactly one"):
        system.substitute_atoms("Si", "Al", count=2, fraction=0.1)
