"""Tests for AtomicSystem.IOMixin (readers, writers, converters)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from cemd import AtomicSystem

from conftest import DATA_DIR, make_water_system


# ---------------------------------------------------------------------------
# from_file: one check per supported extension
# ---------------------------------------------------------------------------


def test_from_file_missing_raises():
    with pytest.raises(FileNotFoundError):
        AtomicSystem.from_file(DATA_DIR / "does_not_exist.data")


def test_from_file_unsupported_extension(tmp_path):
    bogus = tmp_path / "system.xyz"
    bogus.write_text("not a real format")
    with pytest.raises(ValueError, match="Unsupported file format"):
        AtomicSystem.from_file(bogus)


def test_from_file_pdb():
    system = AtomicSystem.from_file(DATA_DIR / "h2o.pdb")
    assert system.num_atoms == 3
    assert system.num_bonds == 2
    assert set(system.atom_types) == {"O", "H"}


def test_from_file_cif():
    system = AtomicSystem.from_file(DATA_DIR / "calcite.cif")
    assert system.num_atoms == 30
    assert set(system.atom_types) == {"C", "Ca", "O"}


def _relabel_axes(structure, perm):
    """Same crystal described with permuted axes (a sign flip keeps it right-handed)."""
    from pymatgen.core import Lattice, Structure

    matrix = structure.lattice.matrix[perm].copy()
    frac = structure.frac_coords[:, perm].copy()
    if np.linalg.det(matrix) < 0:
        matrix[2] *= -1
        frac[:, 2] *= -1
    return Structure(Lattice(matrix), structure.species, frac)


@pytest.mark.parametrize("angles", [(90, 90, 90), (80, 95, 100)])
@pytest.mark.parametrize("perm", [[0, 2, 1], [2, 0, 1], [1, 2, 0]])
def test_pmg_reader_reindex_preserves_crystal(monkeypatch, angles, perm):
    from pymatgen.core import Lattice, Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    from cemd.core._io.formats.pmg import PMGReader

    original = Structure(
        Lattice.from_parameters(6.732, 7.369, 22.68, *angles),
        ["Ca", "O", "O", "Si"],
        [[0.1, 0.2, 0.3], [0.6, 0.1, 0.8], [0.3, 0.7, 0.2], [0.9, 0.5, 0.6]],
    )
    # Simulate a refinement that relabels the axes
    refined = _relabel_axes(original, perm)
    monkeypatch.setattr(
        SpacegroupAnalyzer, "get_refined_structure", lambda self: refined
    )

    with pytest.warns(UserWarning, match="Axes reindexed"):
        topology = PMGReader.read(original, refine=True)

    # The long axis stays along c
    assert topology["box"][2] == pytest.approx(22.68)

    # Same crystal: periodic distances are unchanged
    result = Structure(
        Lattice.from_parameters(*topology["box"]),
        topology["atoms"]["type"].tolist(),
        topology["atoms"][["x", "y", "z"]].to_numpy(),
        coords_are_cartesian=True,
    )
    assert np.allclose(
        np.sort(result.distance_matrix, axis=None),
        np.sort(original.distance_matrix, axis=None),
    )

    # Same handedness: the structure is rotated, not mirrored
    def signed_volume(pos):
        return np.linalg.det(pos[1:4] - pos[0])

    assert signed_volume(topology["atoms"][["x", "y", "z"]].to_numpy()) == pytest.approx(
        signed_volume(original.cart_coords)
    )


def test_from_file_sdf():
    system = AtomicSystem.from_file(DATA_DIR / "ho.sdf")
    assert system.num_atoms == 2
    assert system.num_bonds == 1


def test_from_file_lt():
    system = AtomicSystem.from_file(DATA_DIR / "h2o.lt")
    assert system.num_atoms == 3
    assert system.num_bonds == 2
    assert system.num_angles == 1


def test_from_file_lammps_data_with_type_labels():
    # Regression test: this file uses LAMMPS "Atom Type Labels" (string
    # types directly in the Masses section), which used to crash with
    # `ValueError: invalid literal for int() with base 10: 'C'`.
    system = AtomicSystem.from_file(DATA_DIR / "calcite_ortho.data")
    assert system.num_atoms == 60
    assert set(system.atom_types) == {"C", "Ca", "O"}
    assert system.masses["Ca"] == pytest.approx(40.078, rel=1e-3)


def test_from_file_lammps_data_classic_numeric_types(tmp_path):
    content = (
        "LAMMPS data file\n"
        "2 atoms\n"
        "2 atom types\n"
        "0 10 xlo xhi\n"
        "0 10 ylo yhi\n"
        "0 10 zlo zhi\n"
        "\nMasses\n\n"
        "1 12.011\n"
        "2 15.999\n"
        "\nAtoms # full\n\n"
        "1 1 1 0.5 0.0 0.0 0.0\n"
        "2 1 2 -0.5 1.0 0.0 0.0\n"
    )
    path = tmp_path / "classic.data"
    path.write_text(content)

    system = AtomicSystem.from_file(path)
    # Regression test: total_mass/total_charge/density used to silently
    # come out as 0 for classic numeric LAMMPS types (str/int key mismatch
    # between `atoms["type"]` and the masses/charges dicts).
    assert system.total_mass == pytest.approx(12.011 + 15.999)
    assert system.total_charge == pytest.approx(0.0)
    assert system.density > 0.0


# ---------------------------------------------------------------------------
# write: round-trip through .data and .pdb
# ---------------------------------------------------------------------------


def test_write_unsupported_extension(water_system, tmp_path):
    with pytest.raises(ValueError, match="Unsupported output format"):
        water_system.write(tmp_path / "out.xyz")


def test_write_and_reread_lammps_data_roundtrip(tmp_path):
    # Regression test: the writer emits bond/angle *type labels* directly
    # (e.g. "Hw-Ow") when atom types are strings, but the reader used to
    # assume the type column was always numeric and crashed with
    # `IntCastingNaNError` on its own output.
    system = make_water_system()
    path = tmp_path / "water.data"
    system.write(path)

    reloaded = AtomicSystem.from_file(path)
    assert reloaded.num_atoms == system.num_atoms
    assert reloaded.num_bonds == system.num_bonds
    assert reloaded.num_angles == system.num_angles
    assert set(reloaded.atom_types) == set(system.atom_types)
    assert reloaded.masses["Ow"] == pytest.approx(system.masses["Ow"])
    assert reloaded.charges["Ow"] == pytest.approx(system.charges["Ow"])
    assert set(reloaded.bond_types) == set(system.bond_types)

    np.testing.assert_allclose(
        sorted(reloaded.atoms["x"].tolist()), sorted(system.atoms["x"].tolist())
    )


def test_write_and_reread_pdb_roundtrip(tmp_path):
    system = make_water_system()
    path = tmp_path / "water.pdb"
    system.write(path)

    reloaded = AtomicSystem.from_file(path)
    assert reloaded.num_atoms == system.num_atoms


def test_pdb_roundtrip_does_not_corrupt_charges(tmp_path):
    # Regression test: the standard PDB format has no charge field, but
    # `PDBReader` used to read back the *occupancy* column (always 1.00,
    # written as a fixed placeholder by `PDBWriter`) as "charge" -- so
    # every atom's real charge was silently overwritten with 1.0 on any
    # PDB round-trip, including every Packmol-built system.
    # A plain write()/from_file() round-trip still can't preserve charge
    # (PDB has nowhere to put it), but it must fall back to 0.0, not a
    # fabricated 1.0; the Packmol pipeline restores real charges from the
    # original templates explicitly (see test_build_solution.py).
    system = make_water_system()
    path = tmp_path / "water.pdb"
    system.write(path)

    reloaded = AtomicSystem.from_file(path)
    assert reloaded.charges["Ow"] == pytest.approx(0.0)
    assert reloaded.charges["Hw"] == pytest.approx(0.0)


def test_pdb_roundtrip_preserves_two_letter_elements(tmp_path):
    # Regression test: `PDBWriter` used to fill the element column with the
    # first letter of the atom type, so "Ca" was written as carbon (and "Si"
    # as sulfur). Reading the file back then assigned those wrong masses,
    # and since `elements` is inferred from mass, a calcium silicate came
    # out of any Packmol build as a carbon/sulfur compound. Calcite makes
    # the collision explicit: it holds both "Ca" and "C".
    system = AtomicSystem.from_file(DATA_DIR / "calcite.cif")
    assert system.elements == {"C": "C", "Ca": "Ca", "O": "O"}

    path = tmp_path / "calcite.pdb"
    system.write(path)

    reloaded = AtomicSystem.from_file(path)
    assert reloaded.elements == {"C": "C", "Ca": "Ca", "O": "O"}
    assert reloaded.masses["Ca"] == pytest.approx(system.masses["Ca"])
    assert reloaded.masses["C"] == pytest.approx(system.masses["C"])


def test_pdb_element_column_is_written_upper_case(tmp_path):
    # The PDB standard right-justifies an upper-cased symbol in columns
    # 77-78; the reader normalizes the case on the way back in.
    system = AtomicSystem.from_file(DATA_DIR / "calcite.cif")
    path = tmp_path / "calcite.pdb"
    system.write(path)

    written = {
        line[76:78].strip()
        for line in path.read_text().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    }
    assert written == {"C", "CA", "O"}


# ---------------------------------------------------------------------------
# Converters: to_mda / from_mda / to_pmg / from_smiles
# ---------------------------------------------------------------------------


def test_to_mda_roundtrip_preserves_topology():
    system = make_water_system()
    universe = system.to_mda()
    assert universe.atoms.n_atoms == system.num_atoms
    assert len(universe.bonds) == system.num_bonds

    reloaded = AtomicSystem.from_mda(universe)
    assert reloaded.num_atoms == system.num_atoms
    assert reloaded.num_bonds == system.num_bonds
    assert set(reloaded.atom_types) == set(system.atom_types)


def test_to_pmg_returns_pymatgen_structure():
    from pymatgen.core import Structure

    system = AtomicSystem.from_file(DATA_DIR / "calcite.cif")
    structure = system.to_pmg()
    assert isinstance(structure, Structure)
    assert len(structure) == system.num_atoms


def test_from_smiles_water():
    system = AtomicSystem.from_smiles("O")
    assert system.num_atoms == 3
    assert set(system.atom_types) == {"O", "H"}
    assert system.num_bonds == 2


def test_from_smiles_invalid_raises():
    with pytest.raises(ValueError, match="Invalid SMILES"):
        AtomicSystem.from_smiles("not_a_smiles!!")
