"""Tests for cemd.build.SurfaceBuilder, using the real calcite.cif fixture."""

from __future__ import annotations

import pytest

from cemd import AtomicSystem
from cemd.build import SurfaceBuilder

from conftest import DATA_DIR


@pytest.fixture(scope="module")
def calcite() -> AtomicSystem:
    return AtomicSystem.from_file(DATA_DIR / "calcite.cif")


def test_surface_builder_repr(calcite):
    builder = SurfaceBuilder(calcite)
    text = repr(builder)
    assert "SurfaceBuilder" in text
    assert "sites   : 30" in text


def test_build_generates_slabs_with_vacuum(calcite):
    builder = SurfaceBuilder(calcite)
    slabs, shifts, dipoles, broken = builder.build(
        miller_indices=[1, 0, 4], min_slab_size=10.0, min_vacuum_size=10.0
    )

    assert len(slabs) >= 1
    assert len(slabs) == len(shifts) == len(dipoles)
    assert all(isinstance(s, AtomicSystem) for s in slabs)
    assert slabs[0].num_atoms > calcite.num_atoms
    assert broken >= 0

    # `build()` also stashes its own return value on `.result`.
    assert builder.result == (slabs, shifts, dipoles, broken)


def test_build_refuses_an_edited_crystal(calcite):
    # The structure cached at read time would bring deleted atoms back, and
    # the edited system is no bulk to cut.
    system = calcite.copy()
    system.remove_atoms(system.atoms.index[system.atoms["type"] == "C"].tolist())

    with pytest.raises(ValueError, match="edited since it was read"):
        SurfaceBuilder(system)

    # Asking twice gives the same answer: the stale structure is kept.
    with pytest.raises(ValueError, match="edited since it was read"):
        SurfaceBuilder(system)


def test_build_refuses_a_slab(calcite):
    slabs, *_ = SurfaceBuilder(calcite).build(
        miller_indices=[1, 0, 4], min_slab_size=10.0, min_vacuum_size=10.0
    )
    with pytest.raises(ValueError, match="already a slab"):
        SurfaceBuilder(slabs[0])


def test_build_accepts_a_typed_crystal(calcite):
    from cemd.topology.predefined import CLAYFF_RULES

    # Typing moves nothing, so the structure read stays valid.
    system = calcite.copy()
    system.set_topology(CLAYFF_RULES)
    assert "Oc" in system.atom_types
    assert SurfaceBuilder(system)._structure.composition.reduced_formula == "CaCO3"


def test_structure_of_a_typed_system_uses_elements(calcite, tmp_path):
    from cemd.topology.predefined import CLAYFF_RULES

    system = calcite.copy()
    system.set_topology(CLAYFF_RULES)

    # Read back from a .data file, the system has no cached structure: the
    # builder converts it, and the force-field types are no pymatgen species.
    system.write(tmp_path / "calcite.data")
    system = AtomicSystem.from_file(tmp_path / "calcite.data")
    builder = SurfaceBuilder(system)
    assert builder._structure.composition.reduced_formula == "CaCO3"

    slabs, *_ = builder.build(
        miller_indices=[1, 0, 4], min_slab_size=10.0, min_vacuum_size=10.0
    )
    assert set(slabs[0].atom_types) == {"Ca", "C", "O"}


def test_build_rejects_structure_with_too_few_atoms():
    atoms_df = calcite_single_atom_df()
    tiny = AtomicSystem(
        {
            "atoms": atoms_df,
            "box": [10.0, 10.0, 10.0, 90.0, 90.0, 90.0],
            "masses": {"Ca": 40.078},
            "charges": {},
        }
    )
    builder = SurfaceBuilder(tiny)
    with pytest.raises(ValueError, match="at least 2 atoms"):
        builder.build(miller_indices=[1, 0, 0])


def calcite_single_atom_df():
    import pandas as pd

    return pd.DataFrame(
        {"type": ["Ca"], "charge": [0.0], "x": [0.0], "y": [0.0], "z": [0.0]},
        index=[1],
    )


def _magnesian_calcite(ca_site):
    """Calcite whose Ca sites carry `ca_site`, like COD 9017359 (Ca0.99Mg0.01)."""
    from pymatgen.core import Structure

    structure = Structure.from_file(DATA_DIR / "calcite.cif")
    structure.replace_species({"Ca": ca_site})
    return structure


def test_a_disordered_crystal_is_cut_with_its_main_species():
    # SlabGenerator refuses disordered structures: "species property only
    # supports ordered structures".
    with pytest.warns(UserWarning, match="main species: Ca0.99 Mg0.01"):
        system = AtomicSystem.from_pmg(_magnesian_calcite({"Ca": 0.99, "Mg": 0.01}))

    slabs, *_ = SurfaceBuilder(system).build(miller_indices=(1, 0, 4))

    assert slabs
    types = set(slabs[0].atom_types)
    assert types == {"Ca", "C", "O"}


def test_a_disordered_site_takes_its_majority_species_not_its_first():
    with pytest.warns(UserWarning, match="main species"):
        system = AtomicSystem.from_pmg(_magnesian_calcite({"Mg": 0.3, "Ca": 0.7}))

    assert "Mg" not in system.atom_types
    assert system.get_count("Ca") > 0
