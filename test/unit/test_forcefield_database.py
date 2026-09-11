"""Tests for cemd.forcefield.ForceFieldDatabase against the real bundled
force-field database (no mocking -- these are the actual TOML/lt/prm/frc
files shipped with the package).
"""

from __future__ import annotations

import pandas as pd
import pytest

from cemd.forcefield import ForceFieldDatabase
from cemd.forcefield.models import AtomType, HarmonicAngleParams, HarmonicBondParams


@pytest.fixture(scope="module")
def db() -> ForceFieldDatabase:
    return ForceFieldDatabase()


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def test_database_loads_multiple_models(db):
    names = db.get_model_names()
    assert len(names) > 5
    assert "spc" in names


def test_database_loads_atoms_bonds_angles(db):
    assert len(db.atom) > 0
    assert len(db.bond) > 0
    assert len(db.angle) > 0


def test_custom_db_dir_starts_empty(tmp_path):
    custom = ForceFieldDatabase(db_dir=tmp_path)
    assert custom.atom == {}
    assert custom.get_model_names() == []


# ---------------------------------------------------------------------------
# get_atom_type
# ---------------------------------------------------------------------------


def test_get_atom_type_by_full_name(db):
    atom = db.get_atom_type("spc.ospc")
    assert isinstance(atom, AtomType)
    assert atom.element == "O"
    assert atom.charge == pytest.approx(-0.82)
    assert atom.model == "spc"


def test_get_atom_type_missing_returns_none(db):
    assert db.get_atom_type("not_a_real_type_xyz") is None


def test_get_model(db):
    model = db.get_model("spc")
    assert model is not None
    assert model.name == "SPC"


def test_get_model_missing_returns_none(db):
    assert db.get_model("not_a_real_model") is None


def test_get_atom_types_for_model(db):
    types = db.get_atom_types_for_model("spc")
    assert all(t.startswith("spc.") for t in types)
    assert "spc.ospc" in types
    assert "spc.hspc" in types


# ---------------------------------------------------------------------------
# Pair / triple / quadruple parameter lookups
# ---------------------------------------------------------------------------


def test_get_lj_self_interaction(db):
    params = db.get_lj("spc.ospc", "spc.ospc")
    assert params is not None
    assert params.epsilon == pytest.approx(0.15535)
    assert params.sigma == pytest.approx(3.166)


def test_get_lj_symmetric_short_names_with_explicit_model(db):
    forward = db.get_lj("ospc", "hspc", model="spc")
    # SPC has no o-h cross LJ term defined directly, but a same-order
    # self-pair lookup on a real cross key should still work either way.
    reverse = db.get_lj("hspc", "ospc", model="spc")
    assert forward == reverse


def test_get_lj_missing_returns_none(db):
    assert db.get_lj("spc.ospc", "not_a_real_type", model="spc") is None


def test_get_bond_full_names(db):
    params = db.get_bond("spc.hspc", "spc.ospc")
    assert isinstance(params, HarmonicBondParams)
    assert params.k == pytest.approx(554.1349)
    assert params.r0 == pytest.approx(1.0)


def test_get_bond_reversed_pair_matches(db):
    forward = db.get_bond("hspc", "ospc", model="spc")
    reverse = db.get_bond("ospc", "hspc", model="spc")
    assert forward == reverse
    assert forward is not None


def test_get_angle_short_names_with_model(db):
    params = db.get_angle("hspc", "ospc", "hspc", model="spc")
    assert isinstance(params, HarmonicAngleParams)
    assert params.theta0 == pytest.approx(109.47)


def test_get_angle_missing_returns_none(db):
    assert db.get_angle("zz", "yy", "xx") is None


def test_get_dihedral_and_improper_missing_return_none(db):
    assert db.get_dihedral("a", "b", "c", "d") is None
    assert db.get_improper("a", "b", "c", "d") is None


def test_get_bondbond_bondangle_present_for_iff_cvff(db):
    # Regression coverage for the class2 cross-term accessors specifically
    # (as opposed to just get_bond/get_angle).
    key = next(iter(db.bondbond))
    model, pair = key.split(".", 1)
    t1, t2, t3 = pair.split("-")
    assert db.get_bondbond(t1, t2, t3, model=model) is not None


# ---------------------------------------------------------------------------
# to_dataframes / clear
# ---------------------------------------------------------------------------


def test_to_dataframes_returns_populated_frames(db):
    dfs = db.to_dataframes()
    assert isinstance(dfs["atoms"], pd.DataFrame)
    assert not dfs["atoms"].empty
    assert {"type", "full_type", "element", "model"}.issubset(dfs["atoms"].columns)
    assert not dfs["lj"].empty
    assert not dfs["bonds"].empty


def test_clear_empties_all_collections():
    fresh = ForceFieldDatabase()
    assert len(fresh.atom) > 0
    fresh.clear()
    assert fresh.atom == {}
    assert fresh.lj == {}
    assert fresh.bond == {}
    assert fresh.angle == {}
    assert fresh.models == {}


# ---------------------------------------------------------------------------
# GROMOS parser: scientific notation and per-atom charges
# ---------------------------------------------------------------------------


def test_gromos_angles_are_parsed_in_degrees(db):
    # Regression test: the equilibrium angle is written in scientific
    # notation ("1.2600000000E+02" for 126 deg), but the regex group for it
    # left out "E+-" and so matched only the mantissa. Every GROMOS angle
    # came out a factor of ten or a hundred too small -- 126 deg read as
    # 1.26, the tetrahedral 109.5 as 1.095.
    angles = [v.theta0 for k, v in db.angle.items() if k.startswith("gromos.")]
    assert angles, "the bundled GROMOS database defines angles"

    # Every equilibrium angle must be a real angle, not a mantissa.
    assert min(angles) > 10.0
    assert max(angles) <= 180.0
    # The tetrahedral angle is in there, and it is 109.5, not 1.095.
    assert any(v == pytest.approx(109.5) for v in angles)


def test_gromos_bond_lengths_are_parsed_with_their_exponent(db):
    # Same fragility on the bond side: it happened to be harmless because
    # every r0 exponent in the bundled file is E+00, but the group has to
    # carry the exponent for the value to be trustworthy.
    bonds = [
        v.r0
        for k, v in db.bond.items()
        if k.startswith("gromos.") and not k.endswith("excl")
    ]
    assert bonds
    assert min(bonds) > 0.5
    assert max(bonds) < 5.0


def test_gromos_atom_types_declare_no_per_type_charge(db):
    # GROMOS carries partial charges per atom in the molecular topology,
    # not per atom type. Reporting 0.0 instead of None made
    # `set_ff_from_database()` overwrite a molecule's real partial charges
    # with zeros.
    gromos_atoms = [v for k, v in db.atom.items() if k.startswith("gromos.")]
    assert gromos_atoms
    assert all(a.charge is None for a in gromos_atoms)
