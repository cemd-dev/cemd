"""CGenFF / CHARMM36 support: parser, lookups, typing, charges and LAMMPS output."""

import warnings

import pytest

from cemd.forcefield.forcefield_database import ForceFieldDatabase
from cemd.forcefield.models import (
    CHARMMAngleParams,
    CHARMMDihedralParams,
    CHARMMDihedralSeries,
    LJParams,
)

pytest.importorskip("rdkit")

# Charges of the CGenFF 4.6 model compounds (top_all36_cgenff.rtf), by type
RTF_CHARGES = {
    "CC(=O)[O-]": {"CG331": -0.37, "CG2O3": 0.62, "OG2D2": -0.76, "HGA3": 0.09},  # ACET
    "CC(=O)OC": None,  # MAS, two CG331 with different charges: checked below
    "COCCOC": {"CG331": -0.10, "OG301": -0.34, "CG321": -0.01, "HGA3": 0.09,
               "HGA2": 0.09},  # DME
    "CCC": {"CG331": -0.27, "CG321": -0.18, "HGA3": 0.09, "HGA2": 0.09},  # PRPA
    "CC(C)C": {"CG311": -0.09, "CG331": -0.27, "HGA1": 0.09, "HGA3": 0.09},  # IBUT
    "CC(C)(C)C": {"CG301": 0.0, "CG331": -0.27, "HGA3": 0.09},  # NEOP
}


@pytest.fixture(scope="module")
def db():
    return ForceFieldDatabase()


def build(smiles):
    from cemd.build.cgenff import cgenff_molecule

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cgenff_molecule(smiles, seed=1)


# -- parser and database -----------------------------------------------------


def test_cgenff_and_ion_models_are_loaded(db):
    assert "cgenff" in db.models and "charmm36_ions" in db.models
    ca = db.atom["charmm36_ions.CAL"]
    assert ca.element == "Ca" and ca.charge == 2.0
    assert db.atom["cgenff.OG2D2"].element == "O"


def test_lj_keeps_the_charmm_14_values(db):
    lj = db.lj["cgenff.CG331-CG331"]
    assert lj.epsilon == pytest.approx(0.078)
    assert lj.sigma == pytest.approx(2 * 2.05 / 2 ** (1 / 6))
    assert lj.epsilon14 == pytest.approx(0.01)
    assert lj.sigma14 == pytest.approx(2 * 1.9 / 2 ** (1 / 6))
    assert db.lj["cgenff.OG2D2-OG2D2"].epsilon14 is None


def test_nbfix_section_does_not_overwrite_the_self_terms(db):
    # NBFIX lines "OG2D1 CLGR1 ..." used to be read as LJ lines of OG2D1
    lj = db.lj["cgenff.OG2D1-OG2D1"]
    assert lj.epsilon == pytest.approx(0.12) and lj.sigma == pytest.approx(3.4 / 2 ** (1 / 6))


def test_nbfix_between_two_models(db):
    for pair in (("charmm36_ions.CAL", "cgenff.OG2D2"), ("cgenff.OG2D2", "charmm36_ions.CAL")):
        lj = db.get_lj(*pair)
        assert lj.ref == "NBFIX"
        assert lj.epsilon == pytest.approx(0.12)
        assert lj.sigma == pytest.approx(3.232 / 2 ** (1 / 6))


def test_angles_keep_urey_bradley(db):
    a = db.get_angle("cgenff.HGA3", "cgenff.CG331", "cgenff.HGA3")
    assert isinstance(a, CHARMMAngleParams)
    assert (a.k, a.theta0, a.k_ub, a.r_ub) == (35.5, 108.4, 5.4, 1.802)


def test_multi_term_dihedral_is_a_series_weighted_once(db):
    d = db.get_dihedral("cgenff.OG301", "cgenff.CG321", "cgenff.CG321", "cgenff.OG301")
    assert isinstance(d, CHARMMDihedralSeries)
    assert [(t.k, t.n, t.d) for t in d.terms] == [(0.25, 1, -1), (1.24, 2, 1)]
    assert [t.w for t in d.terms] == [1.0, 0.0]


def test_wildcard_dihedral(db):
    # only "X CG2O3 ..."-type entries exist for some central bonds
    d = db.get_dihedral("cgenff.HGA3", "cgenff.CG331", "cgenff.CG301", "cgenff.CG2O3")
    assert d is not None


def test_improper_found_whatever_the_neighbour_order(db):
    for order in (("OG2D2", "OG2D2", "CG331"), ("CG331", "OG2D2", "OG2D2")):
        p = db.get_improper("cgenff.CG2O3", *[f"cgenff.{t}" for t in order])
        assert p is not None and p.k == 96.0


# -- typing and charges ----------------------------------------------------------


@pytest.mark.parametrize("smiles", [s for s, q in RTF_CHARGES.items() if q])
def test_charges_reproduce_the_cgenff_model_compounds(smiles):
    s = build(smiles)
    for t, q in zip(s.atoms["type"], s.atoms["charge"]):
        assert q == pytest.approx(RTF_CHARGES[smiles][t]), t


def test_methyl_acetate_charges():
    s = build("CC(=O)OC")
    q = sorted(round(x, 2) for x, t in zip(s.atoms["charge"], s.atoms["type"]) if t == "CG331")
    assert q == [-0.31, -0.01]  # acyl and methoxy CH3 of MAS
    by_type = dict(zip(s.atoms["type"], s.atoms["charge"]))
    assert by_type["CG2O2"] == pytest.approx(0.90)
    assert by_type["OG2D1"] == pytest.approx(-0.63)
    assert by_type["OG302"] == pytest.approx(-0.49)


def test_pce_segment_is_complete_and_charged():
    s = build("CC(C)(C(=O)[O-])CC(C)(C)C(=O)OCCOCCOCCOC")
    assert s.atoms["charge"].sum() == pytest.approx(-1.0)
    p = s._ff_params
    for kind in ("bond", "angle", "dihedral", "improper"):
        assert set(getattr(s, f"{kind}_types")) <= set(getattr(p, kind)), kind
    analogues = {label: a for _, label, a in s.cgenff_analogies}
    assert analogues["OG301-CG321-CG321-OG302"] == "OG301-CG321-CG321-OG301"
    assert analogues["CG2O2-CG301"] == "CG2O2-CG311"


def test_unsupported_chemistry_is_refused():
    from cemd.build.cgenff import cgenff_molecule

    with pytest.raises(NotImplementedError):
        cgenff_molecule("CCN")


# -- LAMMPS output ---------------------------------------------------------------


def test_lammps_data_charmm_terms(tmp_path):
    s = build("COCCOC")
    path = tmp_path / "dme.data"
    s.write(str(path))
    text = path.read_text()
    coeffs = text.split("Dihedral Coeffs")[1].split("\n\n")[1].splitlines()
    # the O-C-C-O series is split in two types, the phase written in degrees
    occo = [line for line in coeffs if "OG301-CG321-CG321-OG301" in line]
    assert len(occo) == 2
    fields = [line.split("#")[0].split() for line in occo]
    assert [f[1:] for f in fields] == [["0.2500", "1", "180", "1.0000"],
                                      ["1.2400", "2", "0", "0.0000"]]
    # a dihedral is listed once per term of its series
    n_dih = int(text.split(" dihedrals\n")[0].rsplit("\n", 1)[1])
    params = s._ff_params.dihedral
    expected = sum(
        len(params[t].terms) if isinstance(params[t], CHARMMDihedralSeries) else 1
        for t in s.dihedrals["type"]
    )
    assert n_dih == expected > s.num_dihedrals
    # angles in angle_style charmm form, pair lines with the 1-4 values
    angle_line = text.split("Angle Coeffs")[1].split("\n\n")[1].splitlines()[0]
    assert len(angle_line.split("#")[0].split()) == 5
    pair_lines = text.split("PairIJ Coeffs")[1].split("\n\n")[1].splitlines()
    assert any(len(line.split("#")[0].split()) == 6 for line in pair_lines)


def test_lammps_reader_reads_charmm_terms_back(tmp_path):
    from cemd.core._io.formats.lammps import LAMMPSReader

    s = build("CC(=O)OC")
    path = tmp_path / "mas.data"
    s.write(str(path))
    topo = LAMMPSReader.read(str(path))
    assert all(isinstance(p, CHARMMAngleParams) for p in topo["angle_params"].values())
    assert all(isinstance(p, CHARMMDihedralParams) for p in topo["dihedral_params"].values())
    assert any(isinstance(p, LJParams) and p.epsilon14 is not None
               for p in topo["pair_params"].values())
