#
# This file is part of the CEMD distribution
# Copyright (c) 2022-2026 Jérôme Claverie.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
"""Organic molecules typed for CGenFF from a SMILES string.

Only the chemical families whose CGenFF model compounds give both the types
and the charges are handled -- alkanes, carboxylates, carboxylic esters and
ethers, enough for comb polymers such as methacrylate PCE superplasticizers
(methacrylate backbone, MPEG ester side chains). Any other atom raises an
error instead of being typed by analogy.

Types and charges follow `_data/cgenff_charges.toml`; bonded and LJ terms come
from the CGenFF parameter file of the database (``cgenff.prm``), wildcards
included. Run the result with::

    pair_style      lj/charmmfsw/coul/long 10.0 12.0
    bond_style      harmonic
    angle_style     charmm
    dihedral_style  charmmfsw
    improper_style  harmonic
    special_bonds   lj/coul 0.0 0.0 0.0 dihedral yes

``dihedral yes`` keeps the full 1-4 pair term for atoms three bonds apart that
are not the ends of a dihedral (a mineral framework next to the molecule),
while the dihedrals of the molecule add their own (CHARMM) 1-4 term.
"""

from __future__ import annotations

import tomllib
import warnings
from itertools import combinations, product
from pathlib import Path

import pandas as pd

from ..core._format import canonical_ff_type

_DATA = Path(__file__).parent / "_data" / "cgenff_charges.toml"


def _rules() -> dict:
    with open(_DATA, "rb") as f:
        return tomllib.load(f)


def cgenff_types(mol) -> list[str]:
    """CGenFF atom type of every atom of an RDKit molecule (explicit H)."""
    types = []
    for atom in mol.GetAtoms():
        types.append(_type_of(atom))
    return types


def _type_of(atom) -> str:
    sym = atom.GetSymbol()
    nbrs = list(atom.GetNeighbors())
    hyb = str(atom.GetHybridization())
    if sym == "H":
        (c,) = nbrs
        if c.GetSymbol() != "C" or str(c.GetHybridization()) != "SP3":
            raise NotImplementedError(f"H on {c.GetSymbol()} ({c.GetHybridization()})")
        return f"HGA{c.GetTotalNumHs(includeNeighbors=True)}"
    if sym == "C" and hyb == "SP3":
        if any(n.GetSymbol() not in ("C", "H", "O") for n in nbrs):
            raise NotImplementedError("sp3 C bonded to an element other than C, H, O")
        return {3: "CG331", 2: "CG321", 1: "CG311", 0: "CG301"}[atom.GetTotalNumHs(includeNeighbors=True)]
    if sym == "C" and hyb == "SP2":
        oxygens = [n for n in nbrs if n.GetSymbol() == "O"]
        if len(oxygens) == 2 and len(nbrs) == 3:
            single = [o for o in oxygens if o.GetDegree() == 1 and o.GetFormalCharge() == -1]
            ester = [o for o in oxygens if o.GetDegree() == 2]
            if single:
                return "CG2O3"
            if ester:
                return "CG2O2"
        raise NotImplementedError("sp2 C other than a carboxylate or an ester carbonyl")
    if sym == "O":
        if atom.GetDegree() == 1:
            (c,) = nbrs
            t = _type_of(c)
            if t == "CG2O3":
                return "OG2D2"
            if t == "CG2O2":
                return "OG2D1"
            raise NotImplementedError("terminal O outside a carboxylate or an ester")
        if atom.GetDegree() == 2 and all(n.GetSymbol() == "C" for n in nbrs):
            kinds = sorted(str(n.GetHybridization()) for n in nbrs)
            if kinds == ["SP3", "SP3"]:
                return "OG301"
            if kinds == ["SP2", "SP3"]:
                sp2 = next(n for n in nbrs if str(n.GetHybridization()) == "SP2")
                if _type_of(sp2) == "CG2O2":
                    return "OG302"
        raise NotImplementedError("O other than ether, ester or carboxylate")
    raise NotImplementedError(f"atom {sym} ({hyb}) not handled by the CGenFF rules")


def cgenff_charges(mol, types: list[str]) -> list[float]:
    """CGenFF charges from the base values and the neighbour increments."""
    rules = _rules()
    base, inc = rules["base"], rules["increment"]
    charges = []
    for atom, t in zip(mol.GetAtoms(), types):
        q = base[t]
        if t.startswith("CG3"):
            for n in atom.GetNeighbors():
                q += inc.get(types[n.GetIdx()], 0.0)
        charges.append(round(q, 6))
    return charges


def _impropers(mol, types: list[str]) -> list[tuple[int, ...]]:
    """CGenFF impropers, central atom first, neighbours in the parameter
    file's order: carboxylate ``CG2O3 OG2D2 OG2D2 C`` and ester
    ``CG2O2 C OG2D1 OG302``."""
    out = []
    for atom, t in zip(mol.GetAtoms(), types):
        i = atom.GetIdx()
        nb = {types[n.GetIdx()]: [] for n in atom.GetNeighbors()}
        for n in atom.GetNeighbors():
            nb[types[n.GetIdx()]].append(n.GetIdx())
        if t == "CG2O3":
            (o1, o2), (c,) = nb["OG2D2"], [v[0] for k, v in nb.items() if k != "OG2D2"]
            out.append((i, o1, o2, c))
        elif t == "CG2O2":
            (c,) = [v[0] for k, v in nb.items() if k not in ("OG2D1", "OG302")]
            out.append((i, c, nb["OG2D1"][0], nb["OG302"][0]))
    return out


def cgenff_molecule(smiles: str, seed: int = 0, optimize: bool = True):
    """Build an AtomicSystem typed and parametrized with CGenFF from SMILES.

    Parameters
    ----------
    smiles
        The molecule, hydrogens implicit; a carboxylate is written ``C(=O)[O-]``.
    seed
        Seed of RDKit's conformer embedding.
    optimize
        Relax the conformer with MMFF (geometry only, the parameters stay CGenFF).

    Returns
    -------
    AtomicSystem
        Atom types are the CGenFF types; charges are per atom (CGenFF groups);
        bonds, angles, dihedrals and impropers carry CGenFF parameters.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from ..core.atomic_system import AtomicSystem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) == -1:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) == -1:
            raise RuntimeError(f"RDKit could not embed {smiles}")
    if optimize:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)

    types = cgenff_types(mol)
    charges = cgenff_charges(mol, types)
    total = sum(charges)
    formal = Chem.GetFormalCharge(mol)
    if abs(total - formal) > 1e-6:
        raise RuntimeError(f"CGenFF charges sum to {total:+.4f}, formal charge {formal}")

    conf = mol.GetConformer()
    atoms = pd.DataFrame(
        {
            "type": types,
            "charge": charges,
            "x": [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())],
            "y": [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())],
            "z": [conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())],
        },
        index=pd.RangeIndex(1, mol.GetNumAtoms() + 1, name="id"),
    )

    def frame(rows, kind):
        if not rows:
            return None
        n = len(rows[0])
        records = []
        for idx in rows:
            tt = tuple(types[i] for i in idx)
            if kind == "dihedral":
                label = canonical_ff_type(tt, "dihedral")
                if label != "-".join(tt):
                    idx = idx[::-1]
            elif kind == "improper":
                label = "-".join(tt)  # CHARMM order, central atom first
            else:
                label = canonical_ff_type(tt, kind)
            records.append({"type": label, **{f"atom_{k + 1}": i + 1 for k, i in enumerate(idx)}})
        df = pd.DataFrame(records)
        df.index = pd.RangeIndex(1, len(df) + 1, name="id")
        return df

    adj = {a.GetIdx(): sorted(n.GetIdx() for n in a.GetNeighbors()) for a in mol.GetAtoms()}
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    angles = [(i, j, k) for j in adj for i, k in combinations(adj[j], 2)]
    dihedrals = [
        (i, j, k, l)
        for j, k in bonds
        for i in adj[j] if i != k
        for l in adj[k] if l != j and l != i
    ]

    coords = atoms[["x", "y", "z"]].to_numpy()
    lo, hi = coords.min(axis=0) - 10, coords.max(axis=0) + 10
    topology = {
        "box": ((lo[0], hi[0]), (lo[1], hi[1]), (lo[2], hi[2]), (0.0, 0.0, 0.0)),
        "masses": {},
        "charges": {},
        "atoms": atoms,
        "bonds": frame(bonds, "bond"),
        "angles": frame(angles, "angle"),
        "dihedrals": frame(dihedrals, "dihedral"),
        "impropers": frame(_impropers(mol, types), "improper"),
        "velocities": None,
        "atom_style": "full",
    }
    system = AtomicSystem(topology)
    with warnings.catch_warnings():
        # Missing bonded terms are reported below, after the analogies
        warnings.filterwarnings("ignore", message=".*parameters not found.*")
        system.set_ff_from_database(atom_assignments={t: f"cgenff.{t}" for t in set(types)})
    system.cgenff_analogies = _fill_by_analogy(system)
    return system


def _fill_by_analogy(system) -> list[tuple[str, str, str]]:
    """Give the bonded types left without CGenFF parameters those of the
    closest analogue (``[analogy]`` table). Returns (kind, type, analogue)."""
    from ..forcefield.forcefield_database import ForceFieldDatabase

    db = ForceFieldDatabase()
    analogy = _rules().get("analogy", {})
    params = system._ff_params
    done, missing = [], []
    kinds = [
        ("bond", system.bond_types, params.bond, db.get_bond, system.set_bond_params),
        ("angle", system.angle_types, params.angle, db.get_angle, system.set_angle_params),
        ("dihedral", system.dihedral_types, params.dihedral, db.get_dihedral,
         system.set_dihedral_params),
        ("improper", system.improper_types, params.improper, db.get_improper,
         system.set_improper_params),
    ]
    for kind, labels, store, get, setter in kinds:
        for label in labels:
            if label in store:
                continue
            found = None
            for candidate in _analogues(label.split("-"), analogy):
                found = get(*[f"cgenff.{t}" for t in candidate])
                if found is not None:
                    setter(label, found)
                    done.append((kind, label, "-".join(candidate)))
                    break
            if found is None:
                missing.append(f"{kind} {label}")
    if missing:
        raise RuntimeError(f"CGenFF terms without parameters or analogue: {missing}")
    if done:
        warnings.warn(
            f"{len(done)} CGenFF bonded terms taken by analogy (see "
            f"system.cgenff_analogies): "
            + ", ".join(f"{t} <- {a}" for _, t, a in done),
            UserWarning,
            stacklevel=3,
        )
    return done


def _analogues(types: list[str], analogy: dict[str, list[str]]):
    """Type tuples with some types replaced by their analogues, the fewest
    replacements first, then in the order of the analogue lists."""
    options = [[(t, 0)] + [(a, r + 1) for r, a in enumerate(analogy.get(t, []))] for t in types]
    combos = list(product(*options))
    combos.sort(key=lambda c: (sum(r > 0 for _, r in c), sum(r for _, r in c)))
    for combo in combos[1:]:  # the first one is the unchanged tuple
        yield [t for t, _ in combo]
