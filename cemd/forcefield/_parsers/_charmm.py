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
"""Parser for CHARMM parameter files (.prm) and toppar streams (.str).

Used for CGenFF and the CHARMM36 water/ion stream. Unlike the INTERFACE FF
parser (``_iff_charmm``), it keeps every CHARMM term:

* MASS lines (ATOMS section, or the RTF part of a stream) give the types,
  their masses and elements;
* ANGLES keep their Urey-Bradley term (`CHARMMAngleParams`);
* DIHEDRALS repeated for the same four types with several multiplicities are
  gathered in a `CHARMMDihedralSeries`; wildcards ``X`` are kept as such and
  resolved at lookup time (`ForceFieldDatabase.get_dihedral`);
* NONBONDED keeps the special 1-4 values (``epsilon14``/``sigma14``);
* NBFIX pairs are stored as explicit LJ pairs, which take precedence over
  mixing;
* in a stream, single-atom residues (ions) give the type its charge.

CHARMM LJ is given as (epsilon <= 0, Rmin/2); NBFIX as (epsilon, Rmin). Both
are converted to the positive epsilon and sigma = Rmin / 2^(1/6) used by cemd.
"""

from __future__ import annotations

from ..._constants import MASSES_DICT
from ..models import (
    AtomType,
    CHARMMAngleParams,
    CHARMMDihedralParams,
    CHARMMDihedralSeries,
    HarmonicBondParams,
    HarmonicImproperParams,
    LJParams,
)
from ._base import BaseForceFieldParser, ParseResult

_SIXTH_ROOT_2 = 2 ** (1 / 6)

# Section headers (first word of the line) -> internal section name.
# None: section ignored.
_SECTIONS = {
    "ATOMS": "atoms",
    "BONDS": "bonds",
    "ANGLES": "angles",
    "THETAS": "angles",
    "DIHEDRALS": "dihedrals",
    "PHI": "dihedrals",
    "IMPROPER": "impropers",
    "IMPROPERS": "impropers",
    "IMPHI": "impropers",
    "NONBONDED": "nonbonded",
    "NBONDED": "nonbonded",
    "NBFIX": "nbfix",
    "HBOND": None,
    "CMAP": None,
    "END": None,
}


def _element_from_mass(mass: float) -> str:
    """Element whose standard mass is closest to ``mass`` (within 0.6 u)."""
    best, err = "X", 0.6
    for element, m in MASSES_DICT.items():
        d = abs(float(m) - mass)
        if d < err:
            best, err = element, d
    return best


class CHARMMParser(BaseForceFieldParser):
    """Parser for a CHARMM .prm file or a toppar .str stream."""

    def __init__(self, model_name: str, metadata: dict | None = None):
        super().__init__()
        self.model_name = model_name
        self.metadata = metadata or {}

    def _strip_comments(self, line: str) -> str:
        return line.split("!", 1)[0]

    def parse_file(self, filepath: str) -> ParseResult:
        with open(filepath, encoding="utf-8", errors="replace") as f:
            return self.parse(f.read())

    def parse(self, content: str) -> ParseResult:
        result = ParseResult(model_name=self.model_name)
        result.metadata = {"mixing": "arithmetic", **self.metadata}
        section = None
        in_rtf = False
        residue: list | None = None  # [name, charge, [(atom type, charge)]]
        continuation = False

        for raw in content.split("\n"):
            line = self._strip_comments(raw).strip()
            if not line or line.startswith("*"):
                continue
            words = line.split()
            head = words[0].upper()

            # A header continued on the next line ends with "-"
            # (e.g. "NONBONDED nbxmod 5 atom ... -").
            if continuation:
                continuation = line.endswith("-")
                continue

            # Stream files: "read rtf card" / "read para card" blocks
            if head == "READ" and len(words) > 1:
                in_rtf = words[1].upper().startswith("RTF")
                section = None
                continue

            if head == "MASS" and len(words) >= 4:
                self._parse_mass(words, result)
                continue

            if in_rtf:
                residue = self._parse_rtf_line(words, head, residue, result)
                continue

            if head in _SECTIONS:
                section = _SECTIONS[head]
                continuation = line.endswith("-")
                continue

            handler = {
                "bonds": self._parse_bond,
                "angles": self._parse_angle,
                "dihedrals": self._parse_dihedral,
                "impropers": self._parse_improper,
                "nonbonded": self._parse_nonbonded,
                "nbfix": self._parse_nbfix,
            }.get(section)
            if handler is not None:
                handler(words, result)

        self._close_residue(residue, result)
        return result

    # -- types -------------------------------------------------------------

    def _parse_mass(self, words: list[str], result: ParseResult) -> None:
        # MASS  index  type  mass  [element]
        atom_type, mass = words[2], float(words[3])
        element = words[4] if len(words) > 4 and words[4].isalpha() else None
        element = (element.capitalize() if element else None) or _element_from_mass(mass)
        previous = result.atoms.get(atom_type)
        result.atoms[atom_type] = AtomType(
            element=element,
            charge=previous.charge if previous else None,
            mass=mass,
            environment=f"CHARMM type {atom_type}",
            model=result.model_name,
        )

    def _parse_rtf_line(self, words, head, residue, result):
        """Charges of single-atom residues (ions) of a stream's RTF part."""
        if head in ("RESI", "PRES"):
            self._close_residue(residue, result)
            charge = float(words[2]) if len(words) > 2 and _is_float(words[2]) else 0.0
            return [words[1], charge, []]
        if head == "ATOM" and residue is not None and len(words) >= 4:
            residue[2].append((words[2], float(words[3])))
        return residue

    @staticmethod
    def _close_residue(residue, result: ParseResult) -> None:
        if residue is None or len(residue[2]) != 1:
            return
        atom_type, charge = residue[2][0]
        if atom_type in result.atoms:
            result.atoms[atom_type].charge = charge

    # -- bonded ------------------------------------------------------------

    def _parse_bond(self, w: list[str], result: ParseResult) -> None:
        if len(w) < 4:
            return
        i, j = w[0], w[1]
        key = f"{i}-{j}" if i <= j else f"{j}-{i}"
        result.bonds[key] = HarmonicBondParams(
            k=float(w[2]), r0=float(w[3]), model=result.model_name
        )

    def _parse_angle(self, w: list[str], result: ParseResult) -> None:
        if len(w) < 5:
            return
        k_ub, r_ub = (float(w[5]), float(w[6])) if len(w) >= 7 else (0.0, 0.0)
        a, b, c = w[0], w[1], w[2]
        key = f"{a}-{b}-{c}" if a <= c else f"{c}-{b}-{a}"
        result.angles[key] = CHARMMAngleParams(
            k=float(w[3]), theta0=float(w[4]), k_ub=k_ub, r_ub=r_ub,
            model=result.model_name,
        )

    def _parse_dihedral(self, w: list[str], result: ParseResult) -> None:
        if len(w) < 7:
            return
        types = w[:4]
        if tuple(reversed(types)) < tuple(types):
            types = list(reversed(types))
        key = "-".join(types)
        kchi, n, delta = float(w[4]), int(float(w[5])), float(w[6]) % 360.0
        if abs(delta) < 1e-6:
            d = 1
        elif abs(delta - 180.0) < 1e-6:
            d = -1
        else:
            raise ValueError(f"{key}: dihedral phase {delta} (only 0 and 180 handled)")
        term = CHARMMDihedralParams(k=kchi, d=d, n=n, w=1.0, model=result.model_name)
        previous = result.dihedrals.get(key)
        if previous is None:
            result.dihedrals[key] = term
            return
        terms = previous.terms if isinstance(previous, CHARMMDihedralSeries) else [previous]
        # CHARMM: a repeated multiplicity replaces the earlier one
        terms = [t for t in terms if t.n != n] + [term]
        # 1-4 weight on the first term only (the pair must be counted once)
        terms = [
            CHARMMDihedralParams(k=t.k, d=t.d, n=t.n, w=1.0 if i == 0 else 0.0,
                                 model=t.model)
            for i, t in enumerate(sorted(terms, key=lambda t: t.n))
        ]
        result.dihedrals[key] = CHARMMDihedralSeries(terms=terms, model=result.model_name)

    def _parse_improper(self, w: list[str], result: ParseResult) -> None:
        if len(w) < 7:
            return
        if int(float(w[5])) != 0:
            raise ValueError(f"{'-'.join(w[:4])}: periodic improper not handled")
        result.impropers["-".join(w[:4])] = HarmonicImproperParams(
            k=float(w[4]), chi0=float(w[6]), model=result.model_name
        )

    # -- nonbonded ---------------------------------------------------------

    def _parse_nonbonded(self, w: list[str], result: ParseResult) -> None:
        if len(w) < 4 or not _is_float(w[2]):
            return
        atom = w[0]
        lj = LJParams(
            epsilon=-float(w[2]),
            sigma=2 * float(w[3]) / _SIXTH_ROOT_2,
            model=result.model_name,
        )
        if len(w) >= 7 and _is_float(w[5]) and _is_float(w[6]):
            lj.epsilon14 = -float(w[5])
            lj.sigma14 = 2 * float(w[6]) / _SIXTH_ROOT_2
        result.lj[f"{atom}-{atom}"] = lj
        if atom not in result.atoms:
            result.atoms[atom] = AtomType(
                element="X", mass=None, environment=f"CHARMM type {atom}",
                model=result.model_name,
            )

    def _parse_nbfix(self, w: list[str], result: ParseResult) -> None:
        if len(w) < 4 or not _is_float(w[2]):
            return
        i, j = w[0], w[1]
        key = f"{i}-{j}" if i <= j else f"{j}-{i}"
        result.lj[key] = LJParams(
            epsilon=-float(w[2]),
            sigma=float(w[3]) / _SIXTH_ROOT_2,
            ref="NBFIX",
            model=result.model_name,
        )


def _is_float(token: str) -> bool:
    try:
        float(token)
    except ValueError:
        return False
    return True


class CHARMMLoader:
    """Load a CHARMM .prm/.str file into a `ForceFieldDatabase`."""

    def __init__(self, db, metadata: dict | None = None):
        self.db = db
        self.metadata = metadata or {}

    def load_from_file(self, filepath, model_name: str) -> None:
        result = CHARMMParser(model_name, self.metadata).parse_file(str(filepath))
        # Types known only from NONBONDED (no MASS line) are dropped: they
        # cannot be assigned without a mass.
        for name in [t for t, a in result.atoms.items() if a.mass is None]:
            del result.atoms[name]
        result.to_database(self.db)
