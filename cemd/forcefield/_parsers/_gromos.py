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


import re

from ..._constants import two_letter_element_from_type
from ..forcefield_database import ForceFieldDatabase
from ..models import (
    AtomType,
    FourierDihedralParams,
    FourierTerm,
    HarmonicAngleParams,
    HarmonicBondParams,
    HarmonicImproperParams,
    LJParams,
)
from ._base import BaseForceFieldParser, ParseResult


class GromosLTParser(BaseForceFieldParser):
    """
    Parser for GROMOS files in moltemplate/LAMMPS (.lt) format.
    """

    def __init__(
        self, model_name: str = "gromos", display_name: str = "GROMOS 54A7 ATB"
    ):
        super().__init__()
        self.model_name = model_name  # Key in database
        self.display_name = display_name  # Name displayed in template

    def parse(self, content: str) -> ParseResult:
        """Parse an LT string."""
        lines = self._clean_content(content)

        result = ParseResult(model_name=self.model_name)
        result.metadata = {
            "name": self.display_name,  # Name displayed in template
            "description": "GROMOS 54A7 force field with ATB modifications",
            "ref": "https://doi.org/10.1007/s00249-011-0700-9",
            "tags": ["GROMOS", "54A7", "ATB", "moltemplate"],
        }

        # Parse each section
        self._parse_masses(lines, result)
        self._parse_pair_coeff(lines, result)
        self._parse_bond_coeff(lines, result)
        self._parse_angle_coeff(lines, result)
        self._parse_dihedral_coeff(lines, result)
        self._parse_improper_coeff(lines, result)

        return result

    def parse_file(self, filepath: str) -> ParseResult:
        """Parse an LT file."""
        with open(filepath) as f:
            content = f.read()
        return self.parse(content)

    def _strip_comments(self, line: str) -> str:
        """Remove comments # from LT format."""
        if "#" in line:
            return line[: line.index("#")]
        return line

    def _parse_masses(self, lines: list[str], result: ParseResult) -> None:
        "Parse atomic masses."
        pattern = r"mass\s+@atom:(\w+)\s+([\d.]+)"

        for line in lines:
            match = re.search(pattern, line)
            if match:
                atom_type, mass = match.groups()
                result.atoms[atom_type] = AtomType(
                    element=self._guess_element(atom_type, float(mass)),
                    # GROMOS carries partial charges per atom in the molecule
                    # topology, not per atom type.
                    charge=None,
                    mass=float(mass),
                    environment="",
                    model=result.model_name,
                )

    def _guess_element(self, atom_type: str, mass: float) -> str:
        "Guess the element from the name of the atom type."
        # A two-letter symbol whose known mass agrees with the parsed mass
        # is the element: this reads "CL" as chlorine and "SI" as silicon
        # while leaving one-letter-plus-digit united atoms like "CH2" or
        # "NR1" -- whose first two letters aren't a real element's mass --
        # to the special-cased and first-character fallbacks below.
        two_letter = two_letter_element_from_type(atom_type, mass)
        if two_letter is not None:
            return two_letter

        # Manage special cases
        special = {
            "OM": "O",
            "OA": "O",
            "OE": "O",
            "OW": "O",
            "NT": "N",
            "NL": "N",
            "NR": "N",
            "NZ": "N",
            "NE": "N",
            "CH0": "C",
            "CH1": "C",
            "CH2": "C",
            "CH3": "C",
            "CH4": "C",
            "CH2r": "C",
            "CR1": "C",
            "HC": "H",
            "HS14": "H",
            "CU1+": "Cu",
            "CU2+": "Cu",
            "ZN2+": "Zn",
            "MG2+": "Mg",
            "CA2+": "Ca",
            "NA+": "Na",
            "CL-": "Cl",
            "CLAro": "Cl",
            "CLOpt": "Cl",
            "CLChl": "Cl",
            "CLCl4": "Cl",
            "BROpt": "Br",
        }

        if atom_type in special:
            return special[atom_type]

        # By default, take the first character
        return atom_type[0]

    def _parse_pair_coeff(self, lines: list[str], result: ParseResult) -> None:
        """Parse LJ parameters."""
        pattern = r"pair_coeff\s+@atom:(\w+)\s+@atom:(\w+)\s+([\d.E+-]+)\s+([\d.E+-]+)"

        for line in lines:
            match = re.search(pattern, line)
            if match:
                type1, type2, epsilon, sigma = match.groups()

                # Canonical key (alphabetical order)
                key = f"{type1}-{type2}" if type1 <= type2 else f"{type2}-{type1}"

                # Ignore zero settings (polar hydrogens)
                if float(epsilon) == 0.0 and float(sigma) == 0.0:
                    continue

                result.lj[key] = LJParams(
                    epsilon=float(epsilon), sigma=float(sigma), model=result.model_name
                )

    def _parse_bond_coeff(self, lines: list[str], result: ParseResult) -> None:
        """Parses binding parameters."""
        # Both columns are written in scientific notation, so the exponent
        # has to be part of the second group too -- see `_parse_angle_coeff`.
        pattern = r"bond_coeff\s+@bond:(\w+)\s+([\d.E+-]+)\s+([\d.E+-]+)"

        for line in lines:
            match = re.search(pattern, line)
            if match:
                bond_id, k, r0 = match.groups()
                result.bonds[bond_id] = HarmonicBondParams(
                    k=float(k), r0=float(r0), model=result.model_name
                )

    def _parse_angle_coeff(self, lines: list[str], result: ParseResult) -> None:
        """Parse angle parameters."""
        # The equilibrium angle is written as e.g. "1.2600000000E+02" for
        # 126 deg. Leaving "E+-" out of this group truncated the match at
        # the mantissa, so every GROMOS angle came out a factor of ten or a
        # hundred too small (126 deg read as 1.26, 109.5 as 1.095).
        pattern = r"angle_coeff\s+@angle:(\w+)\s+([\d.E+-]+)\s+([\d.E+-]+)"

        for line in lines:
            match = re.search(pattern, line)
            if match:
                angle_id, k, theta0 = match.groups()
                result.angles[angle_id] = HarmonicAngleParams(
                    k=float(k), theta0=float(theta0), model=result.model_name
                )

    @staticmethod
    def _parse_dihedral_coeff(lines: list[str], result: ParseResult) -> None:
        """Parse dihedral parameters supporting multiple Fourier terms."""

        pattern = r"dihedral_coeff\s+@dihedral:(\w+)\s+(.+)"

        for line in lines:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            match = re.search(pattern, line)
            if not match:
                continue

            dihedral_id, rest = match.groups()

            coeffs = rest.split()
            if not coeffs:
                continue

            try:
                m_terms = int(coeffs[0])
                expected_len = 1 + 3 * m_terms

                if len(coeffs) < expected_len:
                    continue

                terms_list = []
                idx = 1
                for _ in range(m_terms):
                    k_val = float(coeffs[idx])
                    n_val = int(
                        float(coeffs[idx + 1])
                    )  # int(float(...)) handles values like "1.0"
                    delta_val = float(coeffs[idx + 2])

                    terms_list.append(FourierTerm(k=k_val, n=n_val, delta=delta_val))
                    idx += 3

                result.dihedrals[dihedral_id] = FourierDihedralParams(terms=terms_list)

            except (ValueError, TypeError):
                continue

    def _parse_improper_coeff(self, lines: list[str], result: ParseResult) -> None:
        """Parse the parameters as improper."""
        pattern = r"improper_coeff\s+@improper:(\w+)\s+([\d.E+-]+)\s+([\d.E+-]+)"

        for line in lines:
            match = re.search(pattern, line)
            if match:
                improper_id, k, chi0 = match.groups()
                result.impropers[improper_id] = HarmonicImproperParams(
                    k=float(k), chi0=float(chi0), model=result.model_name
                )


class GromosForceFieldLoader:
    """
    Loader to integrate GROMOS settings into ForceFieldDatabase.
    """

    def __init__(self, db: "ForceFieldDatabase"):
        self.db = db
        self.parser = GromosLTParser()

    def load_from_file(self, filepath: str, model_name: str = None) -> None:
        """
        Loads a LT GROMOS file into the database.

        Settings
        ----------
        filepath: str
            Path to LT file
        model_name: str, optional
            Key in database (default: 'gromos')
        """
        if model_name:
            self.parser.model_name = model_name
        # Don't overwrite the display_name
        result = self.parser.parse_file(filepath)
        result.to_database(self.db)

    def load_from_string(self, content: str, model_name: str = None) -> None:
        """
        Loads a LT GROMOS string into the database.

        Settings
        ----------
        content: str
            LT content
        model_name: str, optional
            Key in database (default: 'gromos')
        """
        if model_name:
            self.parser.model_name = model_name
        # Don't overwrite the display_name
        result = self.parser.parse(content)
        result.to_database(self.db)
