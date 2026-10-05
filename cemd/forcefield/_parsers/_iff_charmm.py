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

from typing import Any

from ..._constants import two_letter_element_from_type
from ..forcefield_database import ForceFieldDatabase
from ..models import (
    AtomType,
    CHARMMDihedralParams,
    HarmonicAngleParams,
    HarmonicBondParams,
    HarmonicImproperParams,
    LJParams,
)
from ._base import BaseForceFieldParser, ParseResult


class CHARMMInterfaceParser(BaseForceFieldParser):
    """
    Parser for the CHARMM Interface file (charmm27_interface_v1_5.txt).

    This parser extracts parameters from sections:
    - BONDS: Harmonic bonds
    - ANGLES: Harmonic angles (with Urey-Bradley)
    - DIHEDRALS: Dihedrals
    - IMPROPER: Impropers
    - NONBONDED: LJ parameters (with 1-4 scaling)
    - CMAP: CMAP correction grids (stored as metadata)

    Atom types are extracted from the NONBONDED section.
    """

    def __init__(
        self,
        model_name: str = "iff_charmm",
        display_name: str = "CHARMM INTERFACE FF 1.5",
    ):
        super().__init__()
        self.model_name = model_name  # Key in database
        self.display_name = display_name  # Name displayed in model
        self._current_section = None
        self._line_number = 0

        # Section patterns for detection
        self._section_patterns = {
            "bonds": ["BONDS"],
            "angles": ["ANGLES"],
            "dihedrals": ["DIHEDRALS"],
            "improper": ["IMPROPER"],
            "nonbonded": ["NONBONDED"],
            "cmap": ["CMAP"],
        }

        # Store CMAP data for metadata
        self._cmap_data: list[dict[str, Any]] = []

        # Collect atom types from nonbonded section
        self._atom_types: set[str] = set()

    def parse(self, content: str) -> ParseResult:
        """Parses a CHARMM Interface string."""
        lines = self._clean_content(content)

        result = ParseResult(model_name=self.model_name)
        result.metadata = {
            "name": self.display_name,  # Name displayed in model
            "description": "CHARMM27 Interface force field for clay minerals, silicates, cements, and metals",
            "ref": " https://doi.org/10.1021/la3038846",
            "tags": ["CHARMM", "Interface", "Clay", "Silicate", "Cement", "Metal"],
        }

        self._parse_sections(lines, result)

        # Add CMAP data to metadata if present
        if self._cmap_data:
            result.metadata["cmap"] = self._cmap_data

        return result

    def parse_file(self, filepath: str) -> ParseResult:
        """Parse a CHARMM Interface file."""
        with open(filepath) as f:
            content = f.read()
        return self.parse(content)

    def _strip_comments(self, line: str) -> str:
        """Remove comments from CHARMM format."""
        if "!" in line:
            return line[: line.index("!")]
        return line

    def _detect_section(self, line: str) -> str | None:
        """Detect the section name from the line."""
        line_upper = line.upper()
        for section_name, patterns in self._section_patterns.items():
            for pattern in patterns:
                if pattern.upper() in line_upper:
                    return section_name
        return None

    def _parse_sections(self, lines: list[str], result: ParseResult) -> None:
        """Parses the different sections of the CHARMM file."""
        i = 0
        self._current_section = None

        while i < len(lines):
            line = lines[i]
            self._line_number = i + 1

            # Detect a section change
            if line.startswith("*") or line.startswith("!"):
                # Skip header/comments
                i += 1
                continue

            # Check for section headers
            section = self._detect_section(line)
            if section:
                self._current_section = section
                i += 1
                continue

            # Ignore empty lines
            if not line:
                i += 1
                continue

            # Parse the line if you are in an active section
            if self._current_section is not None and line.strip():
                self._parse_line(line, self._current_section, result)

            i += 1

    def _parse_line(self, line: str, section: str, result: ParseResult) -> None:
        """Parse a line according to the section."""
        handlers = {
            "bonds": self._parse_bond,
            "angles": self._parse_angle,
            "dihedrals": self._parse_dihedral,
            "improper": self._parse_improper,
            "nonbonded": self._parse_nonbonded,
            "cmap": self._parse_cmap,
        }

        handler = handlers.get(section)
        if handler:
            handler(line, result)

    def _parse_bond(self, line: str, result: ParseResult) -> None:
        """
        Parse a bond line.
        Format: i j Kb b0
        Example: CST OST 937.96 1.1600
        """
        parts = line.split()
        # Skip comment lines starting with !
        if not parts or parts[0].startswith("!"):
            return

        if len(parts) >= 4:
            i = parts[0]
            j = parts[1]
            try:
                k = float(parts[2])
                r0 = float(parts[3])
                key = f"{i}-{j}" if i <= j else f"{j}-{i}"
                result.bonds[key] = HarmonicBondParams(
                    k=k, r0=r0, model=result.model_name
                )
                # Collect atom types
                self._atom_types.add(i)
                self._atom_types.add(j)
            except ValueError:
                pass

    def _parse_angle(self, line: str, result: ParseResult) -> None:
        """
        Parse an angle line.
        Format: i j k Ktheta Theta0 [Kub S0]
        Example: CA CA CA 40.000 120.00 35.00 2.41620
        """
        parts = line.split()
        if not parts or parts[0].startswith("!"):
            return

        if len(parts) >= 5:
            i = parts[0]
            j = parts[1]
            k = parts[2]
            try:
                ktheta = float(parts[3])
                theta0 = float(parts[4])
                key = f"{i}-{j}-{k}"
                result.angles[key] = HarmonicAngleParams(
                    k=ktheta, theta0=theta0, model=result.model_name
                )
                # Collect atom types
                self._atom_types.add(i)
                self._atom_types.add(j)
                self._atom_types.add(k)
                # Note: Urey-Bradley parameters (Kub, S0), present on some
                # lines as trailing columns, are not stored as they are not
                # supported by the current angle model.
            except ValueError:
                pass

    def _parse_dihedral(self, line: str, result: ParseResult) -> None:
        """
        Parse a dihedral line.
        Format: i j k l Kchi n delta
        Example: CA CA CA CA 3.1000 2 180.00
        """
        parts = line.split()
        if not parts or parts[0].startswith("!"):
            return

        if len(parts) >= 7:
            i = parts[0]
            j = parts[1]
            k = parts[2]
            l_ = parts[3]
            try:
                kchi = float(parts[4])
                n = int(float(parts[5]))
                delta = float(parts[6])

                if abs(delta) < 1e-5:
                    d_val = 1
                elif abs(delta - 180.0) < 1e-5:
                    d_val = -1
                else:
                    raise ValueError(
                        f"Valeur delta non supportée par LAMMPS dihedral charmm: {delta}"
                    )

                key = f"{i}-{j}-{k}-{l_}"
                # w=1.0: the Interface FF NONBONDED section defines dedicated
                # 1-4 LJ parameters (e14fac=1.0), so 1-4 interactions are kept.
                result.dihedrals[key] = CHARMMDihedralParams(
                    k=kchi, d=d_val, n=n, w=1.0, model=result.model_name
                )
                # Collect atom types
                self._atom_types.add(i)
                self._atom_types.add(j)
                self._atom_types.add(k)
                self._atom_types.add(l_)
            except ValueError:
                pass

    def _parse_improper(self, line: str, result: ParseResult) -> None:
        """
        Parse an improper line.
        Format: i j k l Kpsi ignored psi0
        Example: CPB CPA NPH CPA 20.8000 0 0.0000
        The second-to-last column is unused (kept for CHARMM compatibility).
        """
        parts = line.split()
        if not parts or parts[0].startswith("!"):
            return

        if len(parts) >= 7:
            i = parts[0]
            j = parts[1]
            k = parts[2]
            l_ = parts[3]
            try:
                kpsi = float(parts[4])
                psi0 = float(parts[6])
                key = f"{i}-{j}-{k}-{l_}"
                result.impropers[key] = HarmonicImproperParams(
                    k=kpsi, chi0=psi0, model=result.model_name
                )
                # Collect atom types
                self._atom_types.add(i)
                self._atom_types.add(j)
                self._atom_types.add(k)
                self._atom_types.add(l_)
            except ValueError:
                pass

    def _parse_nonbonded(self, line: str, result: ParseResult) -> None:
        """
        Parse a nonbonded LJ line.
        Format: atom ignored epsilon Rmin/2 ignored eps,1-4 Rmin/2,1-4
        Example: C 0.000000 -0.110000 2.000000
        """
        parts = line.split()
        if not parts or parts[0].startswith("!"):
            return

        # Skip the initial "NONBONDED" header line
        if "NONBONDED" in line.upper():
            return

        if len(parts) >= 4:
            atom = parts[0]
            # Skip comment lines
            if atom.startswith("!"):
                return
            try:
                # CHARMM format: atom, ignored, epsilon, Rmin/2
                # epsilon is negative, Rmin/2 is positive
                epsilon = -float(parts[2])  # Convert to positive
                rmin_half = float(parts[3])
                # CHARMM Rmin = 2 * Rmin/2
                sigma = 2 * rmin_half / (2 ** (1 / 6))  # Convert to LJ sigma

                key = f"{atom}-{atom}"
                result.lj[key] = LJParams(
                    epsilon=epsilon, sigma=sigma, model=result.model_name
                )
                # Collect atom type
                self._atom_types.add(atom)
            except (ValueError, IndexError):
                pass

    def _parse_cmap(self, line: str, result: ParseResult) -> None:
        """
        Parse CMAP grid data.
        CMAP data is stored as metadata for reference.
        Format: C NH1 CT1 C NH1 CT1 C NH1 24
               followed by grid data
        """
        parts = line.split()
        if not parts:
            return

        if len(parts) >= 9 and parts[0] not in ["!", "*"]:
            try:
                # Check if this is a CMAP header
                if parts[0] in ["C", "N"] and len(parts) >= 8:
                    cmap_entry = {
                        "atoms": parts[:8],
                        "grid_size": int(parts[8]) if len(parts) > 8 else 24,
                        "data": [],
                    }
                    self._cmap_data.append(cmap_entry)
                    # Collect atom types from CMAP
                    for atom in parts[:8]:
                        self._atom_types.add(atom)
            except ValueError:
                pass
        elif self._cmap_data and len(parts) >= 24:
            # This is a grid row (24 values)
            try:
                row_data = [float(x) for x in parts]
                self._cmap_data[-1]["data"].append(row_data)
            except ValueError:
                pass

    def _guess_element(self, atom_type: str, mass: float) -> str:
        """Guess the element from atom type name."""
        # CHARMM atom types
        elements = {
            "C": "C",
            "CA": "C",
            "CC": "C",
            "CD": "C",
            "CE1": "C",
            "CE2": "C",
            "CM": "C",
            "CP1": "C",
            "CP2": "C",
            "CP3": "C",
            "CPA": "C",
            "CPB": "C",
            "CPH1": "C",
            "CPH2": "C",
            "CPM": "C",
            "CPT": "C",
            "CS": "C",
            "CST": "C",
            "CT1": "C",
            "CT2": "C",
            "CT3": "C",
            "CY": "C",
            "N": "N",
            "NC2": "N",
            "NH1": "N",
            "NH2": "N",
            "NH3": "N",
            "NP": "N",
            "NPH": "N",
            "NR1": "N",
            "NR2": "N",
            "NR3": "N",
            "NY": "N",
            "O": "O",
            "OB": "O",
            "OC": "O",
            "OH1": "O",
            "OM": "O",
            "OS": "O",
            "OST": "O",
            "OT": "O",
            "S": "S",
            "SM": "S",
            "SS": "S",
            "H": "H",
            "HA": "H",
            "HB": "H",
            "HC": "H",
            "HP": "H",
            "HR1": "H",
            "HR2": "H",
            "HR3": "H",
            "HS": "H",
            "HT": "H",
            "FE": "Fe",
            "ZN": "Zn",
            "NA+": "Na",
            "K+": "K",
            "CLA": "Cl",
            "CAL": "Ca",
            "MG": "Mg",
            "AG": "Ag",
            "AL": "Al",
            "AU": "Au",
            "CU": "Cu",
            "NI": "Ni",
            "PB": "Pb",
            "PD": "Pd",
            "PT": "Pt",
        }

        # Check for interface atoms
        interface_elements = {
            "SY1": "Si",
            "SY2": "Si",
            "AYT1": "Al",
            "AYT2": "Al",
            "AY1": "Al",
            "AY2": "Al",
            "OY1": "O",
            "OY2": "O",
            "OY3": "O",
            "OY4": "O",
            "OY5": "O",
            "OY6": "O",
            "OY7": "O",
            "OY8": "O",
            "OY9": "O",
            "HOY": "H",
            "HOK": "H",
            "SC1": "Si",
            "SC4": "Si",
            "OC1": "O",
            "OC2": "O",
            "OC23": "O",
            "OC24": "O",
            "OC3": "O",
            "OC4": "O",
            "OC5": "O",
            "AC1": "Al",
            "CA++": "Ca",
            "CA+A": "Ca",
            "CA+H": "Ca",
            "PAP": "P",
            "OAP1": "O",
            "OAP2": "O",
            "HOP": "H",
            "HOC": "H",
            "OCL": "O",
            "OHL": "O",
            "OSL": "O",
            "OBL": "O",
            "O2L": "O",
            "PL": "P",
            "SL": "S",
        }

        all_elements = {**elements, **interface_elements}

        # A two-letter symbol whose known mass agrees with the parsed mass
        # is the element: this reads "CL" as carbon (CHARMM's lipid tail
        # carbon) or "NA+" as sodium correctly, deferring to the
        # dictionaries above only when the mass doesn't back up that guess.
        two_letter = two_letter_element_from_type(atom_type, mass)
        if two_letter is not None:
            return two_letter

        # Some keys above, like "NA+" and "CA++", carry their +/- suffix;
        # look those up before stripping it.
        if atom_type in all_elements:
            return all_elements[atom_type]

        # Remove + and - suffixes
        clean_type = atom_type.rstrip("+-")

        if clean_type in all_elements:
            return all_elements[clean_type]

        # Try first character
        if clean_type and clean_type[0] in "CONHPS":
            return clean_type[0]

        return "X"

    def _guess_mass(self, atom_type: str) -> float:
        """Guess the mass from atom type name."""
        masses = {
            # Hydrogens
            "H": 1.008,
            "HA": 1.008,
            "HB": 1.008,
            "HC": 1.008,
            "HP": 1.008,
            "HR1": 1.008,
            "HR2": 1.008,
            "HR3": 1.008,
            "HS": 1.008,
            "HT": 1.008,
            "HOY": 1.008,
            "HOK": 1.008,
            "HOC": 1.008,
            "HOP": 1.008,
            "HOL": 1.008,
            "HAL1": 1.008,
            "HAL2": 1.008,
            "HAL3": 1.008,
            "HEL1": 1.008,
            "HEL2": 1.008,
            "HL": 1.008,
            "HCL": 1.008,
            # Carbons
            "C": 12.011,
            "CA": 12.011,
            "CC": 12.011,
            "CD": 12.011,
            "CE1": 12.011,
            "CE2": 12.011,
            "CM": 12.011,
            "CP1": 12.011,
            "CP2": 12.011,
            "CP3": 12.011,
            "CPA": 12.011,
            "CPB": 12.011,
            "CPH1": 12.011,
            "CPH2": 12.011,
            "CPM": 12.011,
            "CPT": 12.011,
            "CS": 12.011,
            "CST": 12.011,
            "CT1": 12.011,
            "CT2": 12.011,
            "CT3": 12.011,
            "CY": 12.011,
            "CT": 12.011,
            "CL": 12.011,
            "CTL1": 12.011,
            "CTL2": 12.011,
            "CTL3": 12.011,
            "CTL5": 12.011,
            "CEL1": 12.011,
            "CEL2": 12.011,
            "CF1": 12.011,
            "CF2": 12.011,
            "CF3": 12.011,
            "CAP": 12.011,
            "CN": 12.011,
            # Interface carbons
            "SY1": 28.086,
            "SY2": 28.086,
            "SC1": 28.086,
            "SC4": 28.086,
            # Nitrogens
            "N": 14.007,
            "NC2": 14.007,
            "NH1": 14.007,
            "NH2": 14.007,
            "NH3": 14.007,
            "NP": 14.007,
            "NPH": 14.007,
            "NR1": 14.007,
            "NR2": 14.007,
            "NR3": 14.007,
            "NY": 14.007,
            "NS1": 14.007,
            "NS2": 14.007,
            "NH3L": 14.007,
            "NTL": 14.007,
            "NC": 14.007,
            # Oxygens
            "O": 15.999,
            "OB": 15.999,
            "OC": 15.999,
            "OH1": 15.999,
            "OM": 15.999,
            "OS": 15.999,
            "OST": 15.999,
            "OT": 15.999,
            "OY1": 15.999,
            "OY2": 15.999,
            "OY3": 15.999,
            "OY4": 15.999,
            "OY5": 15.999,
            "OY6": 15.999,
            "OY7": 15.999,
            "OY8": 15.999,
            "OY9": 15.999,
            "OC1": 15.999,
            "OC2": 15.999,
            "OC23": 15.999,
            "OC24": 15.999,
            "OC3": 15.999,
            "OC4": 15.999,
            "OC5": 15.999,
            "OAP1": 15.999,
            "OAP2": 15.999,
            "OBL": 15.999,
            "OCL": 15.999,
            "OHL": 15.999,
            "OSL": 15.999,
            "O2L": 15.999,
            "OCA": 15.999,
            # Sulfurs
            "S": 32.065,
            "SM": 32.065,
            "SS": 32.065,
            "SL": 32.065,
            # Metals and others
            "FE": 55.845,
            "ZN": 65.380,
            "PAP": 30.974,
            "PL": 30.974,
            "CA++": 40.078,
            "CA+A": 40.078,
            "CA+H": 40.078,
            "CAL": 40.078,
            "MG": 24.305,
            "NA+": 22.990,
            "SOD": 22.990,
            "K+": 39.098,
            "POT": 39.098,
            "CLA": 35.453,
            "Cl": 35.453,
            "AG": 107.868,
            "AL": 26.982,
            "AU": 196.967,
            "CU": 63.546,
            "NI": 58.693,
            "PB": 207.200,
            "PD": 106.420,
            "PT": 195.084,
            "HE": 4.003,
            "NE": 20.180,
            "DUM": 0.000,
        }

        # Some keys above, like "NA+" and "CA++", carry their +/- suffix;
        # look those up before stripping it.
        if atom_type in masses:
            return masses[atom_type]

        clean_type = atom_type.rstrip("+-")
        if clean_type in masses:
            return masses[clean_type]
        return 0.0

    def finalize(self, result: ParseResult) -> None:
        """
        Add atom types to result from collected atom types.
        This should be called after parsing all sections.
        """
        # Sorted, because `_atom_types` is a set and Python randomizes
        # string hashing per process: iterating it directly made the parsed
        # order -- and therefore the generated docs/_static/ff_data.json --
        # differ on every run, for a 1500-line diff of pure noise.
        for atom_type in sorted(self._atom_types):
            # Skip wildcards
            if atom_type == "X" or atom_type == "*":
                continue

            # Skip if already present
            if atom_type in result.atoms:
                continue

            mass = self._guess_mass(atom_type)
            element = self._guess_element(atom_type, mass)

            result.atoms[atom_type] = AtomType(
                element=element,
                charge=None,  # Charges are not available in PAR file
                mass=mass,
                environment=f"CHARMM Interface atom type {atom_type} (charge not available)",
                model=result.model_name,
            )


class CHARMMInterfaceLoader:
    """Loader to integrate CHARMM Interface parameters into ForceFieldDatabase."""

    def __init__(self, db: "ForceFieldDatabase"):
        self.db = db
        self.parser = CHARMMInterfaceParser()

    def load_from_file(self, filepath: str, model_name: str = None) -> None:
        if model_name:
            self.parser.model_name = model_name
        result = self.parser.parse_file(filepath)
        # Finalize: add atom types from collected atoms
        self.parser.finalize(result)
        result.to_database(self.db)

    def load_from_string(self, content: str, model_name: str = None) -> None:
        if model_name:
            self.parser.model_name = model_name
        result = self.parser.parse(content)
        # Finalize: add atom types from collected atoms
        self.parser.finalize(result)
        result.to_database(self.db)
