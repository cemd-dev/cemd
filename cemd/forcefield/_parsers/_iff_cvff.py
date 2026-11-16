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

from ..forcefield_database import ForceFieldDatabase
from ..models import (
    AtomType,
    Class2AngleAngleParams,
    Class2AngleAngleTorsionParams,
    Class2BondAngleParams,
    Class2BondBondParams,
    CVFFImproperParams,
    HarmonicAngleParams,
    HarmonicBondParams,
    HarmonicDihedralParams,
    LJParams,
    MorseBondParams,
)
from ._base import BaseForceFieldParser, ParseResult


class CVFFInterfaceParser(BaseForceFieldParser):
    """
    Parser for the CVFF Interface file (cvff_interface_v1_5.txt).

    This parser extracts parameters from sections:
    -#atom_types: Definitions of atom types
    -#quadratic_bond: Harmonic bonds
    -#quadratic_angle: Harmonic angles
    -#torsion_1: Dihedrals
    -#out_of_plane: Improper
    -#nonbond(12-6): LJ parameters
    -#bond-bond: Class2 bond-bond terms
    -#bond-angle: Class2 terms bond-angle
    -#equivalence: types borrowing the parameters of another type, per
     kind of term (resolved at lookup time by ForceFieldDatabase)
    """

    def __init__(
        self, model_name: str = "iff_cvff", display_name: str = "CVFF INTERFACE FF 1.5"
    ):
        super().__init__()
        self.model_name = model_name
        self.display_name = display_name
        self._current_section = None
        self._line_number = 0

        # Patterns section for detection.
        # NOTE: order matters - "#out_of_plane-out_of_plane" also contains
        # "#out_of_plane" as a substring, so the more specific "oop-oop"
        # pattern must be checked first (see _detect_section).
        self._section_patterns = {
            "atom_types": ["#atom_types"],
            "equivalence": ["#equivalence"],
            "oop-oop": ["#out_of_plane-out_of_plane"],
            "quadratic_bond": ["#quadratic_bond"],
            "quadratic_angle": ["#quadratic_angle"],
            "torsion_1": ["#torsion_1"],
            "out_of_plane": ["#out_of_plane"],
            "nonbond": ["#nonbond(12-6)"],
            "bond-bond": ["#bond-bond"],
            "bond-angle": ["#bond-angle"],
            "angle-angle-torsion": ["#angle-angle-torsion_1"],
            "angle-angle": ["#angle-angle_1"],
            "morse_bond": ["#morse_bond"],
            "improper": ["#improper"],
        }

    def parse(self, content: str) -> ParseResult:
        """Parses a CVFF Interface string."""
        lines = self._clean_content(content)

        result = ParseResult(model_name=self.model_name)
        result.metadata = {
            "description": "CVFF Interface force field for clay minerals, silicates, and cements",
            "ref": " https://doi.org/10.1021/la3038846",
            "tags": ["CVFF", "Interface", "Clay", "Silicate", "Cement"],
        }

        self._parse_sections(lines, result)
        return result

    def parse_file(self, filepath: str) -> ParseResult:
        """Parse a CVFF Interface file."""
        with open(filepath) as f:
            content = f.read()
        return self.parse(content)

    def _strip_comments(self, line: str) -> str:
        """Remove comments! from CVFF format."""
        if "!" in line:
            return line[: line.index("!")]
        return line

    def _detect_section(self, line: str) -> str | None:
        """Detect the section name from the line."""
        line_lower = line.lower()

        # Skip "_auto" tables: these hold wildcard/equivalence-based
        # parameters generated automatically by msi2lmp, as opposed to the
        # explicit "cvff" tables. The explicit "#equivalence" table is read
        # and resolved at lookup time; the "_auto" one (and the generic
        # parameters it points to) is intentionally ignored, so these
        # duplicate sections (e.g. "#quadratic_bond cvff_auto") do not
        # pollute the explicit ones.
        if "auto" in line_lower:
            return None

        for section_name, patterns in self._section_patterns.items():
            for pattern in patterns:
                if pattern.lower() in line_lower:
                    return section_name
        return None

    def _parse_sections(self, lines: list[str], result: ParseResult) -> None:
        """Parses the different sections of the CVFF file."""
        i = 0
        self._current_section = None

        while i < len(lines):
            line = lines[i]
            self._line_number = i + 1

            # Detect a section change
            if line.startswith("#"):
                section = self._detect_section(line)
                if section:
                    self._current_section = section
                else:
                    self._current_section = None
                i += 1
                continue

            # Ignore comments and empty lines
            if not line or line.startswith("!") or line.startswith(">"):
                i += 1
                continue

            # Parse the line if you are in an active section
            if self._current_section is not None and line.strip():
                self._parse_line(line, self._current_section, result)

            i += 1

    def _parse_line(self, line: str, section: str, result: ParseResult) -> None:
        """Parse a line according to the section."""
        # Handlers are defined as private methods
        handlers = {
            "atom_types": self._parse_atom_type,
            "equivalence": self._parse_equivalence,
            "quadratic_bond": self._parse_bond,
            "quadratic_angle": self._parse_angle,
            "torsion_1": self._parse_torsion,
            "out_of_plane": self._parse_improper,
            "improper": self._parse_improper,
            "nonbond": self._parse_nonbond,
            "bond-bond": self._parse_bondbond,
            "bond-angle": self._parse_bondangle,
            "angle-angle-torsion": self._parse_angleangletorsion,
            "angle-angle": self._parse_angleangle,
            "morse_bond": self._parse_morse_bond,
            # "oop-oop" (out_of_plane-out_of_plane) has no handler: this
            # class2 cross term has no corresponding dataclass in models.py,
            # so its section is recognized (to avoid misparsing) but skipped.
        }

        handler = handlers.get(section)
        if handler:
            handler(line, result)

    def _parse_atom_type(self, line: str, result: ParseResult) -> None:
        """
        Parse a line of atom type.
        Format: Ver Ref Type Mass Element Connections Comment
        """
        parts = line.split()
        if len(parts) >= 6:
            atom_type = parts[2]
            try:
                mass = float(parts[3])
            except ValueError:
                mass = 0.0
            element = parts[4]

            comment = " ".join(parts[6:]) if len(parts) > 6 else ""

            # The .frc file carries no charge: IFF documents it in the
            # comment, and only some types have a single value. Others depend
            # on the environment ("-0.55, -0.783 if next to Al defect",
            # "+1.1 in bulk SiO2 or SiOH, +0.725 in SiO-") or name several
            # ("(+1.45) or magnesium defect (+1.1)"). Those get None (no
            # per-type charge), so that typing rules or the topology set
            # them, rather than a wrong 0.0 or the first value of a list.
            charge = None
            groups = re.findall(r"\(([^()]*)\)", comment)
            numbers = [
                g for g in groups if re.fullmatch(r"\s*[+-]?\d+\.?\d*\s*", g)
            ]
            if len(groups) == 1 and len(numbers) == 1:
                charge = float(numbers[0])

            result.atoms[atom_type] = AtomType(
                element=element,
                charge=charge,
                mass=mass,
                environment=comment,
                model=result.model_name,
            )

    def _parse_equivalence(self, line: str, result: ParseResult) -> None:
        """
        Parse an equivalence line.
        Format: Ver Ref Type NonB Bond Angle Torsion OOP
        """
        parts = line.split()
        if len(parts) < 8:
            return
        atom_type = parts[2]
        columns = ("nonbond", "bond", "angle", "torsion", "oop")
        mapping = dict(zip(columns, parts[3:8]))
        if any(v != atom_type for v in mapping.values()):
            result.equivalence[atom_type] = mapping

    def _parse_bond(self, line: str, result: ParseResult) -> None:
        """
        Parses a harmonic bond line.
        Format: Ver Ref I J R0 K
        """
        parts = line.split()
        if len(parts) >= 6:
            i = parts[2]
            j = parts[3]
            try:
                r0 = float(parts[4])
                k = float(parts[5])
                key = f"{i}-{j}" if i <= j else f"{j}-{i}"
                result.bonds[key] = HarmonicBondParams(
                    k=k, r0=r0, model=result.model_name
                )
            except ValueError:
                pass

    def _parse_angle(self, line: str, result: ParseResult) -> None:
        """
        Parse a harmonic angle line.
        Format: Ver Ref I J K Theta0 K2
        """
        parts = line.split()
        if len(parts) >= 7:
            i = parts[2]
            j = parts[3]
            k = parts[4]
            try:
                theta0 = float(parts[5])
                k2 = float(parts[6])
                key = f"{i}-{j}-{k}"
                result.angles[key] = HarmonicAngleParams(
                    k=k2, theta0=theta0, model=result.model_name
                )
            except ValueError:
                pass

    @staticmethod
    def _phase_to_d(phase: float) -> int:
        """
        Convert a CVFF phase angle (Phi0/Chi0, in degrees) to the +1/-1
        'd' parameter expected by LAMMPS periodic potentials.

        Raises
        ------
        ValueError
            If phase is neither 0 nor 180 degrees (unsupported by LAMMPS).
        """
        if abs(phase) < 1e-5:
            return 1
        if abs(phase - 180.0) < 1e-5:
            return -1
        raise ValueError(f"Valeur de phase non supportée par LAMMPS: {phase}")

    def _parse_torsion(self, line: str, result: ParseResult) -> None:
        """
        Parse a twist line from CVFF and convert to LAMMPS harmonic.
        Format CVFF: Ver Ref I J K L Kphi n Phi0
        """
        parts = line.split()
        if len(parts) >= 9:
            i = parts[2]
            j = parts[3]
            k = parts[4]
            l_ = parts[5]

            try:
                kphi = float(parts[6])
                n = int(float(parts[7]))
                phi0 = float(parts[8])
                d_val = self._phase_to_d(phi0)

                key = f"{i}-{j}-{k}-{l_}"

                result.dihedrals[key] = HarmonicDihedralParams(
                    k=kphi, n=n, d=d_val, ref=parts[1], model=result.model_name
                )
            except ValueError:
                pass

    def _parse_improper(self, line: str, result: ParseResult) -> None:
        """
        Parse an out-of-plane (improper) line.
        Format: Ver Ref I J K L Kchi n Chi0
        CVFF out-of-plane potential: E = Kchi * [1 + cos(n*Chi - Chi0)],
        i.e. the periodic form (LAMMPS improper_style cvff), not harmonic.
        """
        parts = line.split()
        if len(parts) >= 9:
            i = parts[2]
            j = parts[3]
            k = parts[4]
            l_ = parts[5]
            try:
                kchi = float(parts[6])
                n = int(float(parts[7]))
                chi0 = float(parts[8])
                d_val = self._phase_to_d(chi0)

                key = f"{i}-{j}-{k}-{l_}"
                result.impropers[key] = CVFFImproperParams(
                    k=kchi, d=d_val, n=n, model=result.model_name
                )
            except ValueError:
                pass

    def _parse_nonbond(self, line: str, result: ParseResult) -> None:
        """
        Parses a line of LJ parameters.
        Format: Ver Ref I A B
        """
        parts = line.split()
        if len(parts) >= 5:
            atom = parts[2]
            try:
                a = float(parts[3])
                b = float(parts[4])
                key = f"{atom}-{atom}"
                if a > 0 and b > 0:
                    epsilon = b**2 / (4 * a)
                    sigma = (a / b) ** (1 / 6)
                    result.lj[key] = LJParams(
                        epsilon=epsilon, sigma=sigma, model=result.model_name
                    )
                elif b == 0:
                    # Hydrogens without LJ (B = 0, A a token 1e-8): keep them,
                    # with epsilon 0, so they are not reported as missing.
                    result.lj[key] = LJParams(
                        epsilon=0.0, sigma=0.0, model=result.model_name
                    )
            except ValueError:
                pass

    def _parse_bondbond(self, line: str, result: ParseResult) -> None:
        """
        Parses a bond-bond line (class2).
        Format: Ver Ref I J K M
        """
        parts = line.split()
        if len(parts) >= 6:
            i = parts[2]
            j = parts[3]
            k = parts[4]
            try:
                m = float(parts[5])
                key = f"{i}-{j}-{k}"
                result.bondbond[key] = Class2BondBondParams(
                    m=m, r1=0.0, r2=0.0, model=result.model_name
                )
            except ValueError:
                pass

    def _parse_bondangle(self, line: str, result: ParseResult) -> None:
        """
        Parse a bond-angle line (class2).
        Format: Ver Ref I J K N1 [N2]
        N2 is only given when the angle is asymmetric (I != K); for a
        symmetric angle (I == K) a single constant applies to both bonds.
        """
        parts = line.split()
        if len(parts) >= 6:
            i = parts[2]
            j = parts[3]
            k = parts[4]
            try:
                n1 = float(parts[5])
                n2 = float(parts[6]) if len(parts) > 6 else n1
                key = f"{i}-{j}-{k}"
                result.bondangle[key] = Class2BondAngleParams(
                    n1=n1, n2=n2, r1=0.0, r2=0.0, model=result.model_name
                )
            except ValueError:
                pass

    def _parse_morse_bond(self, line: str, result: ParseResult) -> None:
        """
        Parse a line of Morse bond.
        Format: Ver Ref I J R0 D ALPHA
        Only used as a fallback for pairs not already covered by the
        quadratic_bond (harmonic) table, which is the default bond style.
        """
        parts = line.split()
        if len(parts) >= 7:
            i = parts[2]
            j = parts[3]
            try:
                r0 = float(parts[4])
                d = float(parts[5])
                alpha = float(parts[6])
                key = f"{i}-{j}"
                if key not in result.bonds:
                    result.bonds[key] = MorseBondParams(
                        r0=r0, D=d, alpha=alpha, model=result.model_name
                    )
            except ValueError:
                pass

    def _parse_angleangletorsion(self, line: str, result: ParseResult) -> None:
        """
        Parse an angle-angle-torsion cross term line (class2).
        Format: Ver Ref I J K L K(Ang,Ang,Tor)
        The reference angles theta1/theta2 are not given in this table (they
        come from the corresponding quadratic_angle entries) and are left
        at 0.0.
        """
        parts = line.split()
        if len(parts) >= 7:
            i = parts[2]
            j = parts[3]
            k = parts[4]
            l_ = parts[5]
            try:
                m = float(parts[6])
                key = f"{i}-{j}-{k}-{l_}"
                result.angleangletorsion[key] = Class2AngleAngleTorsionParams(
                    m=m, theta1=0.0, theta2=0.0, model=result.model_name
                )
            except ValueError:
                pass

    def _parse_angleangle(self, line: str, result: ParseResult) -> None:
        """
        Parse an angle-angle (improper) cross term line (class2).
        Format: Ver Ref I J K L K(Ang,Ang)
        A single constant is given; msi2lmp uses it for M1=M2=M3. The
        reference angles theta1/theta2/theta3 are not given in this table
        (they come from the corresponding quadratic_angle entries) and are
        left at 0.0.
        """
        parts = line.split()
        if len(parts) >= 7:
            i = parts[2]
            j = parts[3]
            k = parts[4]
            l_ = parts[5]
            try:
                m = float(parts[6])
                key = f"{i}-{j}-{k}-{l_}"
                result.angleangle[key] = Class2AngleAngleParams(
                    m1=m,
                    m2=m,
                    m3=m,
                    theta1=0.0,
                    theta2=0.0,
                    theta3=0.0,
                    model=result.model_name,
                )
            except ValueError:
                pass


class CVFFInterfaceLoader:
    """Loader to integrate CVFF Interface parameters into ForceFieldDatabase."""

    def __init__(self, db: "ForceFieldDatabase"):
        self.db = db
        self.parser = CVFFInterfaceParser()

    def load_from_file(self, filepath: str, model_name: str = None) -> None:
        """
        Load from file.

        Parameters
        ----------
        filepath : str
            Path to the CVFF file
        model_name : str, optional
            Clé dans la base de données (par défaut: 'iff_cvff')
        """
        if model_name:
            self.parser.model_name = model_name
        result = self.parser.parse_file(filepath)
        result.to_database(self.db)

    def load_from_string(self, content: str, model_name: str = None) -> None:
        """
        Load from string.

        Parameters
        ----------
        content : str
            CVFF content
        model_name : str, optional
            Clé dans la base de données (par défaut: 'iff_cvff')
        """
        if model_name:
            self.parser.model_name = model_name
        result = self.parser.parse(content)
        result.to_database(self.db)
