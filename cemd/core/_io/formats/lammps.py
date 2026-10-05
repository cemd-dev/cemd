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

import datetime
from typing import Any, TextIO

import numpy as np
import pandas as pd

from ....forcefield.models import (
    BuckinghamParams,
    CHARMMDihedralParams,
    Class2AngleAngleParams,
    Class2AngleAngleTorsionParams,
    Class2AngleParams,
    Class2BondAngleParams,
    Class2BondBondParams,
    Class2BondParams,
    Class2DihedralParams,
    Class2ImproperParams,
    CVFFImproperParams,
    DistanceImproperParams,
    FourierDihedralParams,
    FourierTerm,
    HarmonicAngleParams,
    HarmonicBondParams,
    HarmonicDihedralParams,
    HarmonicImproperParams,
    LJParams,
    MorseBondParams,
)
from .base import BaseReader, BaseWriter


class LAMMPSReader(BaseReader):
    """Read LAMMPS data files."""

    HEADERS = {
        "Atom Type Labels",
        "Bond Type Labels",
        "Angle Type Labels",
        "Dihedral Type Labels",
        "Improper Type Labels",
        "Pair Coeffs",
        "PairIJ Coeffs",
        "Bond Coeffs",
        "Angle Coeffs",
        "Dihedral Coeffs",
        "Improper Coeffs",
        "BondBond Coeffs",
        "BondAngle Coeffs",
        "MiddleBondTorsion Coeffs",
        "EndBondTorsion Coeffs",
        "AngleTorsion Coeffs",
        "AngleAngleTorsion Coeffs",
        "AngleAngle Coeffs",
        "Masses",
        "Atoms",
        "Velocities",
        "Bonds",
        "Angles",
        "Dihedrals",
        "Impropers",
        "CS-Info",
    }

    @classmethod
    def read(cls, path: str) -> dict:
        """Read LAMMPS data file."""
        lines = cls._read_lines(path)
        sections = cls._parse_sections(lines)
        box = cls._parse_box(path)

        topology = {
            "box": box,
            "atoms": None,
            "bonds": None,
            "angles": None,
            "dihedrals": None,
            "impropers": None,
            "velocities": None,
            "masses": {},
            "charges": {},
            "bond_types": [],
            "angle_types": [],
            "dihedral_types": [],
            "improper_types": [],
            "atom_style": "full",
            "pair_params": {},
            "bond_params": {},
            "angle_params": {},
            "dihedral_params": {},
            "improper_params": {},
            "bondbond_params": {},
            "bondangle_params": {},
            "middlebondtorsion_params": {},
            "endbondtorsion_params": {},
            "angletorsion_params": {},
            "angleangletorsion_params": {},
            "angleangle_params": {},
        }

        cls._parse_section_data(sections, topology)
        return topology

    @staticmethod
    def _read_lines(path: str) -> list[str]:
        """Read and clean lines from file."""
        lines = []
        with open(path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    line = line.partition("#")[-1].strip()
                else:
                    line = line.partition("#")[0].strip()
                if line:
                    lines.append(line)
        return lines

    @staticmethod
    def _parse_sections(lines: list[str]) -> dict[str, list[str]]:
        """Parse sections from lines."""
        starts = [i for i, line in enumerate(lines) if line in LAMMPSReader.HEADERS]
        starts.append(None)
        return {
            lines[j]: lines[j + 1 : starts[i + 1]] for i, j in enumerate(starts[:-1])
        }

    @staticmethod
    def _parse_box(path: str) -> tuple:
        """Parse box parameters from file."""
        xlo = xhi = ylo = yhi = zlo = zhi = 0.0
        xy = xz = yz = 0.0

        with open(path, encoding="utf-8") as f:
            for line in f:
                stripped = line.strip()
                parts = stripped.split()
                if stripped.endswith("xlo xhi"):
                    xlo, xhi = float(parts[0]), float(parts[1])
                elif stripped.endswith("ylo yhi"):
                    ylo, yhi = float(parts[0]), float(parts[1])
                elif stripped.endswith("zlo zhi"):
                    zlo, zhi = float(parts[0]), float(parts[1])
                elif stripped.endswith("xy xz yz"):
                    xy, xz, yz = float(parts[0]), float(parts[1]), float(parts[2])

        return ((xlo, xhi), (ylo, yhi), (zlo, zhi), (xy, xz, yz))

    @classmethod
    def _parse_section_data(
        cls, sections: dict[str, list[str]], topology: dict[str, Any]
    ) -> None:
        """Redirects each section to the appropriate processing method."""

        type_labels = {
            "atom": {},
            "bond": {},
            "angle": {},
            "dihedral": {},
            "improper": {},
        }
        for key, values in sections.items():
            label_key = {
                "Atom Type Labels": "atom",
                "Bond Type Labels": "bond",
                "Angle Type Labels": "angle",
                "Dihedral Type Labels": "dihedral",
                "Improper Type Labels": "improper",
            }.get(key)
            if label_key is not None:
                type_labels[label_key] = {
                    int(parts[0]): parts[1]
                    for value in values
                    if len(parts := value.split()) >= 2
                }

                if label_key == "atom":
                    topology["atom_types"] = [
                        type_labels[label_key][type_id]
                        for type_id in sorted(type_labels[label_key])
                    ]
                else:
                    topology[f"{label_key}_types"] = [
                        type_labels[label_key][type_id]
                        for type_id in sorted(type_labels[label_key])
                    ]

            if key == "Bond Type Labels":
                topology["bond_types"] = [
                    type_labels["bond"][type_id]
                    for type_id in sorted(type_labels["bond"])
                ]
            elif key == "Angle Type Labels":
                topology["angle_types"] = [
                    type_labels["angle"][type_id]
                    for type_id in sorted(type_labels["angle"])
                ]
            elif key == "Dihedral Type Labels":
                topology["dihedral_types"] = [
                    type_labels["dihedral"][type_id]
                    for type_id in sorted(type_labels["dihedral"])
                ]
            elif key == "Improper Type Labels":
                topology["improper_types"] = [
                    type_labels["improper"][type_id]
                    for type_id in sorted(type_labels["improper"])
                ]

        for key, values in sections.items():
            if not values:
                continue
            array = np.array([line.split() for line in values])

            if key == "Atoms":
                cls._parse_atoms(array, topology, type_labels["atom"])
            elif key == "Velocities":
                cls._parse_velocities(array, topology)
            elif key == "Bonds":
                cls._parse_connectivity(
                    array, topology, "bonds", 2, type_labels["bond"]
                )
            elif key == "Angles":
                cls._parse_connectivity(
                    array, topology, "angles", 3, type_labels["angle"]
                )
            elif key == "Dihedrals":
                cls._parse_connectivity(
                    array, topology, "dihedrals", 4, type_labels["dihedral"]
                )
            elif key == "Impropers":
                cls._parse_connectivity(
                    array, topology, "impropers", 4, type_labels["improper"]
                )
            elif key == "Masses":
                cls._parse_masses(array, topology, type_labels["atom"])
            elif key in ("Pair Coeffs", "PairIJ Coeffs"):
                cls._parse_coeffs(
                    array,
                    topology,
                    "pair_params",
                    target_types=topology.get("atom_types"),
                )
            elif key == "Bond Coeffs":
                cls._parse_coeffs(
                    array, topology, "bond_params", target_types=topology["bond_types"]
                )
            elif key == "Angle Coeffs":
                cls._parse_coeffs(
                    array,
                    topology,
                    "angle_params",
                    target_types=topology["angle_types"],
                )
            elif key == "Dihedral Coeffs":
                cls._parse_coeffs(
                    array,
                    topology,
                    "dihedral_params",
                    target_types=topology["dihedral_types"],
                )
            elif key == "Improper Coeffs":
                cls._parse_coeffs(
                    array,
                    topology,
                    "improper_params",
                    target_types=topology["improper_types"],
                )
            elif key == "BondBond Coeffs":
                cls._parse_coeffs(
                    array,
                    topology,
                    "bondbond_params",
                    target_types=topology["angle_types"],
                )
            elif key == "BondAngle Coeffs":
                cls._parse_coeffs(
                    array,
                    topology,
                    "bondangle_params",
                    target_types=topology["angle_types"],
                )
            elif key == "AngleAngleTorsion Coeffs":
                cls._parse_coeffs(
                    array,
                    topology,
                    "angleangletorsion_params",
                    target_types=topology["dihedral_types"],
                )
            elif key == "AngleAngle Coeffs":
                cls._parse_coeffs(
                    array,
                    topology,
                    "angleangle_params",
                    target_types=topology["improper_types"],
                )

    @staticmethod
    def _parse_atoms(
        array: np.ndarray, topology: dict, type_labels: dict[int, str] | None = None
    ) -> None:
        """
        Parse the LAMMPS ``Atoms`` section.

        Supported styles
        ----------------
        atomic
            ``id type x y z``

        charge
            ``id type q x y z``

        full
            ``id mol type q x y z``

        Image flags ``ix iy iz`` are automatically detected and removed
        when the last three values of a line are integers.
        """

        # ------------------------------------------------------------------
        # Detect image flags
        # ------------------------------------------------------------------

        # Detect image flags (ix, iy, iz): present if the last 3 columns
        # are integer-valued.
        n_cols = array.shape[1]

        if n_cols >= 8:
            try:
                last_three = np.asarray(
                    array[:, -3:],
                    dtype=float,
                )

                is_integer = np.all(
                    np.isfinite(last_three) & (last_three == np.round(last_three))
                )

            except (ValueError, TypeError):
                is_integer = False

            if is_integer:
                array = array[:, :-3]
                n_cols -= 3

        # Infer the atom style from the remaining column count (see the
        # supported formats in the docstring above).
        if n_cols == 5:
            atom_style = "atomic"

            array = array[:, [0, 1, 2, 3, 4]]

            columns = [
                "id",
                "type",
                "x",
                "y",
                "z",
            ]

            has_charge = False

        elif n_cols == 6:
            atom_style = "charge"

            array = array[:, [0, 1, 2, 3, 4, 5]]

            columns = [
                "id",
                "type",
                "charge",
                "x",
                "y",
                "z",
            ]

            has_charge = True

        elif n_cols == 7:
            atom_style = "full"

            # The "mol" column (index 1) is ignored.
            array = array[:, [0, 2, 3, 4, 5, 6]]

            columns = [
                "id",
                "type",
                "charge",
                "x",
                "y",
                "z",
            ]

            has_charge = True

        else:
            raise ValueError(
                f"Unsupported LAMMPS Atoms format: {n_cols} columns "
                "after removing image flags. "
                "Expected 5, 6, or 7 columns."
            )

        df = pd.DataFrame(
            array,
            columns=columns,
        )

        df["id"] = pd.to_numeric(
            df["id"],
            errors="raise",
        ).astype(int)

        # Atom types can be numerical or strings.
        try:
            df["type"] = pd.to_numeric(
                df["type"],
                errors="raise",
            ).astype(int)
        except (ValueError, TypeError):
            df["type"] = df["type"].astype(str)

        if type_labels:
            df["type"] = df["type"].map(
                lambda atom_type: type_labels.get(atom_type, atom_type)
            )

        df[["x", "y", "z"]] = df[["x", "y", "z"]].astype(float)

        if has_charge:
            df["charge"] = df["charge"].astype(float)

            topology["charges"] = {
                atom_type: float(charge)
                for atom_type, charge in zip(
                    df["type"],
                    df["charge"],
                )
            }
        else:
            df["charge"] = 0.0
            topology["charges"] = {}

        df.set_index("id", inplace=True)

        topology["atoms"] = df
        topology["atom_style"] = atom_style

    @staticmethod
    def _parse_velocities(array: np.ndarray, topology: dict) -> None:
        """
        Parse Velocities section.

        Only supports the standard format: id vx vy vz.

        Parameters
        ----------
        array : np.ndarray
            Array containing the velocities data.
        topology : dict
            Topology dictionary to update with velocities.

        Notes
        -----
        The Velocities section in a LAMMPS data file has the format:
        atom-ID vx vy vz
        """
        # Take only the first 4 columns (id, vx, vy, vz).
        columns = ["id", "vx", "vy", "vz"]
        array = array[:, :4] if array.shape[1] >= 4 else array

        df = pd.DataFrame(array, columns=columns[: array.shape[1]])

        df["id"] = pd.to_numeric(df["id"], errors="coerce")
        for col in ["vx", "vy", "vz"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        df.set_index("id", inplace=True)
        topology["velocities"] = df

    @staticmethod
    def _parse_connectivity(
        array: np.ndarray,
        topology: dict,
        name: str,
        n_atoms: int,
        type_labels: dict[int, str] | None = None,
    ) -> None:
        """Parse Bonds, Angles, etc."""
        if array.size == 0:
            topology[name] = None
            return

        columns = ["id", "type"] + [f"atom_{i}" for i in range(1, n_atoms + 1)]
        df = pd.DataFrame(array, columns=columns)
        df["id"] = pd.to_numeric(df["id"], errors="coerce")

        # The type column is usually a numeric type id resolved through
        # type_labels, but with LAMMPS "Type Labels" the label can also be
        # spelled out directly in the section itself (this package's own
        # writer does exactly that) -- mirrors the same fallback used for
        # the Atoms section in `_parse_atoms`.
        try:
            df["type"] = pd.to_numeric(df["type"], errors="raise").astype(int)
            if type_labels:
                df["type"] = df["type"].map(
                    lambda connection_type: type_labels.get(
                        connection_type, connection_type
                    )
                )
        except (TypeError, ValueError):
            df["type"] = df["type"].astype(str)

        for col in columns[2:]:
            df[col] = df[col].astype(int)
        df.set_index("id", inplace=True)
        topology[name] = df

    @staticmethod
    def _parse_masses(
        array: np.ndarray,
        topology: dict,
        type_labels: dict[int, str] | None = None,
    ) -> None:
        """Parse Masses section.

        The type column is usually a numeric type id resolved through
        ``type_labels``, but files written with LAMMPS "Type Labels"
        (e.g. round-tripped through this package's own writer) may
        already spell the type out directly (e.g. "C" instead of "1").
        """

        def resolve_type(raw: Any) -> str:
            try:
                return (type_labels or {}).get(int(raw), str(raw))
            except (TypeError, ValueError):
                return str(raw)

        topology["masses"] = {
            resolve_type(row[0]): float(row[1]) for row in array
        }

    @staticmethod
    def _parse_coeffs(
        array: np.ndarray,
        topology: dict[str, Any],
        key: str,
        target_types: list[str] = None,
    ) -> None:
        """Converts text coefficients to appropriate Python dataclasses."""
        params = {}

        for row in array:
            style = None

            if key == "dihedral_params":
                row_text = " ".join(map(str, row))

                if "# fourier" in row_text:
                    style = "fourier"
                elif "# harmonic" in row_text:
                    style = "harmonic"

            type_id = int(row[0])
            label = (
                target_types[type_id - 1]
                if target_types and type_id <= len(target_types)
                else type_id
            )
            coeffs = [float(x) for x in row[1:]]

            # Reconstruction according to the key and the number of coefficients
            if key == "pair_params":
                if len(coeffs) == 2:
                    params[label] = LJParams(epsilon=coeffs[0], sigma=coeffs[1])
                elif len(coeffs) == 3:
                    params[label] = BuckinghamParams(
                        a=coeffs[0], rho=coeffs[1], c=coeffs[2]
                    )
                else:
                    params[label] = coeffs

            elif key == "bond_params":
                if len(coeffs) == 2:
                    params[label] = HarmonicBondParams(k=coeffs[0], r0=coeffs[1])
                elif len(coeffs) == 3:
                    params[label] = MorseBondParams(
                        D=coeffs[0], alpha=coeffs[1], r0=coeffs[2]
                    )
                elif len(coeffs) == 4:
                    params[label] = Class2BondParams(
                        r0=coeffs[0], k2=coeffs[1], k3=coeffs[2], k4=coeffs[3]
                    )
                else:
                    params[label] = coeffs

            elif key == "angle_params":
                if len(coeffs) == 2:
                    params[label] = HarmonicAngleParams(k=coeffs[0], theta0=coeffs[1])
                elif len(coeffs) == 4:
                    params[label] = Class2AngleParams(
                        theta0=coeffs[0], k2=coeffs[1], k3=coeffs[2], k4=coeffs[3]
                    )
                else:
                    params[label] = coeffs

            elif key == "bondbond_params" and len(coeffs) == 3:
                params[label] = Class2BondBondParams(
                    m=coeffs[0], r1=coeffs[1], r2=coeffs[2]
                )

            elif key == "bondangle_params" and len(coeffs) == 4:
                params[label] = Class2BondAngleParams(
                    n1=coeffs[0], n2=coeffs[1], r1=coeffs[2], r2=coeffs[3]
                )

            elif key == "dihedral_params":
                if style == "fourier":
                    m_terms = int(coeffs[0])
                    terms_list = []

                    idx = 1
                    for _ in range(m_terms):
                        k_val = float(coeffs[idx])
                        n_val = int(coeffs[idx + 1])
                        delta_val = float(coeffs[idx + 2])

                        terms_list.append(
                            FourierTerm(k=k_val, n=n_val, delta=delta_val)
                        )
                        idx += 3

                    params[label] = FourierDihedralParams(terms=terms_list)

                elif style == "harmonic":
                    params[label] = HarmonicDihedralParams(
                        k=coeffs[0],
                        d=int(coeffs[1]),
                        n=int(coeffs[2]),
                    )
                elif len(coeffs) == 4:
                    params[label] = CHARMMDihedralParams(
                        k=coeffs[0], n=int(coeffs[1]), d=int(coeffs[2]), w=coeffs[3]
                    )
                elif len(coeffs) == 6:
                    params[label] = Class2DihedralParams(
                        k1=coeffs[0],
                        phi1=coeffs[1],
                        k2=coeffs[2],
                        phi2=coeffs[3],
                        k3=coeffs[4],
                        phi3=coeffs[5],
                    )
                else:
                    params[label] = coeffs

            elif key == "angleangletorsion_params" and len(coeffs) == 3:
                params[label] = Class2AngleAngleTorsionParams(
                    m=coeffs[0], theta1=coeffs[1], theta2=coeffs[2]
                )

            elif key == "improper_params":
                if len(coeffs) == 2:
                    params[label] = HarmonicImproperParams(k=coeffs[0], chi0=coeffs[1])
                elif len(coeffs) == 3:
                    params[label] = CVFFImproperParams(
                        k=coeffs[0], d=int(coeffs[1]), n=int(coeffs[2])
                    )
                else:
                    params[label] = coeffs

            elif key == "angleangle_params" and len(coeffs) == 6:
                params[label] = Class2AngleAngleParams(
                    m1=coeffs[0],
                    m2=coeffs[1],
                    m3=coeffs[2],
                    theta1=coeffs[3],
                    theta2=coeffs[4],
                    theta3=coeffs[5],
                )

            else:
                params[label] = coeffs

        topology[key] = params


class LAMMPSWriter(BaseWriter):
    """Write LAMMPS data files."""

    @classmethod
    def write(
        cls, system, path: str, atom_style: str = "full", oldstyle: bool = False
    ) -> None:
        """Write system to LAMMPS data file.

        Parameters
        ----------
        system : AtomicSystem
            System to write.
        path : str
            Output file path.
        atom_style : str, optional
            LAMMPS atom style ('full' or 'charge').
        oldstyle : bool, optional
            If True, write in old style compatible with VMD/topotools.
            This uses numeric IDs instead of text labels and comments.
        """
        with open(path, "w", encoding="utf-8") as f:
            cls._write_header(f, system)
            cls._write_box(f, system)
            cls._write_labels(f, system, oldstyle)
            cls._write_coeffs(f, system, oldstyle)
            cls._write_masses(f, system, oldstyle)
            cls._write_atoms(f, system, atom_style)
            cls._write_velocities(f, system)
            cls._write_connectivity(f, system, "bonds")
            cls._write_connectivity(f, system, "angles")
            cls._write_connectivity(f, system, "dihedrals")
            cls._write_connectivity(f, system, "impropers")

    @staticmethod
    def _format_coeff_values(params: Any) -> str:
        """Converts a parameter object (Dataclass, list, tuple) into a character string for LAMMPS."""
        if isinstance(params, LJParams):
            return f"{params.epsilon:>12.6e} {params.sigma:>10.5f}"
        elif isinstance(params, BuckinghamParams):
            return f"{params.a:>12.6e} {params.rho:>10.5f} {params.c:>12.6e}"
        elif isinstance(params, HarmonicBondParams):
            return f"{params.k:>12.4f} {params.r0:>10.4f}"
        elif isinstance(params, MorseBondParams):
            return f"{params.D:>12.4f} {params.alpha:>10.4f} {params.r0:>10.4f}"
        elif isinstance(params, Class2BondParams):
            return f"{params.r0:>10.4f} {params.k2:>12.4f} {params.k3:>12.4f} {params.k4:>12.4f}"
        elif isinstance(params, FourierDihedralParams):
            terms_str = " ".join(
                f"{term.k:>12.4f} {term.n:>5d} {term.delta:>10.4f}"
                for term in params.terms
            )
            return f"{params.m:>3d} {terms_str}"
        elif isinstance(params, HarmonicDihedralParams):
            return f"{params.k:>12.4f} {params.d:>3} {params.n:>3}"
        elif isinstance(params, HarmonicAngleParams):
            return f"{params.k:>12.4f} {params.theta0:>10.4f}"
        elif isinstance(params, Class2AngleParams):
            return f"{params.theta0:>10.4f} {params.k2:>12.4f} {params.k3:>12.4f} {params.k4:>12.4f}"
        elif isinstance(params, Class2BondBondParams):
            return f"{params.m:>12.4f} {params.r1:>10.4f} {params.r2:>10.4f}"
        elif isinstance(params, Class2BondAngleParams):
            return f"{params.n1:>12.4f} {params.n2:>12.4f} {params.r1:>10.4f} {params.r2:>10.4f}"
        elif isinstance(params, HarmonicDihedralParams):
            return f"{params.k:>12.4f} {params.d:>3} {params.n:>3}"
        elif isinstance(params, CHARMMDihedralParams):
            return (
                f"{params.k:>12.4f} {params.n:>3} {int(params.d):>3} {params.w:>8.4f}"
            )
        elif isinstance(params, Class2DihedralParams):
            return f"{params.k1:>12.4f} {params.phi1:>10.4f} {params.k2:>12.4f} {params.phi2:>10.4f} {params.k3:>12.4f} {params.phi3:>10.4f}"
        elif isinstance(params, Class2AngleAngleTorsionParams):
            return f"{params.m:>12.4f} {params.theta1:>10.4f} {params.theta2:>10.4f}"
        elif isinstance(params, (HarmonicImproperParams, Class2ImproperParams)):
            return f"{params.k:>12.4f} {params.chi0:>10.4f}"
        elif isinstance(params, CVFFImproperParams):
            return f"{params.k:>12.4f} {params.d:>3} {params.n:>3}"
        elif isinstance(params, Class2AngleAngleParams):
            return f"{params.m1:>12.4f} {params.m2:>12.4f} {params.m3:>12.4f} {params.theta1:>10.4f} {params.theta2:>10.4f} {params.theta3:>10.4f}"
        elif isinstance(params, DistanceImproperParams):
            return f"{params.k2:>12.4f} {params.k4:>12.4f}"
        elif isinstance(params, (list, tuple)):
            return " ".join(
                f"{v:>12.4f}" if isinstance(v, (int, float)) else str(v) for v in params
            )
        else:
            return str(params)

    @classmethod
    def _write_section_coeffs(
        cls,
        f: TextIO,
        section_title: str,
        types_list: list[str],
        params_dict: dict[str, Any],
        oldstyle: bool,
    ) -> None:
        """Writes a standard section of coefficients (Bond Coeffs, Angle Coeffs, Class2, etc.)."""
        if not params_dict:
            return

        if oldstyle:
            f.write(f"\n# {section_title}\n#\n")
            for i, t in enumerate(types_list, 1):
                f.write(f"# {i} {t}\n")
            f.write("\n")
        else:
            f.write(f"\n {section_title}\n\n")
            for i, label in enumerate(types_list, 1):
                params = params_dict.get(label)
                if params is not None:
                    formatted_val = cls._format_coeff_values(params)
                    f.write(f"{i:<3} {formatted_val}   # {label}\n")
            f.write("\n")

    @classmethod
    def _write_coeffs(cls, f: TextIO, system: Any, oldstyle: bool) -> None:
        """Writes all the potential coefficients (pairs and associated connections)."""

        ff_params = getattr(system, "_ff_params", None)
        if ff_params is None:
            return

        if ff_params.pair:
            n = len(system.atom_types)
            actual_entries = len(ff_params.pair)

            all_self = all(k[0] == k[1] for k in ff_params.pair.keys())
            is_pair_coeffs = all_self and actual_entries == n

            section_name = "Pair Coeffs" if is_pair_coeffs else "PairIJ Coeffs"
            f.write(f"\n {section_name}\n\n")

            for i in range(1, n + 1):
                label_i = system.atom_types[i - 1]
                end_j = i if is_pair_coeffs else n

                for j in range(i, end_j + 1):
                    label_j = system.atom_types[j - 1]

                    key1 = f"{label_i}-{label_j}"
                    key2 = f"{label_j}-{label_i}"
                    params = ff_params.pair.get(key1) or ff_params.pair.get(key2)

                    if params is not None:
                        val_str = cls._format_coeff_values(params)
                        prefix = f"{i:<3}" if is_pair_coeffs else f"{i:<3} {j:<3}"
                        f.write(f"{prefix} {val_str}   # {label_i}-{label_j}\n")
            f.write("\n")

        sections_to_write = [
            (
                "Bond Coeffs",
                getattr(system, "bond_types", []),
                ff_params.bond,
            ),
            (
                "Angle Coeffs",
                getattr(system, "angle_types", []),
                ff_params.angle,
            ),
            (
                "Dihedral Coeffs",
                getattr(system, "dihedral_types", []),
                ff_params.dihedral,
            ),
            (
                "Improper Coeffs",
                getattr(system, "improper_types", []),
                ff_params.improper,
            ),
            # Class 2 sections
            (
                "BondBond Coeffs",
                getattr(system, "angle_types", []),
                ff_params.bondbond,
            ),
            (
                "BondAngle Coeffs",
                getattr(system, "angle_types", []),
                ff_params.bondangle,
            ),
            (
                "MiddleBondTorsion Coeffs",
                getattr(system, "dihedral_types", []),
                ff_params.middlebondtorsion,
            ),
            (
                "EndBondTorsion Coeffs",
                getattr(system, "dihedral_types", []),
                ff_params.endbondtorsion,
            ),
            (
                "AngleTorsion Coeffs",
                getattr(system, "dihedral_types", []),
                ff_params.angletorsion,
            ),
            (
                "AngleAngleTorsion Coeffs",
                getattr(system, "dihedral_types", []),
                ff_params.angleangletorsion,
            ),
            (
                "AngleAngle Coeffs",
                getattr(system, "improper_types", []),
                ff_params.angleangle,
            ),
        ]

        for title, types_list, params_dict in sections_to_write:
            cls._write_section_coeffs(f, title, types_list, params_dict, oldstyle)

    @staticmethod
    def _write_labels(f, system, oldstyle: bool) -> None:
        """Write type labels (only if not oldstyle)."""
        if oldstyle or not isinstance(system.atom_types[0], str):
            return

        f.write("\n Atom Type Labels\n\n")
        for i, t in enumerate(system.atom_types):
            f.write(f"{i + 1} {t}\n")

        for name in ["bond", "angle", "dihedral", "improper"]:
            types = getattr(system, f"{name}_types")
            if types:
                f.write(f"\n {name.capitalize()} Type Labels\n\n")
                for i, t in enumerate(types):
                    f.write(f"{i + 1} {t}\n")

    @staticmethod
    def _write_masses(f, system, oldstyle: bool) -> None:
        """Write atomic masses to the file stream."""
        f.write("\n Masses\n\n")

        # Ensure atom_types exists and is not empty
        if not hasattr(system, "atom_types") or not system.atom_types:
            return

        # Determine formatting mode (numeric IDs vs string labels)
        use_numeric_ids = oldstyle or not isinstance(system.atom_types[0], str)

        for i, atype in enumerate(system.atom_types, 1):
            mass = system.masses.get(atype, 0.0)

            if use_numeric_ids:
                f.write(f"{i} {mass}\n")
            else:
                f.write(f"{atype} {mass}\n")

    @staticmethod
    def _write_header(f, system) -> None:
        """Write file header."""
        now = datetime.datetime.now()
        f.write(
            f"LAMMPS data file generated by cemd on {now.strftime('%d/%m/%Y %H:%M:%S')}\n"
        )
        f.write(f" {system.num_atoms} atoms\n")
        if system.num_bonds:
            f.write(f" {system.num_bonds} bonds\n")
        if system.num_angles:
            f.write(f" {system.num_angles} angles\n")
        if system.num_dihedrals:
            f.write(f" {system.num_dihedrals} dihedrals\n")
        if system.num_impropers:
            f.write(f" {system.num_impropers} impropers\n")

        f.write(f" {system.num_atom_types} atom types\n")
        if system.num_bond_types:
            f.write(f" {system.num_bond_types} bond types\n")
        if system.num_angle_types:
            f.write(f" {system.num_angle_types} angle types\n")
        if system.num_dihedral_types:
            f.write(f" {system.num_dihedral_types} dihedral types\n")
        if system.num_improper_types:
            f.write(f" {system.num_improper_types} improper types\n")

    @staticmethod
    def _write_box(f, system) -> None:
        """Write box parameters."""
        xlo, xhi = system._box_lmp[0]
        ylo, yhi = system._box_lmp[1]
        zlo, zhi = system._box_lmp[2]
        f.write(f"   {xlo:>14.6f} {xhi:>14.6f} xlo xhi\n")
        f.write(f"   {ylo:>14.6f} {yhi:>14.6f} ylo yhi\n")
        f.write(f"   {zlo:>14.6f} {zhi:>14.6f} zlo zhi\n")

        xy, xz, yz = system._box_lmp[3]
        if any(abs(v) > 1e-8 for v in (xy, xz, yz)):
            f.write(f"   {xy:>14.6f} {xz:>14.6f} {yz:>14.6f} xy xz yz\n")

    @staticmethod
    def _write_atoms(f, system, atom_style: str) -> None:
        """Write atoms."""
        f.write(f"\n Atoms # {atom_style}\n\n")
        df = system.atoms.copy()

        if atom_style == "full":
            df.insert(0, "molecule", np.ones(len(df), dtype=int))
        # to_string rounds floats to 6 digits by default: with thousands of
        # atoms sharing the same fractional charge (e.g. after
        # neutralize_charge), the rounding error adds up to ~1e-3 e
        formatters = {"charge": lambda q: f"{q:.10f}"} if "charge" in df else None
        f.write(df.to_string(header=False, index_names=False, formatters=formatters))
        f.write("\n")

    @staticmethod
    def _write_velocities(f, system) -> None:
        """Write velocities."""
        if system.velocities is None:
            return
        velocities_sorted = system.velocities.sort_index()
        f.write("\n Velocities\n\n")
        f.write(velocities_sorted.to_string(header=False, index_names=False))
        f.write("\n")

    @staticmethod
    def _write_connectivity(f, system, name: str) -> None:
        """Write bonds, angles, etc."""
        df = getattr(system, name)
        if df is None or df.empty:
            return
        df_sorted = df.sort_index().copy()

        f.write(f"\n {name.capitalize()}\n\n")
        f.write(df_sorted.to_string(header=False, index_names=False))
        f.write("\n")
