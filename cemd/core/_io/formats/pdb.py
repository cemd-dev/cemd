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
import warnings

import numpy as np
import pandas as pd

from ...._constants import MASSES_DICT

# ----------------------------------------------------------------------
# Serial numbers past 99 999
#
# A PDB serial number has five columns, so plain decimal stops at 99 999.
# Writers then disagree. Packmol and VMD switch to hexadecimal ("186A0" is
# 100 000), which cannot be told from decimal field by field: "19000" is
# 19 000 or 0x19000. Hybrid-36 (http://cci.lbl.gov/hybrid_36/) goes on in
# base 36 from "A0000" and is never ambiguous, so `PDBWriter` uses it; the
# reader recognises both.
# ----------------------------------------------------------------------

_WIDTH = 5
_DECIMAL_MAX = 10**_WIDTH
_BLOCK = 26 * 36 ** (_WIDTH - 1)
_LETTER_OFFSET = 10 * 36 ** (_WIDTH - 1)
_DIGITS36 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def hy36decode(field: str) -> int:
    """Read a five-column serial number, decimal or hybrid-36."""
    text = field.strip()
    if not text:
        raise ValueError("empty serial number")
    if text[0].isupper():
        return int(text, 36) - _LETTER_OFFSET + _DECIMAL_MAX
    if text[0].islower():
        return int(text, 36) - _LETTER_OFFSET + _DECIMAL_MAX + _BLOCK
    return int(text)


def hy36encode(value: int) -> str:
    """Write `value` in five columns, decimal up to 99 999, then hybrid-36."""
    if 0 <= value < _DECIMAL_MAX:
        return f"{value:{_WIDTH}d}"
    rest = value - _DECIMAL_MAX
    if not 0 <= rest < 2 * _BLOCK:
        raise ValueError(f"{value} does not fit in a PDB serial number")
    upper = rest < _BLOCK
    rest = (rest if upper else rest - _BLOCK) + _LETTER_OFFSET
    digits = ""
    while rest:
        rest, digit = divmod(rest, 36)
        digits = _DIGITS36[digit] + digits
    return digits if upper else digits.lower()


class _SerialDecoder:
    """Decode the serial numbers of one file, in order.

    Decimal until the first field that is not; then hexadecimal or
    hybrid-36, whichever continues the count, for the rest of the file.
    """

    def __init__(self) -> None:
        self.previous = 0
        self.scheme = None

    def __call__(self, text: str) -> int:
        if self.scheme == "hex":
            value = int(text, 16)
        elif self.scheme is None and text.isdigit():
            value = int(text)
        else:
            value = (
                self._first_non_decimal(text)
                if self.scheme is None
                else (hy36decode(text))
            )
        self.previous = value
        return value

    def _first_non_decimal(self, text: str) -> int:
        candidates = {}
        for scheme, decode in (("hy36", hy36decode), ("hex", lambda t: int(t, 16))):
            try:
                candidates[scheme] = decode(text)
            except ValueError:
                pass
        if not candidates:
            raise ValueError(f"unreadable serial number {text!r}")
        follows = [s for s, v in candidates.items() if v == self.previous + 1]
        self.scheme = follows[0] if follows else next(iter(candidates))
        return candidates[self.scheme]


class PDBReader:
    """Column offset tolerant PDB file reader."""

    @classmethod
    def read(cls, path: str) -> dict:
        atoms_data = []
        bonds_data = []
        box_params = None
        decode_serial = _SerialDecoder()

        with open(path, encoding="utf-8") as f:
            lines = f.readlines()

        for line in lines:
            line_type = line[:6].strip()

            # Analysis of the box (CRYST1)
            if line_type == "CRYST1":
                parts = line.split()
                if len(parts) >= 7:
                    box_params = tuple(float(p) for p in parts[1:7])

            # Parsing atoms (ATOM / HETATM)
            elif line_type in ("ATOM", "HETATM"):
                serial, local_name, x, y, z = cls._parse_atom_line(line)
                atom_id = decode_serial(serial)

                # The standard PDB format has no charge field: columns
                # after x/y/z are occupancy and temperature factor (both
                # written as fixed placeholders by `PDBWriter`, and by
                # Packmol's own output). Reading `floats[3]` back as
                # "charge" was actually reading occupancy (always 1.00),
                # silently overwriting every atom's real charge with 1.0
                # on any PDB round-trip.
                charge = 0.0

                # The element column is written upper-cased ("CA", "SI") by
                # the PDB standard, so it has to be normalized before being
                # looked up in the mass table. An entry that still doesn't
                # match a known element (a truncated symbol from an older
                # file, a non-standard label) falls back to the atom name.
                element = line[76:78].strip() if len(line) >= 78 else ""
                element = element.capitalize()
                if element not in MASSES_DICT:
                    element = cls._guess_element(local_name)

                atoms_data.append(
                    {
                        "id": atom_id,
                        "type": local_name,
                        "charge": charge,
                        "mass": cls._guess_mass(element),
                        "x": x,
                        "y": y,
                        "z": z,
                    }
                )

            # Connectivity analysis (CONECT)
            elif line_type == "CONECT":
                source_id, *targets = cls._parse_conect_line(line)
                for target_id in targets:
                    if source_id < target_id:
                        bonds_data.append((source_id, target_id))

        # Creation of DataFrames
        df_atoms = pd.DataFrame(atoms_data)
        if not df_atoms.empty and not df_atoms["id"].is_unique:
            # Serial numbers that repeat could not be decoded: number the
            # atoms in file order, and drop CONECT records, which would
            # point at the wrong atoms.
            warnings.warn(
                f"{path}: atom serial numbers repeat; the atoms are numbered "
                "in file order and CONECT records are ignored.",
                UserWarning,
                stacklevel=2,
            )
            df_atoms["id"] = range(1, len(df_atoms) + 1)
            bonds_data = []
        if not df_atoms.empty:
            df_atoms.set_index("id", inplace=True)

        df_bonds = None
        if bonds_data:
            pairs = np.asarray(bonds_data, dtype=np.int64)
            types = df_atoms["type"].astype(str)
            t1 = pd.Series(pairs[:, 0]).map(types).fillna("X")
            t2 = pd.Series(pairs[:, 1]).map(types).fillna("X")
            df_bonds = pd.DataFrame(
                {"type": t1 + "-" + t2, "atom_1": pairs[:, 0], "atom_2": pairs[:, 1]},
                index=pd.RangeIndex(1, len(pairs) + 1, name="id"),
            )

        # Configuration of the simulation box
        if box_params:
            a, b, c, _, _, _ = box_params
            lmp_box = ((0.0, a), (0.0, b), (0.0, c), (0.0, 0.0, 0.0))
        elif not df_atoms.empty:
            coords = df_atoms[["x", "y", "z"]].values
            mins = coords.min(axis=0) - 10.0
            maxs = coords.max(axis=0) + 10.0
            lmp_box = (
                (mins[0], maxs[0]),
                (mins[1], maxs[1]),
                (mins[2], maxs[2]),
                (0.0, 0.0, 0.0),
            )
        else:
            lmp_box = ((0.0, 10.0), (0.0, 10.0), (0.0, 10.0), (0.0, 0.0, 0.0))

        charges_dict = {row["type"]: row["charge"] for row in atoms_data}
        masses_dict = {row["type"]: row["mass"] for row in atoms_data}

        return {
            "atoms": df_atoms,
            "bonds": df_bonds,
            "angles": None,
            "dihedrals": None,
            "impropers": None,
            "velocities": None,
            "box": lmp_box,
            "masses": masses_dict,
            "charges": charges_dict,
            "atom_style": "full",
        }

    @staticmethod
    def _parse_atom_line(line: str) -> tuple[str, str, float, float, float]:
        """Serial number, atom name and coordinates of an ATOM/HETATM line.

        The serial number is returned as written; `_SerialDecoder` reads
        it. Read from the standard columns first. Splitting on whitespace broke
        from 10 000 atoms on, where the five-digit serial fills its columns
        and meets the record name ("HETATM10001"): every Packmol build that
        large failed with "invalid literal for int()". Whitespace splitting
        is kept for files whose columns are shifted.
        """
        try:
            serial = line[6:11].strip()
            if not serial:
                raise ValueError("empty serial number")
            return (
                serial,
                line[12:16].strip(),
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
        except ValueError:
            pass

        parts = line.split()
        atom_id = parts[1]
        local_name = parts[2]

        # The 3 coordinates are the first 3 decimal numbers of the line.
        coords_match = re.findall(r"[-+]?\d*\.\d+|\d+", line)
        floats = [float(val) for val in coords_match if "." in val]
        if len(floats) >= 3:
            x, y, z = floats[0], floats[1], floats[2]
        else:
            x, y, z = float(parts[-5]), float(parts[-4]), float(parts[-3])
        return atom_id, local_name, x, y, z

    @staticmethod
    def _parse_conect_line(line: str) -> list[int]:
        """Source atom then bonded atoms of a CONECT line, five columns each."""
        fields = [line[i : i + 5] for i in range(6, 31, 5)]
        try:
            return [hy36decode(f) for f in fields if f.strip()]
        except ValueError:
            return [int(part) for part in line.split()[1:]]

    @staticmethod
    def _guess_element(atom_type: str) -> str:
        cleaned = re.sub(r"[^a-zA-Z]", "", atom_type)
        if not cleaned:
            return "C"
        if len(cleaned) >= 2 and cleaned[1].islower():
            return cleaned[:2].capitalize()
        return cleaned[0].upper()

    @classmethod
    def _guess_mass(cls, element: str) -> float:
        return MASSES_DICT.get(element, 12.011)


class PDBWriter:
    """Writing an AtomicSystem object in strict PDB format (Packmol compatible)."""

    @classmethod
    def write(cls, system, path: str) -> None:
        """
        Exports an AtomicSystem to a .pdb file conforming to the PDB standard.

        Settings
        ----------
        system: AtomicSystem
            Instance of the atomic system to export.
        path : str
            Path of the output file.
        """
        lines = []

        # Box header (CRYST1)
        if hasattr(system, "box") and system.box is not None:
            box = system.box
            a, b, c = box[0], box[1], box[2]
            alpha = box[3] if len(box) > 3 else 90.0
            beta = box[4] if len(box) > 4 else 90.0
            gamma = box[5] if len(box) > 5 else 90.0

            cryst_line = (
                f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}"
                f"{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1"
            )
            lines.append(cryst_line)

        # Atoms section (HETATM)
        atoms_df = system.atoms

        # Elements are inferred from the masses carried by the system, which
        # is the only reliable source here: the type names are force-field
        # labels ("Ow", "Hw", "Oc"...) and truncating them to their first
        # letter turned every two-letter element into another one -- "Ca"
        # was written as carbon, "Si" as sulfur -- so a PDB round-trip (i.e.
        # every Packmol build) silently replaced their masses.
        elements = dict(system.elements)

        # One pass per column instead of `iterrows`, which built a Series
        # per atom: writing 24 000 atoms took 1.7 s.
        types = atoms_df["type"].astype(str).to_numpy()
        names = [f" {t:<3s}" if len(t) < 4 else f"{t[:4]:<4s}" for t in types]

        # Element symbol (ex: C, Ca, Si), upper-cased as the PDB standard
        # requires, right-justified in columns 77-78.
        element_of = {}
        for raw in atoms_df["type"].unique():
            element = elements.get(raw, "") or PDBReader._guess_element(str(raw))
            element_of[str(raw)] = element.upper()[:2]

        # Strict Fortran/PDB column layout required by Packmol:
        # COLS 1-6 : "HETATM"
        # COLS 7-11 : Atom serial number (decimal, then hybrid-36)
        # COLS 13-16 : Atom name (%4s)
        # COLS 18-20 : Residue name (%3s)
        # COL 22 : Chain identifier (%1s)
        # COLS 23-26 : Residue sequence number (%4d)
        # COLS 31-38, 39-46, 47-54 : X, Y, Z (%8.3f)
        # COLS 55-60 : Occupancy (%6.2f)
        # COLS 61-66 : Temp factor (%6.2f)
        # COLS 77-78 : Element symbol (%2s)
        xyz = atoms_df[["x", "y", "z"]].to_numpy(dtype=float)
        for atom_id, name, atom_type, (x, y, z) in zip(
            atoms_df.index, names, types, xyz
        ):
            lines.append(
                f"HETATM{hy36encode(int(atom_id))} {name} "
                f"UNK A   1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"  1.00  0.00          "
                f"{element_of[atom_type]:>2s}"
            )

        # Connections section (CONECT) - optional
        if (
            hasattr(system, "bonds")
            and system.bonds is not None
            and not system.bonds.empty
        ):
            adj = {}
            pairs = system.bonds[["atom_1", "atom_2"]].to_numpy(dtype=np.int64)
            for a1, a2 in pairs.tolist():
                adj.setdefault(a1, []).append(a2)
                adj.setdefault(a2, []).append(a1)

            for source_id in sorted(adj.keys()):
                targets = sorted(adj[source_id])
                for i in range(0, len(targets), 4):
                    chunk = targets[i : i + 4]
                    targets_str = "".join(hy36encode(t) for t in chunk)
                    lines.append(f"CONECT{hy36encode(source_id)}{targets_str}")

        lines.append("END")

        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
