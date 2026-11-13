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

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any

from ..forcefield_database import ForceFieldDatabase
from ..models import (
    AtomType,
    BuckinghamParams,
    Class2AngleAngleParams,
    Class2AngleAngleTorsionParams,
    Class2AngleParams,
    Class2BondAngleParams,
    Class2BondBondParams,
    Class2BondParams,
    DistanceImproperParams,
    ForceFieldModel,
    HarmonicAngleParams,
    HarmonicBondParams,
    HarmonicImproperParams,
    LJParams,
)


@dataclass
class ParseResult:
    """Result of force field analysis."""

    model_name: str
    atoms: dict[str, AtomType] = field(default_factory=dict)
    lj: dict[str, LJParams] = field(default_factory=dict)
    buckingham: dict[str, BuckinghamParams] = field(default_factory=dict)
    bonds: dict[str, HarmonicBondParams | Class2BondParams] = field(
        default_factory=dict
    )
    angles: dict[str, HarmonicAngleParams | Class2AngleParams] = field(
        default_factory=dict
    )
    impropers: dict[str, HarmonicImproperParams | DistanceImproperParams] = field(
        default_factory=dict
    )
    dihedrals: dict[str, Any] = field(default_factory=dict)  # Specific format
    bondbond: dict[str, Class2BondBondParams] = field(default_factory=dict)
    bondangle: dict[str, Class2BondAngleParams] = field(default_factory=dict)
    angleangletorsion: dict[str, Class2AngleAngleTorsionParams] = field(
        default_factory=dict
    )
    angleangle: dict[str, Class2AngleAngleParams] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_database(self, db: "ForceFieldDatabase") -> None:
        """
        Loads this result into a ForceField Database.
        """
        if self.model_name not in db.models:
            db.models[self.model_name] = ForceFieldModel(
                name=self.metadata.get("name", self.model_name),
                description=self.metadata.get("description", ""),
                ref=self.metadata.get("ref", ""),
                tags=self.metadata.get("tags", []),
                category=self.metadata.get("category", ""),
                mixing=self.metadata.get("mixing", ""),
                mixing_note=self.metadata.get("mixing_note", ""),
            )

        for key, value in self.atoms.items():
            full_key = f"{self.model_name}.{key}"
            db.atom[full_key] = value

        for key, value in self.lj.items():
            full_key = f"{self.model_name}.{key}"
            db.lj[full_key] = value

        for key, value in self.buckingham.items():
            full_key = f"{self.model_name}.{key}"
            db.buckingham[full_key] = value

        for key, value in self.bonds.items():
            full_key = f"{self.model_name}.{key}"
            db.bond[full_key] = value

        for key, value in self.angles.items():
            full_key = f"{self.model_name}.{key}"
            db.angle[full_key] = value

        for key, value in self.impropers.items():
            full_key = f"{self.model_name}.{key}"
            db.improper[full_key] = value

        for key, value in self.bondbond.items():
            full_key = f"{self.model_name}.{key}"
            db.bondbond[full_key] = value

        for key, value in self.bondangle.items():
            full_key = f"{self.model_name}.{key}"
            db.bondangle[full_key] = value

        if hasattr(db, "angleangletorsion"):
            for key, value in self.angleangletorsion.items():
                full_key = f"{self.model_name}.{key}"
                db.angleangletorsion[full_key] = value

        if hasattr(db, "angleangle"):
            for key, value in self.angleangle.items():
                full_key = f"{self.model_name}.{key}"
                db.angleangle[full_key] = value

        # Specific dihedral
        if hasattr(db, "dihedral"):
            for key, value in self.dihedrals.items():
                full_key = f"{self.model_name}.{key}"
                db.dihedral[full_key] = value


class BaseForceFieldParser(ABC):
    "Base class for all force fields parsers."

    def __init__(self):
        self.model_name: str = ""
        self.result: ParseResult = ParseResult(model_name="")

    @abstractmethod
    def parse(self, content: str) -> ParseResult:
        """
        Parses the contents of a force field file.

        Settings
        ----------
        content: str
            Content of the file to parse

        Returns
        -------
        ParseResult
            Structured analysis result
        """
        pass

    @abstractmethod
    def parse_file(self, filepath: str) -> ParseResult:
        """
        Parse a force field file.

        Settings
        ----------
        filepath: str
            Path to file

        Returns
        -------
        ParseResult
            Structured analysis result
        """
        pass

    def _clean_content(self, content: str) -> list[str]:
        """Cleans the content (removes comments, empty lines)."""
        lines = []
        for line in content.split("\n"):
            # Remove comments according to format
            line = self._strip_comments(line)
            line = line.strip()
            if line:
                lines.append(line)
        return lines

    def _strip_comments(self, line: str) -> str:
        """
        Remove one-line comments.
        To be overwritten depending on the format.
        """
        if "#" in line:
            return line[: line.index("#")]
        return line
