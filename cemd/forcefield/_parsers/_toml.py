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

import tomllib
from typing import Any

from ..models import (
    AtomType,
    BuckinghamParams,
    Class2AngleParams,
    Class2BondAngleParams,
    Class2BondBondParams,
    Class2BondParams,
    DistanceImproperParams,
    HarmonicAngleParams,
    HarmonicBondParams,
    HarmonicImproperParams,
    LJParams,
)
from ._base import BaseForceFieldParser, ParseResult


class TOMLParser(BaseForceFieldParser):
    """
    Parser pour les fichiers TOML au format CEMD.
    """

    def __init__(self):
        super().__init__()

    def parse(self, content: str, model_name: str | None = None) -> ParseResult:
        """Parse une chaîne TOML."""
        data = tomllib.loads(content)
        return self._parse_data(data, model_name)

    def parse_file(self, filepath: str, model_name: str | None = None) -> ParseResult:
        """Parse un fichier TOML."""
        with open(filepath, "rb") as f:
            data = tomllib.load(f)
        return self._parse_data(data, model_name)

    def _parse_data(
        self, data: dict[str, Any], model_name: str | None = None
    ) -> ParseResult:
        """Parse les données TOML."""
        # The database key: the caller's file-derived name (e.g. "clayff")
        # takes precedence, so ff_keys stay lowercase and filesystem-stable
        # regardless of how [model.name] is capitalized in the TOML.
        display_name = data.get("model", {}).get("name", "unknown")
        model_key = model_name if model_name is not None else display_name

        result = ParseResult(model_name=model_key)

        # Metadata
        if "model" in data:
            result.metadata = {
                "name": display_name,
                "description": data["model"].get("description", ""),
                "ref": data["model"].get("ref", ""),
                "tags": data["model"].get("tags", []),
            }

        # Atoms
        for key, params in data.get("atom", {}).items():
            result.atoms[key] = AtomType(
                element=params["element"],
                # Absent means "this force field has no per-type charge",
                # which must not be applied as a 0.0 over a real charge.
                charge=params.get("charge"),
                environment=params.get("environment", ""),
                ref=params.get("ref", ""),
                mass=params.get("mass"),
                model=model_key,
            )

        # LJ
        for key, params in data.get("lj", {}).items():
            result.lj[key] = LJParams(
                epsilon=params["epsilon"],
                sigma=params["sigma"],
                ref=params.get("ref", ""),
                model=model_key,
            )

        # Buckingham
        for key, params in data.get("buckingham", {}).items():
            result.buckingham[key] = BuckinghamParams(
                a=params["A"],
                rho=params["rho"],
                c=params.get("C", 0.0),
                ref=params.get("ref", ""),
                model=model_key,
            )

        # Harmonic bonds
        for key, params in data.get("bond", {}).get("harmonic", {}).items():
            result.bonds[key] = HarmonicBondParams(
                k=params["k"],
                r0=params["r0"],
                ref=params.get("ref", ""),
                model=model_key,
            )

        # Class2 bonds
        for key, params in data.get("bond", {}).get("class2", {}).items():
            result.bonds[key] = Class2BondParams(
                r0=params["r0"],
                k2=params["k2"],
                k3=params.get("k3", 0.0),
                k4=params.get("k4", 0.0),
                ref=params.get("ref", ""),
                model=model_key,
            )

        # Harmonic angles
        for key, params in data.get("angle", {}).get("harmonic", {}).items():
            result.angles[key] = HarmonicAngleParams(
                k=params["k"],
                theta0=params["theta0"],
                ref=params.get("ref", ""),
                model=model_key,
            )

        # Class2 angles
        for key, params in data.get("angle", {}).get("class2", {}).items():
            result.angles[key] = Class2AngleParams(
                theta0=params["theta0"],
                k2=params["k2"],
                k3=params.get("k3", 0.0),
                k4=params.get("k4", 0.0),
                ref=params.get("ref", ""),
                model=model_key,
            )

        # Harmonic impropers
        for key, params in data.get("improper", {}).get("harmonic", {}).items():
            result.impropers[key] = HarmonicImproperParams(
                k=params["k"],
                chi0=params.get("chi0", 0.0),
                ref=params.get("ref", ""),
                model=model_key,
            )

        # Distance impropers
        for key, params in data.get("improper", {}).get("distance", {}).items():
            result.impropers[key] = DistanceImproperParams(
                k2=params["k2"],
                k4=params.get("k4", 0.0),
                ref=params.get("ref", ""),
                model=model_key,
            )

        # Bondbond
        for key, params in data.get("bondbond", {}).get("class2", {}).items():
            result.bondbond[key] = Class2BondBondParams(
                m=params["m"],
                r1=params["r1"],
                r2=params["r2"],
                ref=params.get("ref", ""),
                model=model_key,
            )

        # Bondangle
        for key, params in data.get("bondangle", {}).get("class2", {}).items():
            result.bondangle[key] = Class2BondAngleParams(
                n1=params["n1"],
                n2=params["n2"],
                r1=params["r1"],
                r2=params["r2"],
                ref=params.get("ref", ""),
                model=model_key,
            )

        # Dihedrals (specific format)
        if "dihedral" in data:
            result.dihedrals = data["dihedral"]

        return result
