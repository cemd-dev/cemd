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

from pathlib import Path
from typing import Any

import pandas as pd

from .models import (
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


class ForceFieldDatabase:
    """
    Force field database that loads all available force fields from the database directory.

    Usage:
        db = ForceFieldDatabase()           # Load all available force fields
        db = ForceFieldDatabase(db_dir="path/to/db")  # Custom directory

    Attributes
    ----------
    atom : dict[str, AtomType]
        Dictionary of atom types {full_name: AtomType}
    lj : dict[str, LJParams]
        Dictionary of LJ parameters {model.pair: LJParams}
    buckingham : dict[str, BuckinghamParams]
        Dictionary of Buckingham parameters {model.pair: BuckinghamParams}
    bond : dict[str, HarmonicBondParams | Class2BondParams]
        Dictionary of bond parameters {model.pair: BondParams}
    angle : dict[str, HarmonicAngleParams | Class2AngleParams]
        Dictionary of angle parameters {model.triple: AngleParams}
    improper : dict[str, HarmonicImproperParams | DistanceImproperParams]
        Dictionary of improper parameters {model.quad: ImproperParams}
    bondbond : dict[str, Class2BondBondParams]
        Dictionary of bond-bond parameters
    bondangle : dict[str, Class2BondAngleParams]
        Dictionary of bond-angle parameters
    angleangletorsion : dict[str, Class2AngleAngleTorsionParams]
        Dictionary of angle-angle-torsion (class2) cross-term parameters
    angleangle : dict[str, Class2AngleAngleParams]
        Dictionary of angle-angle (class2) improper cross-term parameters
    models : dict[str, ForceFieldModel]
        Dictionary of model metadata {name: ForceFieldModel}
    dihedral : dict[str, Any]
        Dictionary of dihedral parameters (format specific)
    """

    def __init__(self, db_dir: str | Path | None = None):
        """
        Initialize the force field database and load all available force fields.

        Parameters
        ----------
        db_dir : str or Path, optional
            Custom directory containing force field database files.
            Defaults to standard package database directory.
        """
        # Determine the forcefields folder path
        self.db_dir = Path(db_dir) if db_dir is not None else self._get_forcefield_dir()

        self.atom: dict[str, AtomType] = {}
        self.lj: dict[str, LJParams] = {}
        self.buckingham: dict[str, BuckinghamParams] = {}
        self.bond: dict[str, HarmonicBondParams | Class2BondParams] = {}
        self.angle: dict[str, HarmonicAngleParams | Class2AngleParams] = {}
        self.improper: dict[str, HarmonicImproperParams | DistanceImproperParams] = {}
        self.bondbond: dict[str, Class2BondBondParams] = {}
        self.bondangle: dict[str, Class2BondAngleParams] = {}
        self.angleangletorsion: dict[str, Class2AngleAngleTorsionParams] = {}
        self.angleangle: dict[str, Class2AngleAngleParams] = {}
        self.models: dict[str, ForceFieldModel] = {}
        self.dihedral: dict[str, Any] = {}

        # Load all available forcefields
        self._load_all()

    @staticmethod
    def _get_forcefield_dir() -> Path:
        """Get the forcefield database directory."""
        current_dir = Path(__file__).parent
        default_dir = current_dir / "db"
        if default_dir.exists():
            return default_dir

        cwd_dir = Path.cwd() / "db"
        if cwd_dir.exists():
            return cwd_dir

        default_dir.mkdir(parents=True, exist_ok=True)
        return default_dir

    def _load_all(self) -> None:
        """Load all supported force field files from the database directory."""
        for filepath in sorted(self.db_dir.iterdir()):
            if not filepath.is_file():
                continue

            filename = filepath.name

            if filename.endswith(".toml"):
                model_name = filepath.stem
                self._load_toml_model(filepath, model_name)

            elif filename.endswith(".lt"):
                self._load_gromos_lt(filepath, "gromos")

            elif filename.endswith(".prm"):
                self._load_charmm_prm(filepath, "iff_charmm")

            elif filename.endswith(".frc"):
                self._load_cvff_frc(filepath, "iff_cvff")

    def _load_toml_model(self, filepath: Path, model_name: str) -> None:
        """Load a TOML force field file using the TOMLParser.

        `model_name` (the file's stem, e.g. "clayff") is the ff_key prefix
        used for every entry, regardless of how [model.name] is capitalized
        in the TOML -- that display name only feeds ForceFieldModel.name.
        """
        from ._parsers._toml import TOMLParser

        parser = TOMLParser()
        parse_result = parser.parse_file(str(filepath), model_name=model_name)
        parse_result.to_database(self)

    def _load_gromos_lt(self, filepath: Path, model_name: str) -> None:
        """Load a GROMOS force field from a moltemplate .lt file."""
        from ._parsers._gromos import GromosForceFieldLoader

        loader = GromosForceFieldLoader(self)
        loader.load_from_file(filepath, model_name)

    def _load_charmm_prm(self, filepath: Path, model_name: str) -> None:
        """Load a CHARMM force field from a .prm file."""
        from ._parsers._iff_charmm import CHARMMInterfaceLoader

        loader = CHARMMInterfaceLoader(self)
        loader.load_from_file(filepath, model_name)

    def _load_cvff_frc(self, filepath: Path, model_name: str) -> None:
        """Load a CVFF force field from a .frc file."""
        from ._parsers._iff_cvff import CVFFInterfaceLoader

        loader = CVFFInterfaceLoader(self)
        loader.load_from_file(filepath, model_name)

    def _extract_short_name(self, name: str) -> tuple[str | None, str]:
        """Extract model and short name from a full name."""
        if "." in name:
            model, short = name.split(".", 1)
            return model, short
        return None, name

    def get_atom_type(self, name: str) -> AtomType | None:
        """Get atom type by full name (model.type) or short name."""
        if name in self.atom:
            return self.atom[name]

        _, short = self._extract_short_name(name)
        for full_key, params in self.atom.items():
            if full_key.endswith(f".{short}"):
                return params
        return None

    def _get_pair_params(
        self, dict_obj: dict, type1: str, type2: str, model: str = None
    ):
        """Generic method to search for pair parameters."""
        _, short1 = self._extract_short_name(type1)
        _, short2 = self._extract_short_name(type2)

        if model is None:
            m1, _ = self._extract_short_name(type1)
            if m1 is not None:
                model = m1

        key = f"{short1}-{short2}"
        key_rev = f"{short2}-{short1}"

        if model:
            for k in [f"{model}.{key}", f"{model}.{key_rev}"]:
                if k in dict_obj:
                    return dict_obj[k]
            return None

        for full_key, params in dict_obj.items():
            if "." in full_key:
                pair_key = full_key.split(".", 1)[1]
                if pair_key == key or pair_key == key_rev:
                    return params
        return None

    def _get_triple_params(
        self, dict_obj: dict, type1: str, type2: str, type3: str, model: str = None
    ):
        """Generic method to search for triple parameters."""
        _, short1 = self._extract_short_name(type1)
        _, short2 = self._extract_short_name(type2)
        _, short3 = self._extract_short_name(type3)

        if model is None:
            m1, _ = self._extract_short_name(type1)
            if m1 is not None:
                model = m1

        key = f"{short1}-{short2}-{short3}"
        key_rev = f"{short3}-{short2}-{short1}"

        if model:
            for k in [f"{model}.{key}", f"{model}.{key_rev}"]:
                if k in dict_obj:
                    return dict_obj[k]
            return None

        for full_key, params in dict_obj.items():
            if "." in full_key:
                triple_key = full_key.split(".", 1)[1]
                if triple_key == key or triple_key == key_rev:
                    return params
        return None

    def _get_quadruple_params(
        self,
        dict_obj: dict,
        type1: str,
        type2: str,
        type3: str,
        type4: str,
        model: str = None,
    ):
        """Generic method to search for quadruple parameters."""
        _, short1 = self._extract_short_name(type1)
        _, short2 = self._extract_short_name(type2)
        _, short3 = self._extract_short_name(type3)
        _, short4 = self._extract_short_name(type4)

        if model is None:
            m1, _ = self._extract_short_name(type1)
            if m1 is not None:
                model = m1

        key = f"{short1}-{short2}-{short3}-{short4}"
        key_rev = f"{short4}-{short3}-{short2}-{short1}"

        if model:
            for k in [f"{model}.{key}", f"{model}.{key_rev}"]:
                if k in dict_obj:
                    return dict_obj[k]
            return None

        for full_key, params in dict_obj.items():
            if "." in full_key:
                quad_key = full_key.split(".", 1)[1]
                if quad_key == key or quad_key == key_rev:
                    return params
        return None

    def get_lj(self, type1: str, type2: str, model: str = None) -> LJParams | None:
        """Get Lennard-Jones 12-6 parameters for a pair."""
        return self._get_pair_params(self.lj, type1, type2, model)

    def get_buckingham(
        self, type1: str, type2: str, model: str = None
    ) -> BuckinghamParams | None:
        """Get Buckingham potential parameters for a pair."""
        return self._get_pair_params(self.buckingham, type1, type2, model)

    def get_bond(
        self, type1: str, type2: str, model: str = None
    ) -> HarmonicBondParams | Class2BondParams | None:
        """Get bond parameters for a pair."""
        return self._get_pair_params(self.bond, type1, type2, model)

    def get_angle(
        self, type1: str, type2: str, type3: str, model: str = None
    ) -> HarmonicAngleParams | Class2AngleParams | None:
        """Get angle parameters for a triple of atom types."""
        return self._get_triple_params(self.angle, type1, type2, type3, model)

    def get_bondbond(
        self, type1: str, type2: str, type3: str, model: str = None
    ) -> Class2BondBondParams | None:
        """Get bondbond parameters for a triple of atom types."""
        return self._get_triple_params(self.bondbond, type1, type2, type3, model)

    def get_bondangle(
        self, type1: str, type2: str, type3: str, model: str = None
    ) -> Class2BondAngleParams | None:
        """Get bondangle parameters for a triple of atom types."""
        return self._get_triple_params(self.bondangle, type1, type2, type3, model)

    def get_angleangletorsion(
        self, type1: str, type2: str, type3: str, type4: str, model: str = None
    ) -> Class2AngleAngleTorsionParams | None:
        """Get angle-angle-torsion parameters for a quadruple of atom types."""
        return self._get_quadruple_params(
            self.angleangletorsion, type1, type2, type3, type4, model
        )

    def get_angleangle(
        self, type1: str, type2: str, type3: str, type4: str, model: str = None
    ) -> Class2AngleAngleParams | None:
        """Get angle-angle parameters for a quadruple of atom types."""
        return self._get_quadruple_params(
            self.angleangle, type1, type2, type3, type4, model
        )

    def get_dihedral(
        self, type1: str, type2: str, type3: str, type4: str, model: str = None
    ) -> Any | None:
        """Get dihedral parameters for a quadruple of atom types."""
        return self._get_quadruple_params(
            self.dihedral, type1, type2, type3, type4, model
        )

    def get_improper(
        self, type1: str, type2: str, type3: str, type4: str, model: str = None
    ) -> HarmonicImproperParams | DistanceImproperParams | None:
        """Get improper parameters for a quadruple of atom types."""
        return self._get_quadruple_params(
            self.improper, type1, type2, type3, type4, model
        )

    def get_model_names(self) -> list[str]:
        """Get list of available model names."""
        return list(self.models.keys())

    def get_model(self, name: str) -> ForceFieldModel | None:
        """Get model metadata by name."""
        return self.models.get(name)

    def get_atom_types_for_model(self, model: str) -> list[str]:
        """Get all atom type names for a specific model."""
        prefix = f"{model}."
        return [key for key in self.atom.keys() if key.startswith(prefix)]

    def clear(self) -> None:
        """Clear all loaded data."""
        self.atom.clear()
        self.lj.clear()
        self.buckingham.clear()
        self.bond.clear()
        self.angle.clear()
        self.improper.clear()
        self.bondbond.clear()
        self.bondangle.clear()
        self.angleangletorsion.clear()
        self.angleangle.clear()
        self.models.clear()
        self.dihedral.clear()

    def edit(self) -> None:
        """
        Launch the Streamlit interactive forcefield editor.
        """
        import subprocess

        editor_path = Path(__file__).parent / "_editor.py"

        if not editor_path.exists():
            raise FileNotFoundError(
                f"Editor script '_editor.py' not found alongside {__file__}."
            )

        print("Launching CEMD Forcefield Editor via Streamlit...")
        try:
            subprocess.run(["streamlit", "run", str(editor_path)], check=True)
        except FileNotFoundError:
            print("Error: 'streamlit' is not installed or not found in your PATH.")
        except subprocess.CalledProcessError as e:
            print(f"Streamlit editor exited with code {e.returncode}")

    def to_dataframes(self) -> dict[str, pd.DataFrame]:
        """Convert to pandas DataFrames for compatibility."""
        dfs = {}

        records = []
        for key, params in self.atom.items():
            short_type = key.split(".")[-1] if "." in key else key
            records.append(
                {
                    "type": short_type,
                    "full_type": key,
                    "element": params.element,
                    "charge": params.charge,
                    "model": params.model,
                    "environment": params.environment,
                    "ref": params.ref,
                    "mass": params.mass,
                }
            )
        dfs["atoms"] = pd.DataFrame(records)

        records = []
        for key, params in self.lj.items():
            pair = key.split(".")[-1] if "." in key else key
            records.append(
                {
                    "pair": pair,
                    "full_key": key,
                    "epsilon": params.epsilon,
                    "sigma": params.sigma,
                    "model": params.model,
                    "ref": params.ref,
                }
            )
        dfs["lj"] = pd.DataFrame(records)

        records = []
        for key, params in self.bond.items():
            bond_id = key.split(".")[-1] if "." in key else key
            if hasattr(params, "k") and hasattr(params, "r0"):
                records.append(
                    {
                        "bond": bond_id,
                        "full_key": key,
                        "type": "harmonic",
                        "k": params.k,
                        "r0": params.r0,
                        "model": params.model,
                        "ref": params.ref,
                    }
                )
            elif hasattr(params, "k2") and hasattr(params, "r0"):
                records.append(
                    {
                        "bond": bond_id,
                        "full_key": key,
                        "type": "class2",
                        "k2": params.k2,
                        "k3": params.k3,
                        "k4": params.k4,
                        "r0": params.r0,
                        "model": params.model,
                        "ref": params.ref,
                    }
                )
            elif hasattr(params, "D") and hasattr(params, "alpha"):
                records.append(
                    {
                        "bond": bond_id,
                        "full_key": key,
                        "type": "morse",
                        "r0": params.r0,
                        "D": params.D,
                        "alpha": params.alpha,
                        "model": params.model,
                        "ref": params.ref,
                    }
                )
        dfs["bonds"] = pd.DataFrame(records)

        records = []
        for key, params in self.buckingham.items():
            pair = key.split(".")[-1] if "." in key else key
            records.append(
                {
                    "pair": pair,
                    "full_key": key,
                    "a": getattr(params, "a", None),
                    "rho": getattr(params, "rho", None),
                    "c": getattr(params, "c", None),
                    "model": params.model,
                    "ref": params.ref,
                }
            )
        dfs["buckingham"] = pd.DataFrame(records)

        # Angles
        records = []
        for key, params in self.angle.items():
            angle_id = key.split(".")[-1] if "." in key else key
            rec = {
                "angle": angle_id,
                "full_key": key,
                "model": params.model,
                "ref": params.ref,
            }
            if hasattr(params, "k") and hasattr(params, "theta0"):
                rec.update({"type": "harmonic", "k": params.k, "theta0": params.theta0})
            elif hasattr(params, "k2") and hasattr(params, "theta0"):
                rec.update(
                    {
                        "type": "class2",
                        "theta0": params.theta0,
                        "k2": params.k2,
                        "k3": getattr(params, "k3", 0.0),
                        "k4": getattr(params, "k4", 0.0),
                    }
                )
            records.append(rec)
        dfs["angles"] = pd.DataFrame(records)

        # Dihedrals
        records = []
        for key, params in self.dihedral.items():
            dih_id = key.split(".")[-1] if "." in key else key
            records.append({"dihedral": dih_id, "full_key": key, "params": str(params)})
        dfs["dihedrals"] = pd.DataFrame(records)

        # Impropers
        records = []
        for key, params in self.improper.items():
            imp_id = key.split(".")[-1] if "." in key else key
            rec = {
                "improper": imp_id,
                "full_key": key,
                "model": getattr(params, "model", ""),
            }
            if hasattr(params, "k") and hasattr(params, "d") and hasattr(params, "n"):
                rec.update(
                    {"k": params.k, "d": params.d, "n": params.n, "type": "periodic"}
                )
            elif hasattr(params, "k"):
                rec.update(
                    {"k": getattr(params, "k"), "chi0": getattr(params, "chi0", 0.0)}
                )
            elif hasattr(params, "k2"):
                rec.update(
                    {"k2": getattr(params, "k2"), "k4": getattr(params, "k4", 0.0)}
                )
            records.append(rec)
        dfs["impropers"] = pd.DataFrame(records)

        return dfs

    def __repr__(self) -> str:
        return f"ForceFieldDatabase(models={list(self.models.keys())})"
