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

from __future__ import annotations

import warnings
from collections.abc import Sequence
from itertools import combinations, combinations_with_replacement
from typing import Any

import numpy as np

from ..forcefield.forcefield_database import ForceFieldDatabase
from ..forcefield.models import LJParams
from ._format import canonical_ff_type


class ForceFieldMixin:
    """Mixin for force field operations on AtomicSystem."""

    def set_ff_keys(
        self,
        atom: dict[str, str] | None = None,
        bond: dict[str, str] | None = None,
        angle: dict[str, str] | None = None,
        dihedral: dict[str, str] | None = None,
        improper: dict[str, str] | None = None,
        overwrite: bool = True,
    ) -> None:
        """Set or update force-field keys for atom and connectivity types.

        Parameters
        ----------
        atom, bond, angle, dihedral, improper
            Mapping between system types and force-field keys.
        overwrite
            If True, existing assignments are replaced.
        """

        mappings = {
            "atom": atom,
            "bond": bond,
            "angle": angle,
            "dihedral": dihedral,
            "improper": improper,
        }

        for kind, values in mappings.items():
            if not values:
                continue

            target = getattr(self._ff_keys, kind)

            if overwrite:
                target.update(values)
            else:
                for system_type, ff_key in values.items():
                    if system_type not in target:
                        target[system_type] = ff_key

    def set_ff_from_database(
        self,
        atom_assignments: dict[str, str] = None,
        bond_assignments: dict[str, str] = None,
        angle_assignments: dict[str, str] = None,
        dihedral_assignments: dict[str, str] = None,
        improper_assignments: dict[str, str] = None,
        overwrite: bool = True,
    ) -> None:

        db = ForceFieldDatabase()

        atom_assignments = atom_assignments or {}
        bond_assignments = bond_assignments or {}
        angle_assignments = angle_assignments or {}
        dihedral_assignments = dihedral_assignments or {}
        improper_assignments = improper_assignments or {}

        resolved_atom_assignments = self._resolve_atom_assignments(
            atom_assignments,
            db,
        )

        self.set_ff_keys(
            atom=resolved_atom_assignments,
            bond=bond_assignments,
            angle=angle_assignments,
            dihedral=dihedral_assignments,
            improper=improper_assignments,
            overwrite=overwrite,
        )

        self._update_masses_and_charges(db)
        self._set_pair_params_from_db(resolved_atom_assignments, db, overwrite)
        self._set_topology_params_from_db("bond", bond_assignments, db, overwrite)
        self._set_topology_params_from_db("angle", angle_assignments, db, overwrite)
        self._set_topology_params_from_db(
            "dihedral", dihedral_assignments, db, overwrite
        )
        self._set_topology_params_from_db(
            "improper", improper_assignments, db, overwrite
        )

        # Class2 cross terms: they apply to an existing angle/dihedral/
        # improper, so they piggyback on that category's ff-key assignments
        # rather than taking their own.
        self._set_topology_params_from_db(
            "bondbond",
            {},
            db,
            overwrite,
            types_attr="angle_types",
            ff_keys_attr="angle",
        )
        self._set_topology_params_from_db(
            "bondangle",
            {},
            db,
            overwrite,
            types_attr="angle_types",
            ff_keys_attr="angle",
        )
        self._set_topology_params_from_db(
            "angleangletorsion",
            {},
            db,
            overwrite,
            types_attr="dihedral_types",
            ff_keys_attr="dihedral",
        )
        self._set_topology_params_from_db(
            "angleangle",
            {},
            db,
            overwrite,
            types_attr="improper_types",
            ff_keys_attr="improper",
        )

    def set_atom_ff_keys(self, value: Sequence[str] | dict[str | int, str]) -> None:
        self._apply_topology_ff_keys("atom", value)

    def set_bond_ff_keys(self, value: Sequence[str] | dict[str, str]) -> None:
        self._apply_topology_ff_keys("bond", value)

    def set_angle_ff_keys(self, value: Sequence[str] | dict[str, str]) -> None:
        self._apply_topology_ff_keys("angle", value)

    def set_dihedral_ff_keys(self, value: Sequence[str] | dict[str, str]) -> None:
        self._apply_topology_ff_keys("dihedral", value)

    def set_improper_ff_keys(self, value: Sequence[str] | dict[str, str]) -> None:
        self._apply_topology_ff_keys("improper", value)

    def explore_ff_database(self, visible_rows: int = 20) -> None:
        """
        Interactive Forcefield Explorer using the prompt_toolkit UI pattern.

        Launches an interactive terminal user interface (TUI) that allows the user
        to browse and select force field parameters for each element and atom type
        present in the system. The explorer reads a force field database from an
        TOML files and presents matching parameter sets for user selection via
        keyboard navigation.

        Parameters
        ----------
        ff_database_dir : str, optional
            Path to the directory containing the force field database.
            The file must contain a sheet named 'list' with force field parameters.
            Defaults to get_ff_db_path() (global constant).

        visible_rows : int, optional
            Number of rows visible in the interactive table. Default is 20.

        Returns
        -------
        None
            The method updates the system's force field parameters in-place via
            successive calls to `set_ff_from_database()` as each selection is made.

        Raises
        ------
        FileNotFoundError
            If the database file does not exist.
        ValueError
            If required sheets are missing from the database.

        Notes
        -----
        The interactive UI provides:
            -Up/Down arrows: Navigate through parameter options
            -Enter/Space: Select the highlighted parameter
            -Ctrl+C/Escape/q: Cancel selection and skip element
            -p: Open the reference DOI in a web browser

        The function filters the database for parameters matching each element
        present in the system (self.elements). For each unique element, the user is
        presented with a list of available parameter sets to choose from.

        Examples
        --------
        >>> system = AtomicSystem(...)
        >>> system.explore_ff_database()
        Assigned: H -> clayff.h*
        Assigned: O -> clayff.o_star

        >>> # Using custom database
        >>> system.explore_ff_database('/path/to/database')
        """

        import webbrowser

        from prompt_toolkit import Application
        from prompt_toolkit.key_binding import KeyBindings
        from prompt_toolkit.layout import HSplit, Layout, ScrollOffsets, Window
        from prompt_toolkit.layout.controls import FormattedTextControl

        db = ForceFieldDatabase()
        dfs = db.to_dataframes()

        if "atoms" not in dfs or dfs["atoms"].empty:
            print("No atom types found in database.")
            return

        df_list = dfs["atoms"]

        assignments = {}

        # `.items()`, not `zip(dict, list)`: iterating the mapping yields its
        # keys, which are the atom *types*, so the lookup below searched the
        # database for element "Ow" or "Hw" and matched nothing. The explorer
        # only ever worked on a structure whose type names happened to be
        # element symbols. Types with no identifiable element are skipped
        # rather than silently misaligned.
        for t, el in self.elements.items():
            subset = df_list[df_list["element"] == el]
            if subset.empty:
                print(f"No parameters for element {el}. Skipping.")
                continue

            types = []
            for _, row in subset.iterrows():
                types.append(
                    {
                        "type": row["type"],
                        "full_type": row["full_type"],
                        "model": row["model"],
                        "environment": row["environment"],
                        "ref": row["ref"],
                    }
                )

            index = 0
            selected = None
            scroll_top = 0
            kb = KeyBindings()

            row_format = "{cursor} {type:<20} {model:<15} | {env:<45.45}"

            def render():
                header = row_format.format(
                    cursor=" ", type="TYPE", model="MODEL", env="ENVIRONMENT"
                )
                separator = "-" * len(header)

                lines = [
                    f"--- Select Forcefield Parameters for type {t} [{index + 1}/{len(types)}] ---",
                    f"Found {len(types)} options  [{index + 1}/{len(types)}]",
                    "↑↓ Move   [Enter] Assign   [p] Open DOI   [q] Quit",
                    "",
                    header,
                    separator,
                ]

                # Show only the visible window
                visible = types[scroll_top : scroll_top + visible_rows]

                for i, r in enumerate(visible):
                    real_index = scroll_top + i
                    cursor = "➜" if real_index == index else " "
                    lines.append(
                        row_format.format(
                            cursor=cursor,
                            type=r["type"],
                            model=r["model"],
                            env=r["environment"][:45] if r["environment"] else "",
                        )
                    )

                # Scroll indicator
                if len(types) > visible_rows:
                    lines.append(
                        f"\n  ↕ {scroll_top + 1}-{min(scroll_top + visible_rows, len(types))} of {len(types)}"
                    )

                return "\n".join(lines)

            text = FormattedTextControl(render)
            window = Window(
                content=text,
                always_hide_cursor=True,
                scroll_offsets=ScrollOffsets(top=1, bottom=1),
                allow_scroll_beyond_bottom=False,
                dont_extend_height=False,
            )

            @kb.add("up")
            def _(e):
                nonlocal index, scroll_top
                if index > 0:
                    index -= 1
                    if index < scroll_top:
                        scroll_top = index
                    e.app.invalidate()

            @kb.add("down")
            def _(e):
                nonlocal index, scroll_top
                if index < len(types) - 1:
                    index += 1
                    if index >= scroll_top + visible_rows:
                        scroll_top = index - visible_rows + 1
                    e.app.invalidate()

            @kb.add("enter")
            def _(e):
                nonlocal selected
                selected = types[index]
                e.app.exit()

            @kb.add("p")
            def _(e):
                ref = types[index].get("ref", "")
                if ref and ref.startswith("http"):
                    webbrowser.open(ref)

            @kb.add("q")
            @kb.add("escape")
            def _(e):
                e.app.exit()

            app = Application(
                layout=Layout(HSplit([window])),
                key_bindings=kb,
                full_screen=True,
            )

            app.run()

            if selected:
                assignments[t] = selected["full_type"]
                print(f"Assigned: {t} -> {selected['full_type']} ({selected['model']})")

                self.set_ff_from_database(atom_assignments=assignments)
            else:
                print(f"Selection cancelled for type {t}. Skipping.")
                continue

    def set_masses(self, value: Sequence[float] | dict[str | int, float]) -> None:
        """Set or update the atomic masses for the system and the atoms DataFrame."""
        self._apply_property_mapping("mass", value, self._masses)

    def set_charges(self, value: Sequence[float] | dict[str | int, float]) -> None:
        """Set or update the charges for the atomic system and the atoms DataFrame."""
        self._apply_property_mapping("charge", value, self._charges)

    def set_pair_params(
        self,
        atom_type1: str | int,
        atom_type2: str | int = None,
        params: Any = None,
    ) -> None:
        """Set pair-interaction parameters for two atom types."""

        if atom_type2 is None:
            atom_type2 = atom_type1

        if atom_type1 not in self.atom_types:
            raise ValueError(f"Atom type '{atom_type1}' does not exist in the system.")

        if atom_type2 not in self.atom_types:
            raise ValueError(f"Atom type '{atom_type2}' does not exist in the system.")

        pair_key = self._normalize_binary_key(atom_type1, atom_type2)
        self._ff_params.pair[pair_key] = params

    def _apply_topology_param(self, kind: str, topo_type: str, params: Any) -> None:
        """Store force-field parameters for a single bond/angle/dihedral/improper type.

        "bond", "angle" and "improper" type strings are canonicalized via
        ``canonical_ff_type`` first, matching the convention the system
        itself uses when generating them (see
        ``TopologyMixin._set_connection``) so a caller-supplied type in a
        different atom order still resolves to the right entry. Dihedrals
        keep their given (chain) order, since they are directional.
        """
        current_types = getattr(self, f"{kind}_types", [])

        if kind in ("bond", "angle", "improper"):
            canonical_type = canonical_ff_type(tuple(topo_type.split("-")), kind)
        else:
            canonical_type = topo_type

        if canonical_type not in current_types:
            warnings.warn(
                f"Target {kind} type '{topo_type}' not present in the system.",
                UserWarning,
            )

        getattr(self._ff_params, kind)[canonical_type] = params

    def set_bond_params(self, bond_type: str, params: Any) -> None:
        self._apply_topology_param("bond", bond_type, params)

    def set_angle_params(self, angle_type: str, params: Any) -> None:
        self._apply_topology_param("angle", angle_type, params)

    def set_dihedral_params(self, dihedral_type: str, params: Any) -> None:
        self._apply_topology_param("dihedral", dihedral_type, params)

    def set_improper_params(self, improper_type: str, params: Any) -> None:
        self._apply_topology_param("improper", improper_type, params)

    def apply_pair_mixing_rules(
        self,
        rule: str = "arithmetic",
        overwrite: bool = False,
    ) -> None:
        """Calculate missing cross-interaction parameters."""

        missing_self_params = []

        for atom_type in self.atom_types:
            key = self._normalize_binary_key(atom_type, atom_type)

            if key not in self._ff_params.pair:
                missing_self_params.append(atom_type)

        if missing_self_params:
            raise ValueError(
                "Missing parameters for self-interaction of the following types: "
                f"{missing_self_params}. "
                "You must set these parameters with 'set_pair_params' "
                "before blending."
            )

        for atom_type in self.atom_types:
            key = self._normalize_binary_key(atom_type, atom_type)
            params = self._ff_params.pair[key]

            if not isinstance(params, LJParams):
                raise ValueError(
                    f"Self-interaction for type '{atom_type}' is not LJ "
                    f"parameters (got {type(params).__name__}). "
                    "Mixing rules only apply to LJ parameters."
                )

        for atom_type1, atom_type2 in combinations(self.atom_types, 2):
            pair_key = self._normalize_binary_key(atom_type1, atom_type2)

            if pair_key in self._ff_params.pair and not overwrite:
                continue

            params1 = self._ff_params.pair[
                self._normalize_binary_key(atom_type1, atom_type1)
            ]
            params2 = self._ff_params.pair[
                self._normalize_binary_key(atom_type2, atom_type2)
            ]

            eps1, sig1 = params1.epsilon, params1.sigma
            eps2, sig2 = params2.epsilon, params2.sigma

            if rule == "arithmetic":
                eps_mixed = (eps1 * eps2) ** 0.5
                sig_mixed = (sig1 + sig2) / 2.0

            elif rule == "geometric":
                eps_mixed = (eps1 * eps2) ** 0.5
                sig_mixed = (sig1 * sig2) ** 0.5

            else:
                raise ValueError(
                    f"Unknown mixing rule: '{rule}'. "
                    "Supported rules: 'arithmetic', 'geometric'."
                )

            self._ff_params.pair[pair_key] = LJParams(
                epsilon=eps_mixed,
                sigma=sig_mixed,
                ref="mixed",
                model="mixed",
            )

    def _apply_property_mapping(
        self,
        property_name: str,
        value: dict[str | int, float],
        storage: dict[str | int, float],
    ) -> None:
        """Validate and apply a mass or charge mapping to the given storage dict."""
        if not isinstance(value, dict):
            raise TypeError(
                f"Unsupported argument type for {property_name}: "
                f"{type(value).__name__}. Expected a dict."
            )

        # Atom types are always exposed as strings (see `atom_types`); keep
        # dict keys consistent so a numeric key (e.g. LAMMPS-style type 1)
        # actually matches instead of silently missing.
        value = {str(k): v for k, v in value.items()}

        current_types = set(self.atom_types)
        missing_types = [t for t in value if t not in current_types]
        if missing_types:
            warnings.warn(
                f"Targeted types missing from current system for {property_name}: {missing_types}",
                UserWarning,
            )

        if property_name == "mass":
            try:
                invalid = any(float(m) <= 0 for m in value.values())
            except (TypeError, ValueError) as e:
                raise ValueError(
                    f"Invalid mass value detected. All masses must be numbers. Details: {e}"
                ) from e
            if invalid:
                raise ValueError("All masses must be strictly positive (> 0).")

        storage.update(value)

        # Charges live per-atom in the atoms DataFrame; masses are per-type only.
        if property_name == "charge" and getattr(self, "atoms", None) is not None:
            mapped_values = self.atoms["type"].map(storage)

            if "charge" in self.atoms.columns:
                self.atoms["charge"] = mapped_values.fillna(self.atoms["charge"])
            else:
                self.atoms["charge"] = mapped_values

            unassigned = self.atoms.loc[self.atoms["charge"].isna(), "type"].unique()
            if len(unassigned):
                warnings.warn(
                    f"Certain types of atoms have no assigned charge: {list(unassigned)}",
                    UserWarning,
                )

    def _apply_topology_ff_keys(
        self, kind: str, value: Sequence[str] | dict[str, str]
    ) -> None:
        """Helper method to set or update force field keys for a given topology kind.

        "atom" and "dihedral" keys are stored as given (an atom type is a
        single name; a dihedral is a directional chain). "bond", "angle"
        and "improper" keys are canonicalized via ``canonical_ff_type`` so
        they match the type strings the system itself generates (see
        ``TopologyMixin._set_connection``), regardless of the atom order
        the caller used.
        """
        current_types_list = getattr(self, f"{kind}_types", [])
        current_types = set(current_types_list)
        target_dict = getattr(self._ff_keys, kind)

        def canonicalize(topo_type: str) -> str:
            if kind in ("bond", "angle", "improper"):
                return canonical_ff_type(tuple(topo_type.split("-")), kind)
            return topo_type

        if isinstance(value, dict):
            canonical_value = {canonicalize(t): ff_key for t, ff_key in value.items()}

            missing_types = [t for t in canonical_value if t not in current_types]
            if missing_types:
                warnings.warn(
                    f"Target {kind} types not present: {missing_types}", UserWarning
                )

            target_dict.update(canonical_value)

        elif isinstance(value, (list, np.ndarray, tuple)):
            if len(value) != len(current_types_list):
                raise ValueError(
                    f"Sequence length ({len(value)}) does not match {kind} types count ({len(current_types_list)})."
                )
            target_dict.update(dict(zip(current_types_list, value)))
        else:
            raise TypeError(f"Unsupported argument type: {type(value).__name__}.")

    def _resolve_pair_parameters(
        self,
        db: ForceFieldDatabase,
        ff_t1: str,
        ff_t2: str,
        orig_t1: str,
        orig_t2: str,
        is_self_interaction: bool,
    ) -> Any:
        """Helper method to find or derive pair parameters from the database."""

        # Search with exact keys (ex: db.lj or db.buckingham)
        params = db.get_lj(ff_t1, ff_t2) or db.get_buckingham(ff_t1, ff_t2)
        if params is not None:
            return params

        # Search with the original keys (without the suffix)
        params = db.get_lj(orig_t1, orig_t2) or db.get_buckingham(orig_t1, orig_t2)
        if params is not None:
            return params

        # If not found, try to derive it from each atom's own
        # self-interaction LJ parameters (Lorentz-Berthelot mixing).
        self1 = db.get_lj(ff_t1, ff_t1) or db.get_lj(orig_t1, orig_t1)
        self2 = db.get_lj(ff_t2, ff_t2) or db.get_lj(orig_t2, orig_t2)

        if self1 is not None and self2 is not None:
            eps_mixed = (self1.epsilon * self2.epsilon) ** 0.5
            sig_mixed = (self1.sigma + self2.sigma) / 2.0
            return LJParams(
                epsilon=eps_mixed,
                sigma=sig_mixed,
                ref="mixed from self-interaction",
                model="mixed",
            )

        # Fallback: if still not found and it is a self-interaction, default to zero
        if is_self_interaction:
            return LJParams(
                epsilon=0.0,
                sigma=0.0,
                ref="zero (not found in DB)",
                model="gromos",
            )

        return None

    def _set_pair_params_from_db(
        self, atom_assignments: dict, db: ForceFieldDatabase, overwrite: bool
    ) -> None:
        """Set pair parameters by searching the database or applying mixing rules."""
        if not self._ff_keys.atom:
            return

        missing_pairs = []

        # Using a "dictionary comprehension" to do this in 1 line
        original_types = {atype: atype.split("_")[0] for atype in self._ff_keys.atom}

        for atom_type1, atom_type2 in combinations_with_replacement(
            self._ff_keys.atom, 2
        ):
            pair_key = self._normalize_binary_key(atom_type1, atom_type2)

            if pair_key in self._ff_params.pair and not overwrite:
                continue

            ff_t1 = self._ff_keys.atom[atom_type1]
            ff_t2 = self._ff_keys.atom[atom_type2]

            orig_t1 = original_types[atom_type1]
            orig_t2 = original_types[atom_type2]

            is_self_interaction = atom_type1 == atom_type2

            # We call our new method tool for doing research
            params = self._resolve_pair_parameters(
                db, ff_t1, ff_t2, orig_t1, orig_t2, is_self_interaction
            )

            if params is not None:
                self._ff_params.pair[pair_key] = params
            elif is_self_interaction:
                missing_pairs.append((atom_type1, ff_t1))

        if missing_pairs:
            warnings.warn(
                f"Missing pair parameters for atomic types:{missing_pairs}",
                category=UserWarning,
                stacklevel=2,
            )

    def _set_topology_params_from_db(
        self,
        topo_type: str,
        assignments: dict,
        db: ForceFieldDatabase,
        overwrite: bool,
        types_attr: str | None = None,
        ff_keys_attr: str | None = None,
    ) -> None:
        """Generic method to retrieve topology parameters (bonds, angles, etc.) from the DB.

        Parameters
        ----------
        types_attr, ff_keys_attr
            By default, both are derived from ``topo_type`` (e.g. "bond" ->
            ``self.bond_types`` / ``self._ff_keys.bond``). Class2 cross
            terms (bondbond, bondangle, angleangletorsion, angleangle) have
            no connectivity or ff-key category of their own: they apply to
            an existing angle/dihedral/improper, so they reuse that
            category's type list and ff-key assignments instead.
        """
        missing_items = []

        types_list = getattr(self, types_attr or f"{topo_type}_types", [])
        ff_keys_dict = getattr(self._ff_keys, ff_keys_attr or topo_type)
        db_attr = getattr(db, topo_type, None)
        db_get_method = getattr(db, f"get_{topo_type}", None)
        params_storage = getattr(self._ff_params, topo_type)

        for topo_str in types_list:
            params = None

            ff_key = ff_keys_dict.get(topo_str)
            if ff_key and db_attr:
                params = db_attr.get(ff_key)

            if params is None and topo_str in assignments:
                target_key = assignments[topo_str]
                if db_attr:
                    params = db_attr.get(target_key)

            if params is None and db_get_method:
                elements = [e.strip() for e in topo_str.split("-")]
                ff_types = [self._get_ff_key_for_type(e) for e in elements]

                if all(ff_types):
                    params = db_get_method(*ff_types)

            if params is not None:
                if topo_str not in params_storage or overwrite:
                    params_storage[topo_str] = params
            else:
                missing_items.append(topo_str)

        if missing_items:
            fr_names = {
                "bond": "liaisons",
                "angle": "d'angles",
                "dihedral": "de dièdres",
                "improper": "impropres",
            }
            warnings.warn(
                f"Parameters {fr_names.get(topo_type, topo_type)} not found for: {missing_items}",
                category=UserWarning,
                stacklevel=2,
            )

    def _resolve_atom_assignments(
        self, assignments: dict, db: ForceFieldDatabase
    ) -> dict:
        """Resolves atom abbreviations (ex: 'h*' -> 'clayff.h*')."""
        resolved = {}
        for sys_type, db_type in assignments.items():
            if "." in db_type:
                resolved[sys_type] = db_type
            else:
                found = next(
                    (
                        full_type
                        for full_type in db.atom.keys()
                        if full_type.endswith(f".{db_type}")
                    ),
                    None,
                )

                if found is None:
                    warnings.warn(
                        f"Atom type '{db_type}' not found in database.",
                        category=UserWarning,
                        stacklevel=3,
                    )
                    continue

                resolved[sys_type] = found
        return resolved

    @staticmethod
    def _normalize_binary_key(type1: str | int, type2: str | int) -> str:
        """Return a canonical key for a two-type interaction."""
        return "-".join(sorted((str(type1), str(type2))))

    def _update_masses_and_charges(self, db: ForceFieldDatabase) -> None:
        """Updates masses and charges from atoms and their ff_key."""

        masses_update = {}
        charges_update = {}
        missing_types = []

        for sys_type, db_type in self.ff_keys.atom.items():
            if db_type is None:
                continue

            atom_type = db.get_atom_type(db_type)

            if atom_type:
                if atom_type.mass is not None:
                    masses_update[sys_type] = atom_type.mass
                # Guarded like the mass: a force field that defines no
                # per-type charge (GROMOS, the CHARMM Interface files)
                # reports None, and applying that as 0.0 silently wiped the
                # partial charges the system already carried -- a molecule
                # solvated from an ATB topology lost all its electrostatics
                # the moment its parameters were looked up.
                if atom_type.charge is not None:
                    charges_update[sys_type] = atom_type.charge
            else:
                missing_types.append((sys_type, db_type))

        if missing_types:
            warnings.warn(
                f"Types introuvables pour les masses/charges : {missing_types}",
                category=UserWarning,
                stacklevel=2,
            )

        if masses_update:
            self.set_masses(masses_update)
        if charges_update:
            self.set_charges(charges_update)

    def _get_ff_key_for_type(self, sys_type: str) -> str | None:
        if hasattr(self._ff_keys, "atom") and self._ff_keys.atom:
            return self._ff_keys.atom.get(sys_type)
        return None

    @staticmethod
    def _remap_interaction_key(
        interaction_type: str, kind: str, type_mapping: dict[str, str]
    ) -> str:
        """Rebuild a canonical interaction key after atom types are renamed."""
        old_types = interaction_type.split("-")
        new_types = tuple(
            type_mapping.get(atom_type, atom_type) for atom_type in old_types
        )

        # For dihedrals, keep the original (chain) order.
        if kind == "dihedral":
            return "-".join(new_types)

        return canonical_ff_type(new_types, kind)

    def _remap_ff_keys(self, type_mapping: dict[str, str]) -> None:
        """Remap force-field keys and already-resolved parameters after
        atom types have been renamed.

        Both ``_ff_keys`` (type -> database key) and ``_ff_params``
        (type -> resolved parameter object, including the class2 cross
        terms) are keyed by interaction-type strings built from atom type
        names. Without this remap, a rename would silently orphan every
        already-resolved parameter: its key would no longer match any
        current bond/angle/dihedral/improper type, so it would simply
        disappear from the active (filtered) view used when writing files.
        """

        # Atom FF keys
        old_atom_keys = dict(self._ff_keys.atom)
        new_atom_keys = {}
        for old_type, ff_key in old_atom_keys.items():
            new_type = type_mapping.get(old_type, old_type)
            new_atom_keys[new_type] = ff_key
        self._ff_keys.atom = new_atom_keys

        # Interaction FF keys
        for kind in ["bond", "angle", "dihedral", "improper"]:
            old_keys = dict(getattr(self._ff_keys, kind))
            new_keys = {
                self._remap_interaction_key(interaction_type, kind, type_mapping): (
                    ff_key
                )
                for interaction_type, ff_key in old_keys.items()
            }
            setattr(self._ff_keys, kind, new_keys)

        # Pair parameters (keyed like a bond: two atom types, symmetric).
        old_pairs = dict(self._ff_params.pair)
        self._ff_params.pair = {
            self._normalize_binary_key(
                *(type_mapping.get(t, t) for t in pair_type.split("-"))
            ): params
            for pair_type, params in old_pairs.items()
        }

        # Interaction parameters, including the class2 cross terms, which
        # piggyback on their parent category's key convention (bondbond and
        # bondangle apply to an angle, angleangletorsion to a dihedral,
        # angleangle to an improper).
        param_kind_for = {
            "bond": "bond",
            "angle": "angle",
            "dihedral": "dihedral",
            "improper": "improper",
            "bondbond": "angle",
            "bondangle": "angle",
            "angleangletorsion": "dihedral",
            "angleangle": "improper",
        }
        for param_attr, kind in param_kind_for.items():
            old_params = dict(getattr(self._ff_params, param_attr))
            new_params = {
                self._remap_interaction_key(interaction_type, kind, type_mapping): (
                    params
                )
                for interaction_type, params in old_params.items()
            }
            setattr(self._ff_params, param_attr, new_params)

    def _validate_and_convert_coeffs(self, coeffs, expected_len, name):
        """Validates and converts a list of coefficients."""
        if len(coeffs) != expected_len:
            raise ValueError(
                f"Invalid number of coefficients for {name}: "
                f"expected {expected_len}, got {len(coeffs)}."
            )
        try:
            return [float(c) for c in coeffs]
        except (TypeError, ValueError) as e:
            raise ValueError(
                f"Invalid coefficient values for {name}: expected numbers, got {coeffs}"
            ) from e
