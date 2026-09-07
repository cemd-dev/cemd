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

import copy
import warnings
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Self

import numpy as np
import pandas as pd

from .._constants import AVOGADRO, INV_MASSES, MASS_KEYS, MASSES_DICT
from ._edit import EditMixin
from ._format import (
    ANGLES_COLUMNS,
    ATOMS_COLUMNS,
    BONDS_COLUMNS,
    DIHEDRALS_COLUMNS,
    IMPROPERS_COLUMNS,
    VELOCITIES_COLUMNS,
    normalize_dataframe,
    normalize_property_to_dict,
)
from ._io import IOMixin
from ._view import view
from .forcefield_classes import ForceFieldKeys, ForceFieldParams
from .forcefield_mixin import ForceFieldMixin
from .topology_mixin import TopologyMixin

if TYPE_CHECKING:
    from ..build import SolutionBuilder


class AtomicSystem(EditMixin, IOMixin, TopologyMixin, ForceFieldMixin):
    """
    Container storing the complete atomic system representation.

    The class stores atomic coordinates, topology information
    (bonds, angles, dihedrals and impropers), simulation box
    parameters, atomic masses and charges, and force field parameters.

    It is designed to represent the content of a LAMMPS data file and
    provides utilities for editing, analysing and exporting atomistic
    systems.

    Parameters
    ----------
    system_dict : dict[str, Any]
        Dictionary containing the complete system data. See
        ``cemd.core._format.SYSTEM_DICT_FORMAT`` for the authoritative,
        machine-checked list of required and optional keys; in summary:

        Required keys:
            - atoms : pandas.DataFrame
                Per-atom data (type, charge, x, y, z), indexed by atom id.
            - box : sequence
                Simulation box, in any format accepted by ``set_box()``.
            - masses, charges : dict
                Per-atom-type masses and charges.

        Optional keys:
            - bonds, angles, dihedrals, impropers : pandas.DataFrame or None
                Topology tables (atom_1..atom_N columns plus a "type" column).
            - velocities : pandas.DataFrame or None
                Atomic velocities (vx, vy, vz).
            - atom_style : str
                LAMMPS atom style (default "full").
            - Force-field keys/params, as produced by
              ``ForceFieldKeys.to_system_dict()`` and
              ``ForceFieldParams.to_system_dict()``.

    Attributes
    ----------
    atoms : pandas.DataFrame
        Atomic information including coordinates, atom types and charges.
    bonds : pandas.DataFrame or None
        Bond topology information.
    angles : pandas.DataFrame or None
        Angle topology information.
    dihedrals : pandas.DataFrame or None
        Dihedral topology information.
    impropers : pandas.DataFrame or None
        Improper topology information.
    velocities : pandas.DataFrame or None
        Atomic velocities when available.
    forcefield_keys : ForceFieldKeys
        Centralized mapping between system types (atoms, bonds, etc.) and database keys.
    forcefield_params : ForceFieldParams
        Centralized storage for all force field parameters (pair, bond, angle, dihedral, improper).
    """

    _atoms: pd.DataFrame
    _bonds: pd.DataFrame | None
    _angles: pd.DataFrame | None
    _dihedrals: pd.DataFrame | None
    _impropers: pd.DataFrame | None
    _velocities: pd.DataFrame | None

    _masses: dict[str | float]
    _charges: dict[str | float]
    _ff_keys: ForceFieldKeys
    _ff_params: ForceFieldParams
    _atom_style: str

    _box_lmp: tuple[
        tuple[float, float],
        tuple[float, float],
        tuple[float, float],
        tuple[float, float, float],
    ]
    _box: np.ndarray
    _box_vectors: tuple[np.ndarray]

    def __init__(self, system_dict: dict[str, Any]) -> None:
        """
        Initialize an atomic system from a dictionnary.

        Parameters
        ----------
        system_dict : dict[str, Any]
            Dictionary containing all atomic, topological and simulation
            box information.
        """
        self._from_system_dict(system_dict)
        self._finalize_data()

    def __repr__(self) -> str:
        """Return a concise string representation of the AtomicSystem."""
        sections = [
            ("bonds", self.num_bonds),
            ("angles", self.num_angles),
            ("dihedrals", self.num_dihedrals),
            ("impropers", self.num_impropers),
        ]

        interactions = ", ".join(
            f"{count} {name}" for name, count in sections if count > 0
        )

        if interactions:
            return f"<AtomicSystem with {self.num_atoms} atoms, {interactions}>"

        return f"<AtomicSystem with {self.num_atoms} atoms>"

    def __str__(self) -> str:
        """
        Return the string representation of the object.

        Returns
        -------
        str
            Equivalent to __repr__.
        """
        return self.__repr__()

    def _from_system_dict(self, system_dict: dict[str, Any]) -> None:
        """Helper to map the system dictionary to class attributes."""

        self._atoms = system_dict["atoms"]
        self._bonds = system_dict.get("bonds")
        self._angles = system_dict.get("angles")
        self._dihedrals = system_dict.get("dihedrals")
        self._impropers = system_dict.get("impropers")
        self._velocities = system_dict.get("velocities")

        self.set_box(system_dict["box"])

        self._masses = normalize_property_to_dict(
            dict(system_dict["masses"]), [str(t) for t in self.atoms.type.unique()]
        )
        self._charges = normalize_property_to_dict(
            dict(system_dict["charges"]), [str(t) for t in self.atoms.type.unique()]
        )
        self._atom_style = system_dict.get("atom_style", "full")

        self._pmg_struct = system_dict.get("_pmg_struct", None)

        self._ff_keys = ForceFieldKeys.from_system_dict(system_dict)
        self._ff_params = ForceFieldParams.from_system_dict(system_dict)

    def _finalize_data(self) -> None:
        """Sorts indices and ensures integer types for atom references."""

        self._atoms = normalize_dataframe(
            self._atoms,
            ATOMS_COLUMNS,
            "atoms",
        )

        self._bonds = normalize_dataframe(
            self._bonds,
            BONDS_COLUMNS,
            "bonds",
        )

        self._angles = normalize_dataframe(
            self._angles,
            ANGLES_COLUMNS,
            "angles",
        )

        self._dihedrals = normalize_dataframe(
            self._dihedrals,
            DIHEDRALS_COLUMNS,
            "dihedrals",
        )

        self._impropers = normalize_dataframe(
            self._impropers,
            IMPROPERS_COLUMNS,
            "impropers",
        )

        self._velocities = normalize_dataframe(
            self._velocities,
            VELOCITIES_COLUMNS,
            "velocities",
        )

        for df in [
            self._atoms,
            self._velocities,
            self._bonds,
            self._angles,
            self._dihedrals,
            self._impropers,
        ]:
            if df is not None:
                df.sort_index(inplace=True)

        topology_dfs = {"bonds": 2, "angles": 3, "dihedrals": 4, "impropers": 4}
        for name, n_atoms in topology_dfs.items():
            df = getattr(self, name)
            if df is not None:
                cols = [f"atom_{i}" for i in range(1, n_atoms + 1)]
                df[cols] = df[cols].astype(int)

    def _replace_internals(self, other: AtomicSystem) -> None:
        """Replace internal content with that of another system."""

        self._atoms = other._atoms.copy()
        self._bonds = other._bonds.copy() if other._bonds is not None else None
        self._angles = other._angles.copy() if other._angles is not None else None
        self._dihedrals = (
            other._dihedrals.copy() if other._dihedrals is not None else None
        )
        self._impropers = (
            other._impropers.copy() if other._impropers is not None else None
        )
        self._velocities = (
            other._velocities.copy() if other._velocities is not None else None
        )

        # Through `set_box`, so that the three internal representations
        # (lattice, LAMMPS bounds, vectors) stay consistent: assigning
        # `_box_lmp` alone left `.box` and `.volume` reporting the
        # *previous* cell after add_layer()/add_droplet()/set_topology().
        self.set_box(other._box_lmp)
        self._masses = dict(other._masses)
        self._charges = dict(other._charges)
        self._atom_style = other._atom_style

        self._ff_params = copy.deepcopy(other._ff_params)
        self._ff_keys = copy.deepcopy(other._ff_keys)

        if other._pmg_struct is not None:
            self._pmg_struct = copy.deepcopy(other._pmg_struct)
        else:
            self._pmg_struct = None

    def _to_system_dict(self) -> dict:
        """Serialize the system to a dictionary for reconstruction."""
        return {
            "atoms": self._atoms.copy(),
            "bonds": self._bonds.copy() if self._bonds is not None else None,
            "angles": self._angles.copy() if self._angles is not None else None,
            "dihedrals": self._dihedrals.copy()
            if self._dihedrals is not None
            else None,
            "impropers": self._impropers.copy()
            if self._impropers is not None
            else None,
            "velocities": self._velocities.copy()
            if self._velocities is not None
            else None,
            "box": self._box_lmp,
            "masses": dict(self._masses),
            "charges": dict(self._charges),
            "atom_style": self._atom_style,
            **self._ff_keys.to_system_dict(),
            **self._ff_params.to_system_dict(),
        }

    def copy(self) -> Self:
        """
        Create an independent deep copy of the atomic system.

        The copied object contains duplicated atomic coordinates,
        topology tables, simulation box information and force field
        parameters.
        """
        new = self.__class__.__new__(self.__class__)

        new._atoms = self._atoms.copy()
        new._bonds = self._bonds.copy() if self._bonds is not None else None
        new._angles = self._angles.copy() if self._angles is not None else None
        new._dihedrals = self._dihedrals.copy() if self._dihedrals is not None else None
        new._impropers = self._impropers.copy() if self._impropers is not None else None
        new._velocities = (
            self._velocities.copy() if self._velocities is not None else None
        )

        new._box_lmp = self._box_lmp
        new._masses = dict(self._masses)
        new._charges = dict(self._charges)
        new._atom_style = self._atom_style
        new._ff_params = copy.deepcopy(self._ff_params)
        new._ff_keys = copy.deepcopy(self._ff_keys)

        new._pmg_struct = (
            copy.deepcopy(self._pmg_struct) if self._pmg_struct is not None else None
        )

        return new

    @property
    def atoms(self) -> pd.DataFrame:
        return self._atoms

    @atoms.setter
    def atoms(self, value: pd.DataFrame) -> None:
        self._atoms = normalize_dataframe(value, ATOMS_COLUMNS, "atoms")

    @property
    def bonds(self) -> pd.DataFrame:
        return self._bonds

    @bonds.setter
    def bonds(self, value: pd.DataFrame) -> None:
        self._bonds = normalize_dataframe(value, BONDS_COLUMNS, "bonds")

    @property
    def angles(self) -> pd.DataFrame:
        return self._angles

    @angles.setter
    def angles(self, value: pd.DataFrame) -> None:
        self._angles = normalize_dataframe(value, ANGLES_COLUMNS, "angles")

    @property
    def dihedrals(self) -> pd.DataFrame:
        return self._dihedrals

    @dihedrals.setter
    def dihedrals(self, value: pd.DataFrame) -> None:
        self._dihedrals = normalize_dataframe(value, DIHEDRALS_COLUMNS, "dihedrals")

    @property
    def impropers(self) -> pd.DataFrame:
        return self._impropers

    @impropers.setter
    def impropers(self, value: pd.DataFrame) -> None:
        self._impropers = normalize_dataframe(value, IMPROPERS_COLUMNS, "dihedrals")

    @property
    def velocities(self) -> pd.DataFrame:
        return self._velocities

    @velocities.setter
    def velocities(self, value: pd.DataFrame) -> None:
        self._velocities = normalize_dataframe(value, VELOCITIES_COLUMNS, "velocities")

    @property
    def box(self) -> np.ndarray:
        """
        Lattice parameters [a, b, c, alpha, beta, gamma].

        Read-only — use set_box() to modify.

        Returns
        -------
        numpy.ndarray
            The box lattice parameters.
        """
        return self._box

    @property
    def volume(self) -> float:
        """
        Return the volume of the simulation box.

        Returns
        -------
        float
            Box volume in Å³.
        """
        v1, v2, v3 = self._box_vectors
        return np.dot(v1, np.cross(v2, v3))

    @property
    def masses(self) -> dict[str | int, float]:
        """Return atomic masses associated with atom types (read-only).

        Returns
        -------
        dict[str | int, float]
            Dictionary matching atom type ID to its mass.
        """
        return MappingProxyType(
            {
                t: float(self._masses.get(t, MASSES_DICT.get(t, 1.0)))
                for t in self.atom_types
            }
        )

    @property
    def charges(self) -> dict[str | int, float]:
        """
        Return atomic charges associated with atom types (read-only).

        Returns
        -------
        dict[str | int, float]
            Dictionary matching atom type ID to its charge.
        """
        return MappingProxyType(
            {t: float(self._charges.get(t, 0)) for t in self.atom_types}
        )

    #: How far an atom type's mass may sit from a real element's before the
    #: type is judged not to be that element. Every genuine element in the
    #: bundled force fields matches to better than 0.05 amu, while a GROMOS
    #: united atom -- a carbon carrying its apolar hydrogens, CH3 at 15.035
    #: -- sits about 1 amu away from anything real. The gap is wide.
    ELEMENT_MASS_TOLERANCE = 0.2

    @property
    def elements(self) -> dict[str | int, str]:
        """Return a mapping of atom types to their elemental symbols (read-only).

        The element is inferred from the atom type's mass. A type whose mass
        matches no element within
        :data:`ELEMENT_MASS_TOLERANCE` is left out of the mapping rather
        than assigned the nearest entry: an element is not something this
        can always know.

        That happens in two situations. A **united-atom** force field folds
        apolar hydrogens into their carbon, so GROMOS ``CH3`` weighs 15.035
        and the nearest element is oxygen -- the information needed to say
        "carbon" is simply not in the mass. And an element absent from the
        internal table (which covers hydrogen through barium) would
        otherwise be reported as the nearest one that is present.

        Returns
        -------
        dict[str | int, str]
            Atom type to element symbol, for the types whose mass identifies
            one. Callers should treat a missing key as "unknown" -- see
            :meth:`~cemd.core.topology_mixin.TopologyMixin.set_types_from_elements`,
            which warns and keeps the original type.
        """
        element_dict: dict[str | int, str] = {}
        unmatched: list[tuple[str | int, float]] = []

        for t, mass_val in self.masses.items():
            if mass_val is None or np.isnan(mass_val):
                continue

            distances = np.abs(MASS_KEYS - mass_val)
            index = distances.argmin()

            if distances[index] > self.ELEMENT_MASS_TOLERANCE:
                unmatched.append((t, float(mass_val)))
                continue

            element_dict[t] = str(INV_MASSES[MASS_KEYS[index]])

        if unmatched:
            detail = ", ".join(f"{t} ({mass:g} amu)" for t, mass in unmatched)
            warnings.warn(
                f"No element matches the mass of: {detail}. These types are "
                "absent from `elements`. A united-atom force field (GROMOS "
                "CH1-CH4) and any element heavier than barium both land here.",
                UserWarning,
                stacklevel=2,
            )

        return MappingProxyType(element_dict)

    @property
    def ff_keys(self) -> ForceFieldKeys:
        """Return force-field keys for active topology types."""
        return self._ff_keys.active_for(self)

    @property
    def ff_params(self) -> ForceFieldParams:
        """Return force-field parameters for active topology types."""
        return self._ff_params.active_for(self)

    @property
    def atom_types(self) -> list[str | int]:
        """
        Return the list of unique atom types.

        Returns
        -------
        list of str or int
            Sorted list of atom types.
        """
        return sorted([str(t) for t in self.atoms.type.unique()])

    @property
    def bond_types(self) -> list[str | int]:
        """
        Return the list of unique bond types.

        Returns
        -------
        list of str or int
            Sorted list of bond types.
        """
        if self.bonds is not None:
            return sorted(self.bonds.type.unique().tolist())
        else:
            return list()

    @property
    def angle_types(self) -> list[str | int]:
        """
        Return the list of unique angle types.

        Returns
        -------
        list of str or int
            Sorted list of angle types.
        """
        if self.angles is not None:
            return sorted(self.angles.type.unique().tolist())
        else:
            return list()

    @property
    def dihedral_types(self) -> list[str | int]:
        """
        Return the list of unique dihedral types.

        Returns
        -------
        list of str or int
            Sorted list of dihedral types.
        """
        if self.dihedrals is not None:
            return sorted(self.dihedrals.type.unique().tolist())
        else:
            return list()

    @property
    def improper_types(self) -> list[str | int]:
        """
        Return the list of unique improper types.

        Returns
        -------
        list of str or int
            Sorted list of improper types.
        """
        if self.impropers is not None:
            return sorted(self.impropers.type.unique().tolist())
        else:
            return list()

    @property
    def num_atoms(self) -> int:
        """
        Return the number of atoms.

        Returns
        -------
        int
            Number of atoms in the system.
        """
        return len(self.atoms)

    @property
    def num_bonds(self) -> int:
        """
        Return the number of bonds.

        Returns
        -------
        int
            Total number of bonds.
        """
        if self.bonds is None:
            return 0
        else:
            return len(self.bonds)

    @property
    def num_angles(self) -> int:
        """
        Return the number of angles.

        Returns
        -------
        int
            Total number of angles.
        """
        if self.angles is None:
            return 0
        else:
            return len(self.angles)

    @property
    def num_dihedrals(self) -> int:
        """
        Return the number of dihedrals.

        Returns
        -------
        int
            Total number of dihedrals.
        """
        if self.dihedrals is None:
            return 0
        else:
            return len(self.dihedrals)

    @property
    def num_impropers(self) -> int:
        """
        Return the number of impropers.

        Returns
        -------
        int
            Total number of impropers.
        """
        if self.impropers is None:
            return 0
        else:
            return len(self.impropers)

    @property
    def num_atom_types(self) -> int:
        """
        Return the number of atom types.

        Returns
        -------
        int
            Count of distinct atom types.
        """
        return len(self.atom_types)

    @property
    def num_bond_types(self) -> int:
        """
        Return the number of bond types.

        Returns
        -------
        int
            Count of distinct bond types.
        """
        return len(self.bond_types)

    @property
    def num_angle_types(self) -> int:
        """
        Return the number of angle types.

        Returns
        -------
        int
            Count of distinct angle types.
        """
        return len(self.angle_types)

    @property
    def num_dihedral_types(self) -> int:
        """
        Return the number of dihedral types.

        Returns
        -------
        int
            Count of distinct dihedral types.
        """
        return len(self.dihedral_types)

    @property
    def num_improper_types(self) -> int:
        """
        Return the number of improper types.

        Returns
        -------
        int
            Count of distinct improper types.
        """
        return len(self.improper_types)

    @property
    def total_charge(self) -> float:
        """
        Return the total system charge.

        Returns
        -------
        float
            Sum of all atomic charges.
        """
        # Summed over the per-atom `charge` column rather than rebuilt from
        # the per-type `charges` mapping: a force field can assign a charge
        # per atom rather than per type (ReaxFF's EEM charges, ATB/GROMOS
        # partial charges), and collapsing those to one representative value
        # per type reported wildly wrong totals -- a neutral ReaxFF-relaxed
        # C-S-H came back as +69 e. `set_charges` keeps this column in sync
        # with `_charges`, so per-type force fields are unaffected.
        if "charge" not in self.atoms:
            return 0.0

        return float(self.atoms["charge"].sum())

    @property
    def total_mass(self) -> float:
        """
        Return the total system mass.

        Returns
        -------
        float
            Total mass calculated from atomic types and counts.
        """
        # See the comment in `total_charge` about the str/raw-dtype mismatch.
        counts = self.atoms["type"].astype(str).value_counts()
        total_mass = sum(
            counts.get(atype, 0) * self.masses.get(atype, 0.0)
            for atype in self.atom_types
        )
        return total_mass

    @property
    def density(self) -> float:
        """
        Return the system density.

        Returns
        -------
        float
            Density in g/cm³.
        """
        return self.total_mass / AVOGADRO / self.volume / 1e-24

    def _get_type_summary(self) -> pd.DataFrame:
        """
        Return a summarized DataFrame of atom types, numbers, masses and charges.
        """
        df_atoms = self.atoms.copy()
        # Keep types as strings, matching `atom_types`/`masses`/`charges`,
        # so the `_masses` lookup below doesn't miss on a dtype mismatch
        # (e.g. numeric LAMMPS types).
        df_atoms["type"] = df_atoms["type"].astype(str)

        df_atoms["number"] = df_atoms.groupby("type")["type"].transform("size")

        red_df = df_atoms.drop_duplicates(subset="type")[["type", "number", "charge"]]

        red_df["sort_key"] = red_df["type"].apply(
            lambda x: (not str(x).isdigit(), int(x) if str(x).isdigit() else x)
        )
        red_df = red_df.sort_values("sort_key").drop(columns=["sort_key"])

        red_df["mass"] = red_df["type"].apply(
            lambda t: float(self._masses.get(t, MASSES_DICT.get(t, 1.0)))
        )

        return red_df

    def summary(self) -> None:
        """
        Print a complete summary of the AtomicSystem.

        """
        lines = []

        # Header
        output = f"<AtomicSystem with {self.num_atoms} atoms"

        sections = {
            "bonds": self.num_bonds,
            "angles": self.num_angles,
            "dihedrals": self.num_dihedrals,
            "impropers": self.num_impropers,
        }

        active_sections = [
            f"{count} {name}" for name, count in sections.items() if count > 0
        ]

        if active_sections:
            output += ", " + ", ".join(active_sections)

        lines.extend([output + ">", ""])

        # Box
        lines.append("Box")

        df_box = pd.DataFrame(
            np.reshape(self.box.T, (1, 6)),
            columns=["a (Å)", "b (Å)", "c (Å)", "α (°)", "β (°)", "γ (°)"],
        )

        lines.extend(
            [
                df_box.to_string(
                    index=False,
                    float_format="%.2f",
                ),
                "",
            ]
        )

        # Atoms
        lines.append("Atoms")

        df_atoms = self._get_type_summary()

        column_order = ["type", "number", "mass", "charge"]

        lines.extend(
            [
                df_atoms[column_order].to_string(index=False),
                "",
            ]
        )

        # Interactions
        def append_interaction_info(
            lines: list[str],
            df: pd.DataFrame | None,
            name: str,
        ) -> None:
            if df is None or df.empty:
                return

            lines.append(name)

            summary = df.groupby("type").size().rename("number").reset_index()

            lines.extend(
                [
                    summary.to_string(index=False),
                    "",
                ]
            )

        append_interaction_info(lines, self.bonds, "Bonds")
        append_interaction_info(lines, self.angles, "Angles")
        append_interaction_info(lines, self.dihedrals, "Dihedrals")
        append_interaction_info(lines, self.impropers, "Impropers")

        # Footer
        lines.extend(
            [
                f"Total charge: {self.total_charge:.3f} e",
                f"Volume: {self.volume / 1e3:.2f} nm³",
                f"Density: {self.density:.2f} g/cm³",
            ]
        )

        print("\n".join(lines))

    def add_structure(
        self,
        structure_to_add: AtomicSystem,
        distance: float = 2.0,
        axis: str = "z",
        vacuum: float = 10.0,
    ) -> Self:
        """
        Add a structure on top of this system.

        This method adds any structure (droplet, liquid layer, molecule, etc.)
        on top of the current system along the specified axis. The structure is
        aligned by center of mass in the transverse directions and placed at the
        specified distance from the surface.

        Parameters
        ----------
        structure_to_add : AtomicSystem or str
            The structure to add. Can be an AtomicSystem object or a path to a
            file that can be loaded by AtomicSystem.from_file().
        distance : float, default=2.0
            Distance between the current system surface and the structure
            in Ångströms.
        axis : str, default='z'
            Axis along which to add the structure ('x', 'y', or 'z').
        vacuum : float, default=10.0
            Vacuum space added above the structure in Ångströms.

        Returns
        -------
        Self
            The updated system with the structure added.

        Examples
        --------
        >>> from cemd import AtomicSystem
        >>>
        >>> # Add a droplet from a file
        >>> surface = AtomicSystem.from_file("surface.lmp")
        >>> droplet = AtomicSystem.from_file("droplet.lmp")
        >>> system = surface.add_structure(droplet, distance=2.0, vacuum=10.0)
        >>>
        >>> # Add a structure from a file path (auto-loads)
        >>> system = surface.add_structure("droplet.lmp", distance=2.0)
        >>>
        >>> # Add a custom structure with different axis
        >>> system = surface.add_structure(molecule, axis='y', distance=3.0)
        """
        from ..build import _add_structure

        new_system = _add_structure(
            solid_system=self,
            structure_to_add=structure_to_add,
            distance=distance,
            axis=axis,
            vacuum=vacuum,
        )
        self._replace_internals(new_system)

        return self

    def add_liquid_layer(
        self,
        blueprint: SolutionBuilder,
        thickness: float,
        distance: float = 2.0,
        vacuum: float = 10.0,
        axis: str = "z",
    ) -> Self:
        """
        Add a liquid layer on top of this system.

        This method creates a liquid layer with the composition defined by the
        blueprint and places it on top of the current system along the specified
        axis. The layer is aligned to match the transverse dimensions of the
        current system.

        Parameters
        ----------
        blueprint : SolutionBuilder
            Solution blueprint defining the liquid composition (density,
            solutes, etc.).
        thickness : float
            Thickness of the liquid layer in Ångströms.
        distance : float, default=2.0
            Distance between the current system surface and the liquid layer
            in Ångströms.
        vacuum : float, default=10.0
            Vacuum space added above the liquid layer in Ångströms.
        axis : str, default='z'
            Axis along which to add the layer ('x', 'y', or 'z').

        Returns
        -------
        Self
            The updated system with the liquid layer added.

        Examples
        --------
        >>> from cemd import AtomicSystem
        >>> from cemd.build import SolutionBuilder
        >>>
        >>> surface = AtomicSystem.from_file("surface.lmp")
        >>> blueprint = SolutionBuilder(
        ...     density=1.0,
        ...     molarities={'NaCl': 0.1}
        ... )
        >>>
        >>> # Add a 30 Å water layer on the surface
        >>> system = surface.add_liquid_layer(blueprint, thickness=30.0, distance=2.0)
        """
        from ..build import _add_liquid_layer

        new_system = _add_liquid_layer(
            self, blueprint, thickness, distance, vacuum, axis
        )
        self._replace_internals(new_system)

        return self

    def add_droplet(
        self,
        blueprint: SolutionBuilder,
        radius: float,
        distance: float = 2.0,
        vacuum: float = 10.0,
        axis: str = "z",
    ) -> Self:
        """
        Add a hemispherical liquid droplet on top of this system.

        This method creates a hemispherical droplet with the composition defined
        by the blueprint and places it on top of the current system along the
        specified axis. The droplet sits on the surface with a flat bottom.

        Parameters
        ----------
        blueprint : SolutionBuilder
            Solution blueprint defining the droplet composition (density,
            solutes, etc.).
        radius : float
            Radius of the hemispherical droplet in Ångströms.
        distance : float, default=2.0
            Distance between the current system surface and the droplet
            in Ångströms.
        vacuum : float, default=10.0
            Vacuum space added above the droplet in Ångströms.
        axis : str, default='z'
            Axis along which the droplet sits ('x', 'y', or 'z').

        Returns
        -------
        Self
            The updated system with the droplet added.

        Examples
        --------
        >>> from cemd import AtomicSystem
        >>> from cemd.build import SolutionBuilder
        >>>
        >>> surface = AtomicSystem.from_file("surface.lmp")
        >>> blueprint = SolutionBuilder(
        ...     density=1.0,
        ...     molarities={'NaCl': 0.1}
        ... )
        >>>
        >>> # Add a 15 Å radius water droplet on the surface
        >>> system = surface.add_droplet(blueprint, radius=15.0, distance=2.0)
        """
        from ..build import _add_droplet

        new_system = _add_droplet(self, blueprint, radius, distance, vacuum, axis)
        self._replace_internals(new_system)

        return self

    def get_count(self, symbol: str | int) -> int:
        """
        Calculate the exact number of atoms for a specific type.

        Parameters
        ----------
        symbol : str or int
            The atom type identifier to count.

        Returns
        -------
        int
            The number of atoms of the specified type.

        Raises
        ------
        ValueError
            If the symbol is not found.
        """
        summary = self._get_type_summary()

        mask = summary["type"].astype(str) == symbol

        if not mask.any():
            res = 0

        else:
            res = int(summary.loc[mask, "number"].sum())

        return res

    def get_center_of_mass(self) -> np.ndarray:
        """
        Calculates the center of mass (COM) of the system.

        Returns:
            np.ndarray: [x, y, z] coordinates of the center of mass.
        """
        atom_masses = self.atoms["type"].map(
            lambda t: self._masses.get(t, MASSES_DICT.get(t, 1.0))
        )

        total_mass = atom_masses.sum()

        if total_mass <= 0:
            return np.zeros(3)

        weighted_pos = self.atoms[["x", "y", "z"]].multiply(atom_masses, axis=0)

        return weighted_pos.sum().values / total_mass

    def view(self, trajectory: str = None) -> None:
        """
        Visualize the current system in VMD, with an optional trajectory.

        Parameters
        ----------
        trajectory : str, optional
            Path to a trajectory file to overlay onto this
            system's topology. All formats supported by MDAnalysis can be used. Defaults to None.

        Examples
        --------
        >>> system = AtomicSystem("input.data")
        >>> system.view()
        >>> system.view(trajectory="production.dcd")
        """

        view(self, trajectory=trajectory)
