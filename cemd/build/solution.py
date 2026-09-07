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

import tempfile
import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field
from functools import cache
from typing import TYPE_CHECKING

import numpy as np

from .._constants import AVOGADRO, MASSES_DICT
from ._packmol import (
    PackmolInput,
    PackmolStructure,
    get_structure_path,
    run_packmol,
)
from .base import require_program

if TYPE_CHECKING:
    from ..core.atomic_system import AtomicSystem

WATER_MOLAR_MASS: float = 2 * MASSES_DICT["H"] + MASSES_DICT["O"]

# Extensions probed by `_packmol._resolve_structure`, in its own order.
_BUNDLED_EXTENSIONS: tuple[str, ...] = (".lt", ".pdb", ".sdf")


def _bundled_structure_names() -> set[str]:
    """Names of the polyatomic species shipped with CEMD."""
    from ._structures import STRUCTURES_DIR

    return {
        path.stem.upper()
        for extension in _BUNDLED_EXTENSIONS
        for path in STRUCTURES_DIR.glob(f"*{extension}")
    }


@cache
def _bundled_structure_mass(species: str) -> float | None:
    """Total mass of a species shipped with CEMD, or None if it has none."""
    from ..core.atomic_system import AtomicSystem
    from ._structures import STRUCTURES_DIR

    for extension in _BUNDLED_EXTENSIONS:
        path = STRUCTURES_DIR / f"{species.lower()}{extension}"
        if path.exists():
            return float(AtomicSystem.from_file(path).total_mass)

    return None


def _concentration2count(molarity: float | int, volume: float) -> tuple[int, float]:
    """Calculates the integer particle count and relative error for a single molarity.

    Parameters
    ----------
    molarity : float | int
        Target concentration in mol/L (M).
    volume : float
        Volume of the simulation box in cubic Angstroms (Å³).

    Returns
    -------
    tuple[int, float]
        -Particle count (integer).
        -Relative discretization error percentage (float).
    """
    if molarity == 0:
        return 0, 0.0

    theoretical_n = molarity * AVOGADRO * volume * 1e-27
    final_count = int(round(theoretical_n))

    if theoretical_n > 0:
        error_pct = abs(final_count - theoretical_n) / theoretical_n
    else:
        error_pct = 0.0

    return final_count, error_pct


@dataclass
class SolutionBuilder:
    """
    Blueprint for creating a solution.

    This blueprint defines WHAT the solution contains (composition),
    not HOW it is built (geometry). The geometry is determined by
    the builder or method that uses this blueprint.

    Parameters
    ----------
    density : float, default=1.0
        Target density in g/cm³
    molarities : dict, optional
        Species molarities in mol/L
    counts : dict, optional
        Species explicit counts (number of molecules/ions)
    structures : dict, optional
        Custom AtomicSystem structures for complex molecules

    Examples
    --------
    >>> # Define a salt solution
    >>> blueprint = SolutionBuilder(
    ...     density=1.0,
    ...     molarities={'NaCl': 0.1, 'KCl': 0.05}
    ... )

    >>> # Use in different contexts
    >>> solution = blueprint.build(box=[30, 30, 30])  # Standalone
    >>> system = surface.add_liquid_layer(blueprint, thickness=30.0)  # Layer
    >>> system = surface.add_droplet(blueprint, radius=15.0)  # Droplet

    >>> # Pure water
    >>> water = SolutionBuilder.from_water()

    >>> # With explicit counts
    >>> blueprint = SolutionBuilder(
    ...     density=1.0,
    ...     counts={'Na': 50, 'Cl': 50}
    ... )
    """

    density: float = 1.0
    molarities: dict[str, float] = field(default_factory=dict)
    counts: dict[str, int] = field(default_factory=dict)
    structures: dict[str, AtomicSystem] = field(default_factory=dict)

    def __post_init__(self):
        """Validate the blueprint."""
        self._validate()

    def __repr__(self) -> str:
        """
        Return a styled string representation of the SolutionBuilder.

        Returns
        -------
        str
            Summary of the solution composition and density.
        """
        lines = ["<SolutionBuilder>"]
        lines.append("")

        # ====== Composition ======
        lines.append("┌─ Composition")
        lines.append(f"│   density:  {self.density:.2f} g/cm³")

        # Solutes
        if self.molarities:
            lines.append("│   type:     molarities")
            for species, molarity in self.molarities.items():
                lines.append(f"│             {species}: {molarity:.3f} M")
        elif self.counts:
            lines.append("│   type:     counts")
            for species, count in self.counts.items():
                lines.append(f"│             {species}: {count}")
        else:
            lines.append("│   type:     pure water")

        # ====== Custom structures ======
        if self.structures:
            lines.append("│")
            lines.append("├─ Custom structures")
            for species, struct in self.structures.items():
                lines.append(f"│   {species}: {struct.num_atoms} atoms")

        return "\n".join(lines)

    def _validate(self):
        """Check that the blueprint is valid."""
        if self.density <= 0:
            raise ValueError(f"Density must be positive, got {self.density}")

        # Validate molarities
        for species, value in self.molarities.items():
            if not isinstance(value, (int, float)):
                raise TypeError(
                    f"Molarity for '{species}' must be a number, "
                    f"got {type(value).__name__}"
                )
            if value < 0:
                raise ValueError(f"Molarity must be >= 0 for '{species}'")

        # Validate counts
        for species, value in self.counts.items():
            if not isinstance(value, int):
                raise TypeError(
                    f"Count for '{species}' must be an integer, "
                    f"got {type(value).__name__}"
                )
            if value < 0:
                raise ValueError(f"Count must be >= 0 for '{species}'")

        # Validate custom structures
        for species, struct in self.structures.items():
            from ..core.atomic_system import AtomicSystem

            if not isinstance(struct, AtomicSystem):
                raise TypeError(
                    f"Structure for '{species}' must be an AtomicSystem, "
                    f"got {type(struct).__name__}"
                )

    @classmethod
    def from_water(
        cls,
        density: float = 1.0,
    ) -> SolutionBuilder:
        """
        Create a blueprint for pure water.

        Parameters
        ----------
        density : float, default=1.0
            Target density in g/cm³

        Returns
        -------
        SolutionBuilder
            A blueprint configured for pure water.

        Examples
        --------
        >>> water = SolutionBuilder.from_water()
        >>> solution = water.build(box=[30, 30, 30])
        """
        return cls(
            density=density,
            molarities={},
            counts={},
            structures={},
        )

    def to_counts(self, volume: float) -> dict[str, int]:
        """
        Convert molarities to explicit counts for a given volume and
        combine them with explicitly specified counts.
        """
        counts = self.counts.copy()

        for species, molarity in self.molarities.items():
            if species in counts:
                raise ValueError(
                    f"Species '{species}' is specified in both "
                    "`counts` and `molarities`."
                )

            count, error = _concentration2count(molarity, volume)
            counts[species] = count

            if error > 0.1:
                warnings.warn(
                    f"Species '{species}' has {error:.1%} discretization error. "
                    f"Using {count} molecules (target: {molarity} M)."
                )

        return counts

    def get_solute_mass(self, volume: float) -> float:
        """
        Calculate total mass of solutes in atomic mass units for a given volume.

        Parameters
        ----------
        volume : float
            Volume in Å³

        Returns
        -------
        float
            Total solute mass in amu

        Raises
        ------
        ValueError
            If a species is not found in the mass database.
        """
        counts = self.to_counts(volume)
        total_mass = 0.0

        for species, count in counts.items():
            if species in self.structures:
                total_mass += count * self.structures[species].total_mass
                continue

            # Same resolution order as the packer (`_resolve_structure`):
            # a bundled file wins over the element table, so that the mass
            # used here describes whatever actually gets packed. Without
            # this, the polyatomic species shipped with CEMD -- "HO",
            # "CO3", "SO4" -- could be packed by name but not weighed, and
            # had to be passed again through `structures=`.
            bundled = _bundled_structure_mass(species)
            if bundled is not None:
                total_mass += count * bundled
            elif species in MASSES_DICT:
                total_mass += count * MASSES_DICT[species]
            else:
                raise ValueError(
                    f"Species '{species}' not found in mass database "
                    f"and no custom structure provided.\n"
                    f"Available: {list(MASSES_DICT.keys())} "
                    f"and the bundled structures "
                    f"{sorted(_bundled_structure_names())}"
                )

        return total_mass

    def get_water_count(self, volume: float) -> int:
        """
        Calculate number of water molecules for a given volume.

        Parameters
        ----------
        volume : float
            Volume in Å³

        Returns
        -------
        int
            Number of water molecules

        Raises
        ------
        ValueError
            If the density or volume is too low for the solutes.
        """
        mass_total_target_g = self.density * volume * 1e-24
        mass_solutes_g = self.get_solute_mass(volume) / AVOGADRO
        mass_water_g = mass_total_target_g - mass_solutes_g

        if mass_water_g <= 0:
            raise ValueError(
                f"Target density ({self.density} g/cm³) or volume ({volume:.1f} Å³) "
                f"is too low to hold the solutes.\n"
                f"Solute mass: {mass_solutes_g:.4f} g\n"
                f"Target total mass: {mass_total_target_g:.4f} g"
            )

        return int(round((mass_water_g * AVOGADRO) / WATER_MOLAR_MASS))

    def build(
        self,
        box: Sequence[float],
        margin: float = 0.95,
    ) -> AtomicSystem:
        """
        Build a standalone solution system from this blueprint.

        Parameters
        ----------
        box : Sequence[float]
            Box dimensions [a, b, c] in Å.
        margin : float, default=0.95
            Safety scaling factor to prevent atoms from sitting directly
            on box edges.

        Returns
        -------
        AtomicSystem
            The built solution system.
        """
        from ..core.atomic_system import AtomicSystem

        require_program("packmol")

        box = list(box)
        boxa, boxb, boxc = box[:3]
        volume = boxa * boxb * boxc

        solute_counts = self.to_counts(volume)
        num_water = self.get_water_count(volume)

        inside_box = (
            0.0,
            0.0,
            0.0,
            boxa * margin,
            boxb * margin,
            boxc * margin,
        )

        with tempfile.TemporaryDirectory(dir=".") as tmp:
            h2o_path = get_structure_path("H2O", tmp)
            h2o = AtomicSystem.from_file(h2o_path)
            h2o.set_types({"H1": "Hw", "H2": "Hw", "O1": "Ow"})

            structures = [
                PackmolStructure(
                    structure=h2o,
                    number=num_water,
                    inside_box=inside_box,
                )
            ]

            for species, count in solute_counts.items():
                if count <= 0:
                    continue

                structure = self.structures.get(species, species)

                structures.append(
                    PackmolStructure(
                        structure=structure,
                        number=count,
                        center=True,
                        inside_box=inside_box,
                    )
                )

            packmol_input = PackmolInput(
                tolerance=2.0,
                output="solution.pdb",
                filetype="pdb",
                structures=structures,
            )

            data = run_packmol(packmol_input)

        final_box = box[:3] + [90.0, 90.0, 90.0] if len(box) == 3 else box
        data.set_box(final_box)

        data._solution_metadata = {
            "density": self.density,
            "solutes": solute_counts,
            "num_water": num_water,
            "volume": volume,
            "blueprint": self,
        }

        return data

    def build_hemisphere(
        self,
        radius: float,
        axis: str = "z",
    ) -> AtomicSystem:
        """
        Build a hemispherical solution system from this blueprint.

        The hemisphere is centered at the origin, with its flat surface
        lying in the plane perpendicular to ``axis``.

        Parameters
        ----------
        radius : float
            Hemisphere radius in Å.
        axis : {"x", "y", "z"}, default="z"
            Axis normal to the flat surface of the hemisphere.

        Returns
        -------
        AtomicSystem
            The built hemispherical solution system.

        Raises
        ------
        ValueError
            If ``radius`` is not positive or ``axis`` is invalid.
        """
        from ..core.atomic_system import AtomicSystem

        require_program("packmol")

        if radius <= 0:
            raise ValueError("The radius must be positive.")

        if axis not in {"x", "y", "z"}:
            raise ValueError(f"Invalid axis: {axis!r}. Expected 'x', 'y', or 'z'.")

        # Hemisphere volume
        volume = 2.0 * np.pi * radius**3 / 3.0

        # Number of solute molecules
        solute_counts = self.to_counts(volume)

        # Number of water molecules
        num_water = self.get_water_count(volume)

        # Packmol half-space definition
        above_plane = {
            "x": (1.0, 0.0, 0.0, 0.0),
            "y": (0.0, 1.0, 0.0, 0.0),
            "z": (0.0, 0.0, 1.0, 0.0),
        }[axis]
        inside_sphere = (0.0, 0.0, 0.0, radius)

        with tempfile.TemporaryDirectory(dir=".") as tmp:
            h2o_path = get_structure_path("H2O", tmp)
            h2o = AtomicSystem.from_file(h2o_path)
            # Match `build()`'s naming convention: downstream topology
            # rules (e.g. CLAYFF_RULES) and force-field assignment key off
            # "Ow"/"Hw", not the raw template types "O1"/"H1"/"H2".
            h2o.set_types({"H1": "Hw", "H2": "Hw", "O1": "Ow"})

            structures = [
                PackmolStructure(
                    structure=h2o,
                    number=num_water,
                    inside_sphere=inside_sphere,
                    above_plane=above_plane,
                )
            ]

            # Solutes
            for species, count in solute_counts.items():
                if count <= 0:
                    continue

                structure = self.structures.get(species, species)

                structures.append(
                    PackmolStructure(
                        structure=structure,
                        number=count,
                        center=True,
                        inside_sphere=inside_sphere,
                        above_plane=above_plane,
                    )
                )

            packmol_input = PackmolInput(
                tolerance=2.0,
                output="solution_hemisphere.pdb",
                filetype="pdb",
                structures=structures,
            )

            data = run_packmol(packmol_input)
        # Simulation box
        box_size = 2.5 * radius
        final_box = [
            box_size,
            box_size,
            1.5 * radius,
            90.0,
            90.0,
            90.0,
        ]

        data.set_box(final_box)

        data._solution_metadata = {
            "density": self.density,
            "solutes": solute_counts,
            "num_water": num_water,
            "volume": volume,
            "radius": radius,
            "axis": axis,
            "blueprint": self,
        }

        return data
