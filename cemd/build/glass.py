# cemd/build/glass.py

from __future__ import annotations

import tempfile
import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from .._constants import AVOGADRO, CHARGES_DICT, MASSES_DICT
from ..core.atomic_system import AtomicSystem
from ._packmol import PackmolInput, PackmolStructure, get_structure_path, run_packmol
from .base import require_program

if TYPE_CHECKING:
    pass


@dataclass
class GlassBuilder:
    """
    Blueprint for creating a glass structure.

    This blueprint defines WHAT the glass contains (composition),
    not HOW it is built (geometry). The geometry is determined by
    the builder or method that uses this blueprint.

    Parameters
    ----------
    density : float
        Target density in g/cm³.
    composition : dict[str, float]
        Dictionary mapping formula strings to stoichiometric coefficients.
        Can be elements (e.g., 'Si', 'Al'), oxides (e.g., 'SiO2', 'Al2O3'),
        or molecules (e.g., 'H2O').
    structures : dict, optional
        Custom AtomicSystem structures for complex molecules.

    Examples
    --------
    >>> # Define a glass by elements
    >>> blueprint = GlassBuilder(
    ...     density=2.3,
    ...     composition={'Si': 3, 'Al': 4, 'Na': 4, 'O': 14}
    ... )
    >>> glass = blueprint.build(box=[20, 20, 20])
    >>>
    >>> # Define a glass by oxides
    >>> blueprint = GlassBuilder(
    ...     density=2.3,
    ...     composition={'SiO2': 3, 'Al2O3': 2, 'Na2O': 2}
    ... )
    >>> glass = blueprint.build(box=[20, 20, 20])
    >>>
    >>> # With water
    >>> blueprint = GlassBuilder(
    ...     density=2.3,
    ...     composition={'SiO2': 3, 'Al2O3': 2, 'Na2O': 2, 'H2O': 6}
    ... )
    """

    density: float = field(default=None)
    composition: dict[str, float] = field(default_factory=dict)
    structures: dict[str, AtomicSystem] = field(default_factory=dict)

    def __post_init__(self):
        """Validate the blueprint."""
        self._validate()

    def __repr__(self) -> str:
        """
        Return a styled string representation of the GlassBuilder.

        Returns
        -------
        str
            Summary of the glass composition and density.
        """
        lines = ["<GlassBuilder>"]
        lines.append("")

        # ====== Composition ======
        lines.append("┌─ Composition")
        lines.append(f"│   density:  {self.density:.2f} g/cm³")

        for formula, coeff in self.composition.items():
            lines.append(f"│             {formula}: {coeff:.3f}")

        # ====== Structures personnalisées ======
        if self.structures:
            lines.append("│")
            lines.append("├─ Custom structures")
            for species, struct in self.structures.items():
                lines.append(f"│   {species}: {struct.num_atoms} atoms")

        return "\n".join(lines)

    def _validate(self):
        """Check that the blueprint is valid."""
        if self.density is None or self.density <= 0:
            raise ValueError(f"density must be positive, got {self.density}")

        if not self.composition:
            raise ValueError("composition must not be empty")

        self._validate_composition()

        for species, structure in self.structures.items():
            if not isinstance(structure, AtomicSystem):
                raise TypeError(
                    f"Structure for '{species}' must be an AtomicSystem, "
                    f"got {type(structure).__name__}"
                )

        # Check charge neutrality.
        charge = self.get_total_charge()

        if not np.isclose(charge, 0.0):
            warnings.warn(
                f"Composition is not charge neutral: "
                f"total formal charge = {charge:+.3f}",
                UserWarning,
                stacklevel=2,
            )

    def _validate_composition(self):
        """Validate the chemical composition."""
        from pymatgen.core import Composition

        for formula, coefficient in self.composition.items():
            if not isinstance(formula, str) or not formula:
                raise TypeError(
                    f"Composition keys must be non-empty strings, got {formula!r}"
                )

            if not isinstance(coefficient, (int, float, np.integer, np.floating)):
                raise TypeError(
                    f"Coefficient for '{formula}' must be a number, "
                    f"got {type(coefficient).__name__}"
                )

            if coefficient <= 0:
                raise ValueError(
                    f"Coefficient for '{formula}' must be positive, got {coefficient}"
                )

            # Custom structure: no need to parse the formula.
            if formula in self.structures:
                continue

            try:
                composition = Composition(formula)
            except Exception as exc:
                raise ValueError(
                    f"Invalid chemical formula '{formula}'. "
                    "Expected an element or a valid chemical formula "
                    "(e.g. 'Si', 'SiO2', 'Al2O3', 'Na2O')."
                ) from exc

            if not composition:
                raise ValueError(
                    f"Chemical formula '{formula}' does not contain any elements."
                )

    def get_mass_per_formula_unit(self) -> float:
        """
        Calculate mass per formula unit.

        A formula unit is defined by the stoichiometric coefficients
        in the composition dictionary.

        Returns
        -------
        float
            Mass per formula unit in amu.

        Raises
        ------
        ValueError
            If a species is not found in the mass database and no custom
            structure is provided.
        """
        total_mass = 0.0

        for formula, coeff in self.composition.items():
            if formula in self.structures:
                mass = self.structures[formula].total_mass
            else:
                # Essayer de décomposer la formule avec pymatgen
                try:
                    from pymatgen.core import Composition

                    comp = Composition(formula)
                    mass = sum(
                        MASSES_DICT.get(el, 0) * amt
                        for el, amt in comp.get_el_amt_dict().items()
                    )
                except Exception:
                    # Si pymatgen échoue, essayer comme élément simple
                    if formula in MASSES_DICT:
                        mass = MASSES_DICT[formula]
                    else:
                        raise ValueError(
                            f"Formula '{formula}' not found in mass database "
                            f"and no custom structure provided.\n"
                            f"Available: {list(MASSES_DICT.keys())}"
                        )

            total_mass += mass * coeff

        return total_mass

    def get_elemental_composition(self) -> dict[str, float]:
        """Get the elemental composition from the blueprint."""
        from pymatgen.core import Composition

        elemental: dict[str, float] = {}

        for formula, coeff in self.composition.items():
            if formula in self.structures:
                structure = self.structures[formula]
                elem_map = structure.elements

                # Count per atom (not per unique atom type), so a custom
                # structure like H2O contributes 2 H and 1 O per coeff,
                # matching the pymatgen branch below.
                for atype in structure.atoms["type"]:
                    symbol = elem_map.get(atype, str(atype))
                    elemental[symbol] = elemental.get(symbol, 0.0) + coeff

                continue

            composition = Composition(formula)

            for symbol, amount in composition.get_el_amt_dict().items():
                elemental[symbol] = elemental.get(symbol, 0.0) + amount * coeff

        return elemental

    def get_num_formula_units(self, volume: float) -> int:
        """
        Calculate number of formula units for a given volume.

        Parameters
        ----------
        volume : float
            Volume in Å³

        Returns
        -------
        int
            Number of formula units.
        """
        mass_per_unit = self.get_mass_per_formula_unit()
        num_units = int(
            np.round((self.density * volume * 1e-24 * AVOGADRO) / mass_per_unit)
        )
        return max(1, num_units)

    def get_actual_density(self, volume: float) -> float:
        """
        Calculate actual density for a given volume.

        Parameters
        ----------
        volume : float
            Volume in Å³

        Returns
        -------
        float
            Actual density in g/cm³.
        """
        num_units = self.get_num_formula_units(volume)
        mass_per_unit = self.get_mass_per_formula_unit()
        return (num_units * mass_per_unit) / (AVOGADRO * volume * 1e-24)

    def get_total_charge(self) -> float:
        """
        Calculate the total formal charge of the composition.

        Returns
        -------
        float
            Total formal charge.
        """
        elemental = self.get_elemental_composition()

        return sum(
            amount * CHARGES_DICT[element]
            for element, amount in elemental.items()
            if element in CHARGES_DICT
        )

    def build(
        self,
        box: Sequence[float],
        margin: float = 0.95,
    ) -> AtomicSystem:
        """
        Build a standalone glass system from this blueprint.

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
            The built glass system.

        Examples
        --------
        >>> blueprint = GlassBuilder(
        ...     density=2.3,
        ...     composition={"SiO2": 3, "Al2O3": 2, "Na2O": 2}
        ... )
        >>> glass = blueprint.build(box=[20, 20, 20])
        """

        require_program("packmol")

        box = list(box)
        boxa, boxb, boxc = box[:3]
        volume = boxa * boxb * boxc

        num_units = self.get_num_formula_units(volume)

        actual_density = self.get_actual_density(volume)
        density_error = (actual_density - self.density) / self.density

        if abs(density_error) > 0.05:
            warnings.warn(
                f"Density error: {density_error:.1%}. "
                f"Target: {self.density:.3f} g/cm³, "
                f"Actual: {actual_density:.3f} g/cm³"
            )

        elemental_comp = self.get_elemental_composition()

        inside_box = (
            0.0,
            0.0,
            0.0,
            boxa * margin,
            boxb * margin,
            boxc * margin,
        )

        with tempfile.TemporaryDirectory(dir=".") as tmp:
            structures = []

            for element, amount in elemental_comp.items():
                count = int(np.round(amount * num_units))

                if count <= 0:
                    continue

                structure = self.structures.get(element, element)

                if isinstance(structure, str):
                    structure = get_structure_path(structure, tmp)

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
                output="glass.pdb",
                filetype="pdb",
                structures=structures,
            )

            data = run_packmol(packmol_input)

        final_box = box[:3] + [90.0, 90.0, 90.0] if len(box) == 3 else box
        data.set_box(final_box)

        data._glass_metadata = {
            "density": self.density,
            "actual_density": actual_density,
            "composition": self.composition,
            "num_units": num_units,
            "volume": volume,
            "blueprint": self,
        }

        return data

    @classmethod
    def from_stoichiometry(
        cls,
        density: float,
        composition: dict[str, float],
        structures: dict[str, AtomicSystem] | None = None,
    ) -> GlassBuilder:
        """
        Create a blueprint from stoichiometric coefficients.

        Parameters
        ----------
        density : float
            Target density in g/cm³.
        composition : dict
            Dictionary mapping formula strings to stoichiometric coefficients.
        structures : dict, optional
            Custom structures for species.

        Returns
        -------
        GlassBuilder
            The blueprint.

        Examples
        --------
        >>> # By elements
        >>> blueprint = GlassBuilder.from_stoichiometry(
        ...     density=2.3,
        ...     composition={'Si': 3, 'Al': 4, 'Na': 4, 'O': 14}
        ... )
        >>>
        >>> # By oxides
        >>> blueprint = GlassBuilder.from_stoichiometry(
        ...     density=2.3,
        ...     composition={'SiO2': 3, 'Al2O3': 2, 'Na2O': 2}
        ... )
        """
        return cls(
            density=density, composition=composition, structures=structures or {}
        )

    @classmethod
    def from_oxides(
        cls,
        density: float,
        oxides: dict[str, float],
        structures: dict[str, AtomicSystem] | None = None,
    ) -> GlassBuilder:
        """
        Create a blueprint from oxide composition (alias for from_stoichiometry).

        Parameters
        ----------
        density : float
            Target density in g/cm³.
        oxides : dict
            Dictionary mapping oxide formulas to stoichiometric coefficients.
        structures : dict, optional
            Custom structures for species.

        Returns
        -------
        GlassBuilder
            The blueprint.
        """
        return cls.from_stoichiometry(density, oxides, structures)

    @classmethod
    def from_elements(
        cls,
        density: float,
        elements: dict[str, float],
        structures: dict[str, AtomicSystem] | None = None,
    ) -> GlassBuilder:
        """
        Create a blueprint from elemental composition.

        Parameters
        ----------
        density : float
            Target density in g/cm³.
        elements : dict
            Dictionary mapping element symbols to stoichiometric coefficients.
        structures : dict, optional
            Custom structures for species.

        Returns
        -------
        GlassBuilder
            The blueprint.

        Examples
        --------
        >>> blueprint = GlassBuilder.from_elements(
        ...     density=2.3,
        ...     elements={'Si': 3, 'Al': 4, 'Na': 4, 'O': 14}
        ... )
        """
        return cls.from_stoichiometry(density, elements, structures)

    def is_pure(self) -> bool:
        """Check if the blueprint contains only one component."""
        return len(self.composition) == 1

    def get_components(self) -> list[str]:
        """Get list of components in the composition."""
        return list(self.composition.keys())

    def get_coefficient(self, formula: str) -> float:
        """Get the stoichiometric coefficient for a specific formula."""
        return self.composition.get(formula, 0.0)

    def normalize(self) -> GlassBuilder:
        """
        Normalize the composition so that the sum of coefficients equals 1.

        Returns
        -------
        GlassBuilder
            A new builder with normalized composition.
        """
        total = sum(self.composition.values())
        if total == 0:
            raise ValueError("Cannot normalize empty composition")

        normalized = {k: v / total for k, v in self.composition.items()}
        return GlassBuilder(
            density=self.density, composition=normalized, structures=self.structures
        )
