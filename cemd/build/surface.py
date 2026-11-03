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

import logging
from collections.abc import Sequence
from typing import TYPE_CHECKING

from ..core._io.formats.pmg import order_by_majority
from ..core.atomic_system import AtomicSystem
from .base import BaseBuilder

if TYPE_CHECKING:
    from pymatgen.core import Structure

logger = logging.getLogger(__name__)


class SurfaceBuilder(BaseBuilder):
    """Builder for generating surfaces."""

    def __init__(self, system: AtomicSystem):
        super().__init__(system)
        self._structure = self._get_structure()
        self.result = None

    def __repr__(self) -> str:
        """
        Return a styled string representation of the SurfaceBuilder.

        Returns
        -------
        str
            Summary of the builder state with box and atom information.
        """

        a, b, c = self._structure.lattice.abc
        alpha, beta, gamma = self._structure.lattice.angles

        lines = ["<SurfaceBuilder>"]
        lines.append("")

        # ====== System ======
        lines.append("┌─ Structure")
        lines.append(f"│   sites   : {len(self._structure)}")
        lines.append(f"│   formula : {self._structure.formula}")
        lines.append(f"│   lattice : a={a:.2f}  b={b:.2f}  c={c:.2f} Å")
        lines.append(f"│             α={alpha:.1f}  β={beta:.1f}  γ={gamma:.1f} °")
        lines.append(f"│   ordered : {'✓' if self._structure.is_ordered else '✗'}")

        return "\n".join(lines)

    def _get_structure(self) -> Structure:
        """Get the structure for surface generation."""
        from pymatgen.core.surface import Slab

        # A system read from a format that carries no crystal (.data, .xyz)
        # is taken as it is: nothing tells whether it is still a bulk.
        if self.system._pmg_struct is None:
            return self.system.to_pmg()

        if isinstance(self.system._pmg_struct, Slab):
            raise ValueError(
                "A slab can only be cut from a bulk crystal, and this system "
                "is already a slab. Cut the new orientation from the bulk "
                "structure."
            )

        # The structure the system was read from (Miller indices refer to its
        # axes), as long as the system still matches it. A copy, since `build`
        # adds oxidation states to it, ordered the way the reader ordered the
        # atoms: SlabGenerator refuses a disordered structure ("species
        # property only supports ordered structures"), e.g. the Ca0.99Mg0.01
        # sites of a magnesian calcite (COD 9017359).
        cached = self.system._cached_pmg()
        if cached is None:
            raise ValueError(
                "A slab can only be cut from a bulk crystal, and this system "
                "has been edited since it was read (atoms moved, added or "
                "removed). Cut the slab from the structure as read, then edit "
                "the slab."
            )
        return order_by_majority(cached)

    def build(
        self,
        miller_indices: Sequence[int],
        min_slab_size: float = 25.0,
        min_vacuum_size: float = 15.0,
        max_broken_bonds: int = 20,
        bonds_dict: dict | None = None,
    ) -> tuple[list[AtomicSystem], list[float], list[float], int]:
        """
        Generate crystalline surfaces from the structure.

        Parameters
        ----------
        miller_indices : Sequence[int]
            Miller indices for the surface cut (e.g., (1, 1, 1)).
        min_slab_size : float, default=25.0
            Minimum thickness of the slab in Angstroms.
        min_vacuum_size : float, default=15.0
            Minimum size of the added vacuum in Angstroms.
        max_broken_bonds : int, default=20
            Maximum number of broken bonds to accept.
        bonds_dict : dict, optional
            Dictionary of bond lengths for different element pairs.
            If None, uses default values.

        Returns
        -------
        tuple
            (slabs_list, shift_list, dipole_list, actual_broken)

        Raises
        ------
        RuntimeError
            If no structure is available for surface generation.
        ValueError
            If the structure is not suitable for surface generation.
        """

        from pymatgen.core.surface import SlabGenerator

        struct = self._structure

        if len(struct) < 2:
            raise ValueError(
                f"Structure has only {len(struct)} atom(s). "
                "Surface generation requires at least 2 atoms."
            )

        if not struct.is_ordered:
            logger.warning(
                "Structure is not fully ordered. Surface generation may fail."
            )

        # Check that the structure is not too large
        if len(struct) > 10000:
            logger.warning(
                f"Structure has {len(struct)} atoms. "
                "Surface generation may be slow or fail."
            )

        try:
            struct.add_oxidation_state_by_guess()
            logger.info("Oxidation states guessed successfully.")
        except Exception as e:
            logger.warning(f"Could not guess oxidation states: {e}")

        if bonds_dict is None:
            bonds_dict = {
                # With oxidation states
                ("Al4+", "O2-"): 1.8,
                ("Si4+", "O2-"): 1.8,
                ("C4+", "O2-"): 1.6,
                ("H+", "O2-"): 1.1,
                # Without oxidation states (fallback)
                ("Al", "O"): 1.8,
                ("Si", "O"): 1.8,
                ("C", "O"): 1.6,
                ("H", "O"): 1.1,
            }

        try:
            slabgen = SlabGenerator(
                struct,
                miller_indices,
                min_slab_size,
                min_vacuum_size,
                primitive=True,
                lll_reduce=True,
            )
            logger.info(f"SlabGenerator created for {miller_indices}")
        except Exception as e:
            raise ValueError(
                f"Cannot create SlabGenerator: {e}\n"
                f"Miller indices: {miller_indices}\n"
                "Make sure the structure is periodic and valid."
            ) from e

        slabs_list_pmg = []
        actual_broken = -1
        for i in range(max_broken_bonds + 1):
            try:
                slabs = slabgen.get_slabs(bonds=bonds_dict, max_broken_bonds=i)
                if len(slabs) > 0:
                    slabs_list_pmg = slabs
                    actual_broken = i
                    logger.info(f"Found {len(slabs)} slabs with {i} broken bonds.")
                    break
            except Exception as e:
                logger.debug(f"Failed with {i} broken bonds: {e}")
                continue

        # Fallback: try without bonds
        if len(slabs_list_pmg) == 0:
            try:
                slabs_list_pmg = slabgen.get_slabs()
                if len(slabs_list_pmg) > 0:
                    actual_broken = -1
                    logger.warning(
                        f"Found {len(slabs_list_pmg)} slabs using fallback "
                        "(no bond constraints)."
                    )
                else:
                    raise RuntimeError(
                        f"No surfaces could be generated with up to {max_broken_bonds} broken bonds. "
                        "The structure may not be suitable for surface generation."
                    )
            except Exception as e:
                raise RuntimeError(
                    f"Surface generation failed: {e}\n"
                    "Try adjusting min_slab_size, min_vacuum_size, or max_broken_bonds."
                ) from e

        slabs_list_as = []
        dipole_list = []
        shift_list = []
        for i, raw_slab in enumerate(slabs_list_pmg):
            slab = raw_slab.get_orthogonal_c_slab().get_sorted_structure()
            shift_list.append(slab.shift)
            dipole_list.append(0.0)
            as_system = AtomicSystem.from_pmg(slab, refine=False).wrap()
            as_system._cache_pmg(slab)
            slabs_list_as.append(as_system)

        if len(slabs_list_as) == 0:
            raise RuntimeError("No slabs could be converted to AtomicSystem.")

        logger.info(f"Successfully generated {len(slabs_list_as)} surfaces.")

        self.result = (slabs_list_as, shift_list, dipole_list, actual_broken)

        return slabs_list_as, shift_list, dipole_list, actual_broken

    def explore(self) -> AtomicSystem | None:
        """
        Interactive surface explorer.

        Allows the user to interactively select Miller indices,
        slab thickness, and vacuum size, then generates and displays
        the available surfaces.

        Returns
        -------
        AtomicSystem or None
            The selected surface, or None if the user cancels.
        """
        import questionary
        from prompt_toolkit.application import Application
        from prompt_toolkit.key_binding import KeyBindings
        from prompt_toolkit.layout import Layout
        from prompt_toolkit.layout.containers import HSplit, Window
        from prompt_toolkit.layout.controls import FormattedTextControl

        miller_str = questionary.text(
            "Enter Miller indices (h k l):", default="1 0 0"
        ).ask()
        if not miller_str:
            return None

        try:
            miller = [int(x) for x in miller_str.split()]
            if len(miller) != 3:
                raise ValueError("Need exactly 3 indices")
        except ValueError as e:
            print(f"Invalid Miller indices: {e}")
            return None

        min_slab = questionary.text("Min slab size (Å):", default="25.0").ask()
        if not min_slab:
            return None

        min_vac = questionary.text("Min vacuum size (Å):", default="15.0").ask()
        if not min_vac:
            return None

        try:
            slabs, shifts, dipoles, n_broken = self.build(
                miller_indices=miller,
                min_slab_size=float(min_slab),
                min_vacuum_size=float(min_vac),
            )
        except Exception as e:
            print(f"Error generating surfaces: {e}")
            return None

        if not slabs:
            print("\nNo surfaces generated.")
            return None

        # Format display
        header_format = "{cursor} {idx:>4} | {shift:>10} | {dipole:>12} | {natoms:>8}"
        row_format = (
            "{cursor} {idx:>4} | {shift:>10.4f} | {dipole:>12.4f} | {natoms:>8}"
        )
        index = 0

        def render():
            header = header_format.format(
                cursor=" ",
                idx="#",
                shift="Shift",
                dipole="Dipole (D)",
                natoms="N atoms",
            )

            if n_broken == -1:
                generated_str = f"\n✓ Generated {len(slabs)} surfaces"
            else:
                generated_str = (
                    f"\n✓ Generated {len(slabs)} surfaces with {n_broken} broken bonds"
                )

            lines = [
                generated_str,
                "↑↓ Move   [Enter] Select   [q] Quit",
                "",
                header,
                "-" * len(header),
            ]
            for i, (slab, shift, dipole) in enumerate(zip(slabs, shifts, dipoles)):
                cursor = "➜" if i == index else " "
                lines.append(
                    row_format.format(
                        cursor=cursor,
                        idx=i + 1,
                        shift=shift,
                        dipole=dipole,
                        natoms=slab.num_atoms,
                    )
                )
            return "\n".join(lines)

        selected_idx = [0]
        kb = KeyBindings()

        @kb.add("up")
        def _(e):
            nonlocal index
            if index > 0:
                index -= 1
            e.app.invalidate()

        @kb.add("down")
        def _(e):
            nonlocal index
            if index < len(slabs) - 1:
                index += 1
            e.app.invalidate()

        @kb.add("enter")
        def _(e):
            selected_idx[0] = index
            e.app.exit()

        @kb.add("q")
        @kb.add("escape")
        def _(e):
            selected_idx[0] = -1
            e.app.exit()

        window = Window(
            content=FormattedTextControl(render),
            always_hide_cursor=True,
        )
        app = Application(
            layout=Layout(HSplit([window])), key_bindings=kb, full_screen=False
        )
        app.run()

        if selected_idx[0] == -1:
            return None

        return slabs[selected_idx[0]]
