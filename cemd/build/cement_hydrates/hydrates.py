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
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ...analysis import analyze_silicates
from ...core.atomic_system import AtomicSystem
from .._structures import STRUCTURES_DIR
from ..base import BaseBuilder, require_program
from ._silicate_helpers import substitute_si_by_al

if TYPE_CHECKING:
    pass


# =============================================================================
# CSHBuilder
# =============================================================================


@dataclass
class CSHBuilder(BaseBuilder):
    """
    Builder for creating C-S-H (Calcium-Silicate-Hydrate) structures.

    C-S-H is built from a tobermorite model: bridging silicates are removed
    to reach the requested Ca/Si ratio, then the interlayers are filled with
    water and charge-balancing calcium.

    Parameters
    ----------
    cs_ratio : float
        Calcium/Silicate ratio.
    ws_ratio : float
        Water/Silicate ratio, counted *before* the interlayer calcium: the
        builder adds ``round(n_Si * ws_ratio)`` water molecules plus one
        more per Ca²⁺ inserted to reach a high ``cs_ratio``. Up to
        Ca/Si ~ 1.33 the target is met by removing bridging silicates
        alone, no calcium is inserted, and the H2O/Si reported by
        :meth:`analyze` comes back exactly equal to ``ws_ratio``. Past that
        threshold ``min_mcl`` (see :meth:`build`) forbids further
        vacancies, so calcium is inserted instead, one extra water molecule
        per Ca²⁺. ``neutralize_csh_charge`` (``_silicate_helpers.py``) then
        strips two hydrogens per inserted Ca²⁺ -- one each from up to two
        water molecules, picked at random -- to neutralize the +2 charge,
        which removes exactly as many hydrogens as that extra water added.
        So the measured H2O/Si still comes back to ``ws_ratio`` (checked at
        Ca/Si 1.4: 1.150).
    progress_callback : Callable, optional
        Progress callback function.

    Examples
    --------
    >>> builder = CSHBuilder(cs_ratio=1.5, ws_ratio=1.0)
    >>> system = builder.build(
    ...     supercell=[4,1,1],
    ...     model='tob11a_merlino.cif'
    ... )
    """

    # Required parameters
    cs_ratio: float = field(default=None)
    ws_ratio: float = field(default=None)

    # Internal state
    _system: AtomicSystem | None = field(default=None, init=False, repr=False)
    _analysis: dict | None = field(default=None, init=False, repr=False)
    _is_built: bool = field(default=False, init=False, repr=False)

    def __post_init__(self):
        """Initialize and validate the builder."""
        super().__init__(None)
        if self.cs_ratio is not None and self.cs_ratio <= 0:
            raise ValueError(f"cs_ratio must be positive, got {self.cs_ratio}")
        if self.ws_ratio is not None and self.ws_ratio < 0:
            raise ValueError(f"ws_ratio must be >= 0, got {self.ws_ratio}")

    def __repr__(self) -> str:
        """Return a styled string representation."""
        lines = ["<CSHBuilder>"]
        lines.append("")
        lines.append("┌─ Configuration")

        cs_display = self.cs_ratio
        if (
            cs_display is None
            and self._analysis is not None
            and "Ca/Si" in self._analysis
        ):
            cs_display = self._analysis["Ca/Si"]

        # Determine ws_ratio to display
        ws_display = self.ws_ratio
        if (
            ws_display is None
            and self._analysis is not None
            and "H2O/Si" in self._analysis
        ):
            ws_display = self._analysis["H2O/Si"]

        lines.append(
            f"│   Ca/Si ratio  : {cs_display:.3f}"
            if cs_display is not None
            else "│   Ca/Si ratio  : None"
        )
        lines.append(
            f"│   H2O/Si ratio : {ws_display:.3f}"
            if cs_display is not None
            else "│   H2O/Si ratio : None"
        )

        if self._system is not None:
            lines.append(f"│   System       : {self._system.num_atoms} atoms")

            # Automatic analysis if it does not exist yet
            if self._analysis is None:
                try:
                    self._analysis = analyze_silicates(self._system)
                except Exception as e:
                    self._analysis = {"error": str(e)}

            if self._analysis is not None and "error" not in self._analysis:
                lines.append("│")
                lines.append("├─ Analysis")
                lines.append(f"│   Ca/(Si+Al)   : {self._analysis['Ca/(Si+Al)']:.3f}")
                lines.append(f"│   Al/Si        : {self._analysis['Al/Si']:.3f}")
                lines.append(f"│   H2O/(Si+Al)  : {self._analysis['H2O/(Si+Al)']:.3f}")
                mcl = self._analysis["MCL"]
                mcl_str = f"{mcl:.2f}" if mcl != float("inf") else "∞"
                lines.append(f"│   MCL          : {mcl_str}")
                qn = self._analysis["Qn_distribution"]
                lines.append(
                    f"│   Qⁿ dist      : Q0={qn[0]:.1f}% Q1={qn[1]:.1f}% Q2={qn[2]:.1f}% Q3={qn[3]:.1f}% Q4={qn[4]:.1f}%"
                )
            elif self._analysis is not None and "error" in self._analysis:
                lines.append("│")
                lines.append("├─ Analysis")
                lines.append(f"│   ⚠️ {self._analysis['error']}")
        else:
            lines.append("│   System       : None")

        lines.append("└─")
        return "\n".join(lines)

    @classmethod
    def from_system(cls, system: AtomicSystem) -> CSHBuilder:
        """
        Create a builder from an existing AtomicSystem.

        This is useful for analyzing or modifying existing C-S-H structures.

        Parameters
        ----------
        system : AtomicSystem
            The existing C-S-H system.

        Returns
        -------
        CSHBuilder
            A builder instance configured with the given system.

        Examples
        --------
        >>> system = AtomicSystem.from_file("csh.lmp")
        >>> builder = CSHBuilder.from_system(system)
        >>> analysis = builder.analyze()
        >>> cash_system, ratio = builder.to_cash(as_ratio=0.1)
        """
        builder = cls()
        builder._system = system
        try:
            builder._analysis = analyze_silicates(system)
        except Exception as e:
            builder._analysis = {"error": str(e)}
        return builder

    def analyze(self, types_map: dict | None = None, cutoff: float = 1.85) -> dict:
        """
        Analyze the silicate network of the current system.

        This method uses `analyze_silicates` to compute:
        - Stoichiometric ratios (Ca/(Si+Al), Al/Si, H2O/(Si+Al))
        - Q^n distribution (polymerization)
        - Mean Chain Length (MCL)

        Parameters
        ----------
        types_map : dict, optional
            Custom mapping of atom types for Si, O, Al, and Ca.
            Defaults to TYPES_PRESET.
        cutoff : float, default=1.85
            Distance cutoff for neighbor selection in Å.

        Returns
        -------
        dict
            Analysis results containing 'Ca/(Si+Al)', 'Al/Si', 'H2O/(Si+Al)',
            'MCL', and 'Qn_distribution'.

        Examples
        --------
        >>> builder = CSHBuilder.from_system(system)
        >>> analysis = builder.analyze()
        >>> print(f"MCL: {analysis['MCL']:.2f}")
        >>> print(f"Qn distribution: {analysis['Qn_distribution']}")
        """
        if self._system is None:
            raise RuntimeError(
                "No system available. Build one first or use from_system()."
            )

        system = self._system

        self._analysis = analyze_silicates(system, types_map=types_map, cutoff=cutoff)
        return self._analysis

    def _update_analysis(self) -> None:
        """Update analysis for the current system."""
        if self._system is not None:
            try:
                self._analysis = analyze_silicates(self._system)
            except Exception as e:
                self._analysis = {"error": str(e)}

    # ------------------------------------------------------------------
    # Classic method (tobermorite)
    # ------------------------------------------------------------------

    def build(
        self,
        supercell: Sequence[int] | None = None,
        model: str = "tob11a_merlino.cif",
        min_mcl: float = 3.0,
        symmetry: bool = True,
        _progress_callback: Callable[[int, str], None] | None = None,
    ) -> AtomicSystem:
        """
        Build C-S-H using the classic tobermorite-based approach.

        Parameters
        ----------
        supercell : Sequence[int], optional
            Supercell dimensions [a, b, c].
            Defaults: [3,5,1] for 'tob11a_hamid.cif', [4,1,1] for 'tob11a_merlino.cif'
        model : str, default='tob11a_merlino.cif'
            Tobermorite model ('tob11a_hamid.cif' or 'tob11a_merlino.cif').
        min_mcl : float, default=3.0
            Minimum mean chain length.
        symmetry : bool, default=True
            Remove bridging silicates symmetrically.

        Returns
        -------
        AtomicSystem
            The built C-S-H structure.

        Examples
        --------
        >>> builder = CSHBuilder(cs_ratio=1.5, ws_ratio=1.0)
        >>> system = builder.build(
        ...     supercell=[4,1,1],
        ...     model='tob11a_merlino.cif'
        ... )
        """

        require_program("packmol")

        from ._interlayer_helpers import (
            distribute_species_in_layers,
            fill_csh_interlayers,
        )
        from ._silicate_helpers import (
            calculate_csh_modifiers,
            find_dense_layer_seam_offset,
            neutralize_csh_charge,
            remove_bridging_silicates,
        )

        # Validate classic-specific parameters
        if model not in ["tob11a_hamid.cif", "tob11a_merlino.cif"]:
            raise ValueError(f"Unsupported model: {model}")
        if min_mcl < 2:
            raise ValueError(f"min_mcl must be >= 2, got {min_mcl}")

        system = AtomicSystem.from_file(STRUCTURES_DIR / model)

        # Apply supercell
        if model == "tob11a_hamid.cif":
            supercell = supercell or [3, 5, 1]
            system.replicate(supercell)
            system.unskew()
        elif model == "tob11a_merlino.cif":
            system.orthogonalize()
            supercell = supercell or [4, 1, 1]
            system.replicate(supercell)

        univ = system.to_mda()
        box = univ.dimensions

        # Re-origin the z-axis so the periodic seam falls inside a dense
        # silicate sheet rather than inside an interlayer gap -- see
        # find_dense_layer_seam_offset for why this matters to Packmol.
        seam_offset = find_dense_layer_seam_offset(univ)
        if seam_offset is not None:
            # AtomGroup.positions returns a copy (fancy indexing via .ix),
            # so it must be reassigned as a whole array -- slice-assigning
            # into the getter's return value silently no-ops.
            positions = univ.atoms.positions.copy()
            positions[:, 2] = (positions[:, 2] - seam_offset) % box[2]
            univ.atoms.positions = positions

        # Calculate modifications
        nsi = len(univ.select_atoms("type Si"))
        nca = len(univ.select_atoms("type Ca"))
        nsi_to_remove, nca_to_add, vacancy_fraction = calculate_csh_modifiers(
            nsi, nca, self.cs_ratio, min_mcl
        )

        if _progress_callback:
            _progress_callback(10, "Removing bridging silicates...")

        univ = remove_bridging_silicates(univ, nsi_to_remove, symmetry)

        # Update Si count
        si_sel = univ.select_atoms("type Si")
        nsi = len(si_sel)

        # Calculate water count
        num_water = round(nsi * self.ws_ratio) + nca_to_add

        if num_water == 0:
            final_system = AtomicSystem.from_mda(univ)
            final_system.set_box(box)
            return final_system

        # Distribute species in layers
        si_layers_pos, nw_layers, nca_layers, nlayers = distribute_species_in_layers(
            num_water,
            nca_to_add,
            si_sel.positions[:, 2],
            supercell[2] if supercell else 1,
            box[2],
        )

        with tempfile.TemporaryDirectory(dir=".") as tmp:
            final_pdb = fill_csh_interlayers(
                univ,
                box,
                si_layers_pos,
                nw_layers,
                nca_layers,
                nca_to_add,
                nlayers,
                tmp,
                _progress_callback,
            )

            if nca_to_add != 0:
                neutralize_csh_charge(final_pdb)

            final_system = AtomicSystem.from_file(final_pdb)

        final_system.set_box(box)

        if _progress_callback:
            _progress_callback(100, "Structure complete")

        self._system = final_system
        self._update_analysis()
        return final_system

    # ------------------------------------------------------------------
    # Conversion to C-A-S-H
    # ------------------------------------------------------------------

    def to_cash(
        self, as_ratio: float, supercell_factor: int = 1
    ) -> tuple[AtomicSystem, float, dict | None, dict | None]:
        """Convert C-S-H to C-A-S-H by substituting Si with Al."""
        if self._system is None:
            raise RuntimeError(
                "No system available. Build one first or use from_system()."
            )

        system = self._system

        universe = system.to_mda()
        si_sel = universe.select_atoms("type Si")

        if len(si_sel) == 0:
            raise ValueError("No silicon atoms found in the system.")

        nsi = len(si_sel)
        nal_target = round(as_ratio * nsi)

        # `substitute_si_by_al`'s 3rd positional parameter is `symmetry`
        # (bool), not a supercell factor -- passing `supercell_factor`
        # there silently mapped any nonzero value to `symmetry=True`
        # regardless of its actual value. There is currently no
        # supercell-aware substitution logic to forward it to.
        universe, substituted_ids = substitute_si_by_al(
            universe, nal_target, symmetry=True
        )

        system.set_type2atoms(substituted_ids, "Al")
        system.keep_connection_types(
            bond_types=["Hsi-Osih", "Ha-Oah", "Hh-Oh", "Hw-Ow"],
            angle_types=["Hw-Ow-Hw"],
        )

        self._system = system
        self._update_analysis()
        return system


@dataclass
class AFBuilder(BaseBuilder):
    """
    Builder for creating AFt (Ettringite) or AFm structures.

    This builder creates AFt/AFm structures by adding water molecules to
    a sulfate-containing framework.

    Parameters
    ----------
    ws_ratio : float
        Water/Sulfate ratio.
    progress_callback : Callable, optional
        Progress callback function.

    Examples
    --------
    >>> # Create AFt
    >>> builder = AFBuilder(ws_ratio=10.0)
    >>> system = builder.build_aft(
    ...     supercell=[2,2,1]
    ... )
    >>>
    >>> # Create AFm
    >>> builder = AFBuilder(ws_ratio=8.0)
    >>> system = builder.build_afm(
    ...     supercell=[3,1,1]
    ... )
    """

    # Required parameters
    ws_ratio: float = field(default=None)

    # Optional parameter (common)
    progress_callback: Callable[[int, str], None] | None = None

    # Internal state
    _result: AtomicSystem | None = field(default=None, init=False, repr=False)

    def __post_init__(self):
        """Initialize and validate the builder."""
        super().__init__(None)
        if self.ws_ratio is None or self.ws_ratio < 0:
            raise ValueError(f"ws_ratio must be >= 0, got {self.ws_ratio}")

    def __repr__(self) -> str:
        """Return a styled string representation."""
        lines = ["<AFBuilder>"]
        lines.append("")
        lines.append("┌─ Configuration")
        lines.append(f"│   ws_ratio   : {self.ws_ratio:.3f}")
        lines.append("└─")
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Build methods
    # ------------------------------------------------------------------

    def build(
        self,
        structure_type: str = "aft",
        supercell: Sequence[int] | None = None,
        model_file: str | None = None,
        **kwargs,
    ) -> AtomicSystem:
        """Build an AFt or AFm structure, dispatching on ``structure_type``.

        Satisfies the ``BaseBuilder`` interface; prefer :meth:`build_aft`
        or :meth:`build_afm` directly when the structure type is known.
        """
        if structure_type == "aft":
            return self.build_aft(
                supercell=supercell, model_file=model_file or "aft_moore.cif"
            )
        elif structure_type == "afm":
            return self.build_afm(
                supercell=supercell, model_file=model_file or "afm.cif", **kwargs
            )
        else:
            raise ValueError(
                f"Unknown structure_type: {structure_type!r}. Expected 'aft' or 'afm'."
            )

    def build_aft(
        self, supercell: Sequence[int] | None = None, model_file: str = "aft_moore.cif"
    ) -> AtomicSystem:
        """
        Build AFt (Ettringite) structure.

        Parameters
        ----------
        supercell : Sequence[int], optional
            Supercell dimensions [a, b, c].
        model_file : str, default='aft_moore.cif'
            CIF file for the AFt structure.

        Returns
        -------
        AtomicSystem
            The built AFt structure.

        Examples
        --------
        >>> builder = AFBuilder(ws_ratio=10.0)
        >>> system = builder.build_aft(supercell=[2,2,1])
        """
        return self._build_af("aft", supercell, model_file)

    def build_afm(
        self,
        supercell: Sequence[int] | None = None,
        model_file: str = "afm.cif",
        _progress_callback: Callable[[int, str], None] | None = None,
    ) -> AtomicSystem:
        """
        Build AFm structure.

        Parameters
        ----------
        supercell : Sequence[int], optional
            Supercell dimensions [a, b, c].
        model_file : str, default='afm.cif'
            CIF file for the AFm structure.

        Returns
        -------
        AtomicSystem
            The built AFm structure.

        Examples
        --------
        >>> builder = AFBuilder(ws_ratio=8.0)
        >>> system = builder.build_afm(supercell=[3,1,1])
        """

        require_program("packmol")

        return self._build_af("afm", supercell, model_file, _progress_callback)

    def _build_af(
        self,
        structure_type: str,
        supercell: Sequence[int] | None = None,
        model_file: str = "aft_moore.cif",
        _progress_callback: Callable[[int, str], None] | None = None,
    ) -> AtomicSystem:
        """Internal method to build AFt/AFm structures."""
        require_program("packmol")

        from .._packmol import (
            PackmolInput,
            PackmolStructure,
            get_structure_path,
            run_packmol,
            write_mda_pdb,
        )

        if supercell is not None and len(supercell) != 3:
            raise ValueError(f"supercell must have 3 elements, got {len(supercell)}")

        # Load base structure
        input_system = AtomicSystem.from_file(STRUCTURES_DIR / model_file)

        if supercell is not None:
            input_system.replicate(supercell)

        input_system.orthogonalize()

        # Get information with MDAnalysis
        univ = input_system.to_mda()
        box = univ.dimensions

        # Count sulfur atoms
        s_sel = univ.select_atoms("type S")
        num_s = len(s_sel)
        num_water = round(num_s * self.ws_ratio)

        # Get reference atom for Packmol
        sel = univ.select_atoms("all")
        idmin = np.argmin(sel.positions[:, 2]) + 1

        if _progress_callback:
            _progress_callback(10, f"Building {structure_type.upper()}...")

        with tempfile.TemporaryDirectory(dir=".") as tmp:
            tmp_path = Path(tmp)
            tempfile_ipdb = tmp_path / "itmp.pdb"
            write_mda_pdb(sel, tempfile_ipdb)

            # get_structure_path may resolve to a .lt (moltemplate) file,
            # which Packmol cannot read directly; load it as an
            # AtomicSystem so run_packmol rewrites it to PDB first.
            h2o_path = get_structure_path("h2o", tmp)
            h2o = AtomicSystem.from_file(h2o_path)

            inside_box = (0.0, 0.0, 0.0, box[0], box[1], box[2])

            structures = [
                PackmolStructure(
                    structure=str(tempfile_ipdb),
                    number=1,
                    inside_box=inside_box,
                    extra_instructions=[
                        f"atoms {idmin}",
                        "fixed 0 0 0 0 0 0",
                        "end atoms",
                    ],
                ),
                PackmolStructure(
                    structure=h2o,
                    number=num_water,
                    inside_box=inside_box,
                ),
            ]

            if _progress_callback:
                _progress_callback(50, f"Adding {num_water} water molecules...")

            packmol_input = PackmolInput(
                tolerance=2.0,
                output="output.pdb",
                filetype="pdb",
                structures=structures,
            )
            output_system = run_packmol(packmol_input)

        output_system.set_box(box)

        if _progress_callback:
            _progress_callback(100, f"{structure_type.upper()} complete")

        self._result = output_system
        return self._result
