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
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ...._constants import CHARGES_DICT, MASSES_DICT
from .base import BaseReader

if TYPE_CHECKING:
    from pymatgen.core import Structure


class PMGReader(BaseReader):
    """Read from Pymatgen Structure."""

    @classmethod
    def read(cls, structure: Structure, refine: bool = True) -> dict:
        """
        Convert a Pymatgen Structure to a topology dictionary.

        Parameters
        ----------
        structure : Structure
            Pymatgen Structure object
        refine : bool
            If True, refine the structure using SpacegroupAnalyzer.

        Returns
        -------
        dict
            Topology dictionary ready for AtomicSystem initialization
        """
        from itertools import permutations

        from pymatgen.core import Lattice
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        # Store initial structure parameters
        original_abc = list(structure.lattice.abc)

        # Refine if requested
        if refine:
            try:
                analyzer = SpacegroupAnalyzer(structure)
                refined_structure = analyzer.get_refined_structure()
            except Exception as e:
                print(f"Warning: Could not refine structure: {e}")
                refined_structure = structure
        else:
            refined_structure = structure

        # Extract parameters from the refined structure
        abc = list(refined_structure.lattice.abc)
        angles = list(refined_structure.lattice.angles)
        positions = refined_structure.cart_coords.copy()

        # Reindex the axes of the refined structure
        # to match the order of the initial structure
        best_mapping = [0, 1, 2]
        min_diff = float("inf")

        for p in permutations([0, 1, 2]):
            # Compare the order of the axes of the refined structure with the original
            current_diff = sum(abs(abc[i] - original_abc[j]) for j, i in enumerate(p))
            if current_diff < min_diff:
                min_diff = current_diff
                best_mapping = list(p)

        # Apply reindexing if order has changed
        if best_mapping != [0, 1, 2]:
            # Permute the lattice vectors and the fractional coordinates together,
            # so that the crystal itself is left unchanged (only the labels a, b, c
            # move). An odd permutation makes the lattice left-handed: flip the
            # last axis to stay right-handed and avoid mirroring the structure.
            matrix = refined_structure.lattice.matrix[best_mapping].copy()
            frac = refined_structure.frac_coords[:, best_mapping].copy()
            if np.linalg.det(matrix) < 0:
                matrix[2] *= -1
                frac[:, 2] *= -1
            permuted = Lattice(matrix)
            reindexed_abc = list(permuted.abc)
            reindexed_angles = list(permuted.angles)

            # Express the positions in the standard orientation of the new box
            # (a along x, b in the xy plane), as expected by the box parameters.
            standard = Lattice.from_parameters(*reindexed_abc, *reindexed_angles)
            positions = standard.get_cartesian_coords(frac)

            warnings.warn(
                f"Axes reindexed with mapping {best_mapping}: "
                f"original {original_abc} -> refined {abc} -> reindexed {reindexed_abc}",
                category=UserWarning,
                stacklevel=2,
            )
            abc = reindexed_abc
            angles = reindexed_angles

        # Extract Atom Types
        types = [site.species.elements[0].name for site in refined_structure]
        unique_types = sorted(set(types))

        # Create the DataFrames
        ids = np.arange(1, len(positions) + 1)
        charges = np.array([CHARGES_DICT.get(t, 0.0) for t in types])

        df_atoms = pd.DataFrame(
            {
                "id": ids,
                "type": types,
                "charge": charges,
                "x": positions[:, 0],
                "y": positions[:, 1],
                "z": positions[:, 2],
            }
        )

        # Force data types
        df_atoms["id"] = df_atoms["id"].astype(int)
        df_atoms[["charge", "x", "y", "z"]] = df_atoms[
            ["charge", "x", "y", "z"]
        ].astype(float)
        df_atoms.set_index("id", inplace=True)

        # Masses and charges per type
        masses_dict = {t: MASSES_DICT.get(t, 1.0) for t in unique_types}
        charges_dict = {t: CHARGES_DICT.get(t, 0.0) for t in unique_types}

        # Assemble the topology
        topology = {
            "box": abc + angles,
            "masses": masses_dict,
            "charges": charges_dict,
            "atoms": df_atoms,
            "bonds": None,
            "angles": None,
            "dihedrals": None,
            "impropers": None,
            "velocities": None,
            "atom_style": "full",
        }

        # Store the original structure (Miller indices of surfaces refer to its axes)
        topology["_pmg_struct"] = structure

        return topology
