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
from typing import TYPE_CHECKING, Self

import numpy as np

from .._constants import MASSES_DICT
from ._format import BoxFormat, lattice2vectors, normalize_box

if TYPE_CHECKING:
    from .atomic_system import AtomicSystem


class EditMixin:
    """Mixin class containing editing and modification methods
    for the AtomicSystem class.
    """

    def add_atoms(
        self,
        atypes: list[str | int] | np.ndarray,
        positions: list[list[float]] | np.ndarray,
        charges: list[float] | np.ndarray | None = None,
        masses: list[float] | np.ndarray | None = None,
    ) -> None:
        """Add multiple atoms to the system in a single operation.

        Parameters
        ----------
        atypes : list of str or int, or np.ndarray
            Atom types (e.g. ``['H', 'Ow', 'Ca']``).
        positions : array-like of shape (n, 3)
            Cartesian coordinates in Angstroms.
        charges : array-like of shape (n,), optional
            Partial charges in elementary charge units. Default is ``0.0``
            for each atom.
        masses : list of float, optional
            Atomic masses in g/mol. If ``None``, masses are looked up in
            ``MASSES_DICT`` or default to ``1.0`` for unknown types.
        """

        import pandas as pd

        atypes = list(atypes)
        positions = np.asarray(positions, dtype=float)

        if positions.ndim == 1:
            positions = positions.reshape(1, 3)

        n = len(atypes)
        charges = np.zeros(n) if charges is None else np.asarray(charges, dtype=float)
        masses = [None] * n if masses is None else list(masses)

        new_id = 1 if self.atoms.empty else self.atoms.index.max() + 1
        new_ids = np.arange(new_id, new_id + n)

        new_rows = []
        for atype, pos, charge, mass in zip(atypes, positions, charges, masses):
            new_rows.append(
                {
                    "type": str(atype),
                    "x": pos[0],
                    "y": pos[1],
                    "z": pos[2],
                    "charge": float(charge),
                }
            )
            if atype not in self._masses:
                if mass is not None:
                    self._masses[atype] = float(mass)
                elif atype in MASSES_DICT:
                    self._masses[atype] = MASSES_DICT[atype]
                else:
                    self._masses[atype] = 1.0

        new_df = pd.DataFrame(new_rows, index=new_ids)
        self.atoms = pd.concat([self.atoms, new_df])

    def add_atom(
        self,
        atype: str | int,
        position: list[float],
        charge: float = 0.0,
        mass: float = None,
    ) -> None:
        """Add a single atom. See :meth:`add_atoms` for batch insertion."""
        self.add_atoms([atype], [position], [charge], [mass])

    def protonate_atoms(
        self, atom_indices: list[int], bond_length: float = 1.0
    ) -> None:
        """Add protons to multiple atoms in a single operation.

        Parameters
        ----------
        atom_indices : list[int]
            Atom indices, as used everywhere else in this class (i.e. the
            `atoms` DataFrame's label index / atom id) -- not a 0-based
            positional row offset.
        """

        u = self.to_mda()
        p_mass = MASSES_DICT.get("H", 1.008)

        new_atoms = []

        for atom_index in atom_indices:
            target_atom = self.atoms.loc[atom_index]
            pos_target = target_atom[["x", "y", "z"]].values.astype(float)

            # `to_mda()` tags atoms with an `ids` TopologyAttr matching
            # `self.atoms.index`, so select by that real id directly rather
            # than MDAnalysis's own 0-based positional "index" -- using the
            # latter here previously mismatched `.loc` by one atom (or
            # raised IndexError for the last atom in the system).
            neighbors = u.select_atoms(f"around 2.2 id {atom_index}")
            if len(neighbors) > 0:
                direction = pos_target - neighbors.center_of_mass()
            else:
                direction = np.array([0.0, 0.0, 1.0])

            norm = np.linalg.norm(direction)
            direction = direction / norm if norm > 1e-5 else np.array([0.0, 0.0, 1.0])

            new_atoms.append(
                {
                    "atype": "H",
                    "position": pos_target + direction * bond_length,
                    "charge": 1.0,
                    "mass": p_mass,
                }
            )

        for atom in new_atoms:
            self.add_atom(**atom)

    def protonate_atom(self, atom_index: int, bond_length: float = 1.0) -> AtomicSystem:
        """Add a proton to a single atom."""
        self.protonate_atoms([atom_index], bond_length)

    def remove_atoms(self, indices: Sequence[int] | int) -> None:
        """Remove specified atoms from the system and update topology accordingly.

        Parameters
        ----------
        indices : Sequence[int] or int
            Index or sequence of indices of the atoms to remove.

        Notes
        -----
        This method automatically reindexes the remaining atoms, updates velocities,
        removes unused atom types from internal mass and charge storage, and updates
        or removes any bonds, angles, dihedrals, or impropers containing the removed atoms.
        """
        # Normalize input to a list of integers
        if isinstance(indices, int):
            target_indices = [indices]
        else:
            target_indices = list(indices)

        if not hasattr(self, "atoms") or self.atoms is None or self.atoms.empty:
            warnings.warn("System contains no atoms to remove.", UserWarning)
            return

        # Filter out indices that do not exist in the current system
        valid_indices = [idx for idx in target_indices if idx in self.atoms.index]

        if not valid_indices:
            warnings.warn(
                f"None of the target indices {target_indices} were found in the system.",
                UserWarning,
            )
            return

        # Store unique atom types before removal
        old_types = set(self.atom_types)

        # Remove selected atoms from the DataFrame
        df_atoms = self.atoms.drop(index=valid_indices).copy()

        # Create mapping from old atom IDs to new contiguous IDs (1 to N)
        old_ids = df_atoms.index
        new_ids = np.arange(1, len(df_atoms) + 1)
        df_atoms.index = new_ids
        self.atoms = df_atoms

        # Synchronize velocities DataFrame if present
        if hasattr(self, "velocities") and self.velocities is not None:
            df_vel = self.velocities.drop(index=valid_indices).copy()
            df_vel.index = new_ids
            self.velocities = df_vel

        # Clean up unused atom types from internal storage (dictionaries)
        new_types = set(self.atom_types)
        removed_types = old_types - new_types

        for atype in removed_types:
            if hasattr(self, "_masses") and isinstance(self._masses, dict):
                self._masses.pop(atype, None)
            if hasattr(self, "_charges") and isinstance(self._charges, dict):
                self._charges.pop(atype, None)

        # Update topology tables (bonds, angles, dihedrals, impropers)
        id_map = dict(zip(old_ids, new_ids))

        topology_names = [
            ("bonds", 2),
            ("angles", 3),
            ("dihedrals", 4),
            ("impropers", 4),
        ]

        for name, n_cols in topology_names:
            df = getattr(self, name, None)
            if df is None or df.empty:
                continue

            atom_cols = [f"atom_{i}" for i in range(1, n_cols + 1)]

            # Filter out any interaction involving at least one removed atom
            mask_keep = ~df[atom_cols].isin(valid_indices).any(axis=1)
            df = df.loc[mask_keep].copy()

            if df.empty:
                setattr(self, name, None)
                continue

            # Remap old atom IDs to new IDs
            for col in atom_cols:
                df[col] = df[col].map(id_map)

            # Reset row index from 1 to N
            df.index = np.arange(1, len(df) + 1)
            setattr(self, name, df)

        # Clear internal cache

    def remove_atom(self, index: int) -> None:
        """Remove a single atom. See :meth:`remove_atoms` for batch removal."""
        self.remove_atoms([index])

    def set_box(self, new_box: Sequence[float] | np.ndarray) -> None:
        """Assign a new box to the system and compute all representations.

        Parameters
        ----------
        new_box : Sequence[float] or np.ndarray
            The box in any valid format (Lattice parameters, LAMMPS bounds/tilts, or 3x3 matrix).
        """
        # Clear the cache since the box geometry changes

        # Conversion/Normalization vers les 3 formats internes
        self._box = normalize_box(new_box, target=BoxFormat.LATTICE)
        self._box_lmp = normalize_box(new_box, target=BoxFormat.LAMMPS)
        self._box_vectors = normalize_box(new_box, target=BoxFormat.VECTORS)

    def set_atom_position(self, index: int, position: Sequence[float]) -> AtomicSystem:
        """Modify the coordinates of a single atom using a vector (x, y, z).

        Parameters
        ----------
        index : int
            The index of the atom to modify.
        position : Sequence[float]
            The new (x, y, z) coordinates.

        Returns
        -------
        AtomicSystem
            The modified system.
        """
        if len(position) != 3:
            raise ValueError("Position must be an iterable of 3 floats (x, y, z)")

        if index not in self.atoms.index:
            print(f"Warning: Index {index} not found.")
            return self

        x, y, z = position

        self.atoms.at[index, "x"] = float(x)
        self.atoms.at[index, "y"] = float(y)
        self.atoms.at[index, "z"] = float(z)

        return self

    def replicate(self, factors: Sequence[int]) -> Self:
        """Replicate the system by integer factors along the lattice vectors (a, b, c).

        Parameters
        ----------
        factors : Sequence[int]
            Number of copies along each lattice vector (nx, ny, nz).

        Returns
        -------
        Self
            The replicated system.
        """

        import pandas as pd

        nx, ny, nz = factors
        if nx == 1 and ny == 1 and nz == 1:
            return self

        # Prepare the basis vectors
        # lattice2vectors returns a (3, 3) matrix where each row is a vector v1, v2, v3
        vecs = np.array(self._box_vectors)
        v1, v2, v3 = vecs[0], vecs[1], vecs[2]

        original_atoms = self.atoms.copy()
        original_num_atoms = len(original_atoms)

        all_atoms = []

        # Dictionaries to store new interactions
        new_interactions = {"bonds": [], "angles": [], "dihedrals": [], "impropers": []}

        # Replication loop
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if i == 0 and j == 0 and k == 0:
                        all_atoms.append(original_atoms)
                        continue

                    # Calculate the shift for this cell
                    shift = i * v1 + j * v2 + k * v3

                    # Copy atoms
                    new_atoms = original_atoms.copy()
                    new_atoms["x"] += shift[0]
                    new_atoms["y"] += shift[1]
                    new_atoms["z"] += shift[2]

                    # Updating indices (IDs) to avoid duplicates
                    # shift the indexes by (number_atoms_orig *unique_multiplier)
                    offset = (i * ny * nz + j * nz + k) * original_num_atoms
                    new_atoms.index += offset
                    all_atoms.append(new_atoms)

                    # Replicate the interactions (Bonds, Angles, etc.)
                    for name in ["bonds", "angles", "dihedrals", "impropers"]:
                        df = getattr(self, name)
                        if df is not None and not df.empty:
                            new_df = df.copy()
                            # shift the index type of the interaction itself
                            new_df.index += (i * ny * nz + j * nz + k) * len(df)

                            # shift the columns atom_1, atom_2, etc.
                            atom_cols = [c for c in df.columns if c.startswith("atom_")]
                            for col in atom_cols:
                                new_df[col] += offset

                            new_interactions[name].append(new_df)

        # Merge Dataframes
        self.atoms = pd.concat(all_atoms)

        if self.velocities is not None:
            self.velocities = None

        for name in ["bonds", "angles", "dihedrals", "impropers"]:
            if new_interactions[name]:
                full_list = [getattr(self, name)] + new_interactions[name]
                setattr(self, name, pd.concat(full_list))

        new_box = self.box.copy()
        new_box[0] *= nx
        new_box[1] *= ny
        new_box[2] *= nz
        self.set_box(new_box)

        return self

    def wrap(self) -> None:
        """Wraps all atoms back into the simulation box [0, L] using Periodic Boundary Conditions.

        Handles both orthogonal and triclinic boxes.
        """
        if self.atoms.empty:
            return self

        # Get the pass matrix (H_matrix)
        h_matrix = np.array(self._box_vectors)
        inv_h_matrix = np.linalg.inv(h_matrix)

        # Extract current coordinates (N, 3)
        coords = self.atoms[["x", "y", "z"]].values

        # Conversion to fractional coordinates (0 to 1)
        # s = r . inv(H)
        frac_coords = np.dot(coords, inv_h_matrix)

        # Application of modulo 1.0
        # This brings everything into the range [0, 1[
        frac_coords %= 1.0

        # Return to Cartesian coordinates
        # r_new = s_new . H
        new_coords = np.dot(frac_coords, h_matrix)

        self.atoms[["x", "y", "z"]] = new_coords

        return self

    def orthogonalize(self, tolerance: float = 0.1, max_replica: int = 10) -> bool:
        """
        Attempt to transform a triclinic (skewed) box into an orthogonal one.

        This method searches for new lattice vectors that align with the Cartesian
        axes (X, Y, Z) by exploring linear combinations of the current basis
        vectors. It then re-populates the new box by shifting and wrapping
        existing atom coordinates.

        Parameters
        ----------
        tolerance
            The maximum deviation allowed from the Cartesian axes (in Å) for
            a candidate vector to be considered orthogonal.
        max_replica
            The search range for linear combinations (m, n, o) of the original
            lattice vectors. Higher values increase the chance of finding an
            orthogonal cell but significantly increase computation time.

        Raises
        ------
        ValueError
            If no orthogonal representation could be found within ``max_replica``
            replications.

        Notes
        -----
        - This operation will change the total number of atoms in the system.
        - Velocities and interactions (bonds, etc.) are typically lost or
          invalidated during this specific geometric transformation.
        - The resulting box will have tilt factors (xy, xz, yz) set to zero.
        """

        import pandas as pd

        h_matrix = np.array(lattice2vectors(self.box))
        new_vectors = []

        for i in range(3):
            best_v = None
            best_len = float("inf")

            for m in range(-max_replica, max_replica + 1):
                for n in range(-max_replica, max_replica + 1):
                    for o in range(-max_replica, max_replica + 1):
                        if m == 0 and n == 0 and o == 0:
                            continue

                        v_cand = m * h_matrix[0] + n * h_matrix[1] + o * h_matrix[2]
                        # check the alignment with the i axis
                        others = [v_cand[j] for j in range(3) if j != i]

                        if all(abs(comp) < tolerance for comp in others):
                            v_len = abs(v_cand[i])  # We want the length on the axis
                            if 0.1 < v_len < best_len:  # 0.1 to avoid zero vector
                                best_len = v_len
                                best_v = v_cand

            if best_v is None:
                raise ValueError(
                    f"Could not find an orthogonal representation for axis {i} "
                    f"within {max_replica} replications. "
                    f"Try increasing 'max_replica' or use 'unskew()' instead."
                )
            new_vectors.append(best_v)

        diag_dim = [
            abs(new_vectors[0][0]),
            abs(new_vectors[1][1]),
            abs(new_vectors[2][2]),
        ]
        new_h_matrix = np.diag(diag_dim)

        search_range = np.arange(-max_replica, max_replica + 1)
        m_grid, n_grid, o_grid = np.meshgrid(search_range, search_range, search_range)
        translations = np.vstack([m_grid.ravel(), n_grid.ravel(), o_grid.ravel()]).T
        translation_vectors = translations @ h_matrix

        all_replicas = []
        eps = 1e-5

        for vec in translation_vectors:
            temp_df = self.atoms.copy()
            temp_df["x"] += vec[0]
            temp_df["y"] += vec[1]
            temp_df["z"] += vec[2]

            mask = (
                (temp_df["x"] >= -eps)
                & (temp_df["x"] < diag_dim[0] - eps)
                & (temp_df["y"] >= -eps)
                & (temp_df["y"] < diag_dim[1] - eps)
                & (temp_df["z"] >= -eps)
                & (temp_df["z"] < diag_dim[2] - eps)
            )
            if mask.any():
                all_replicas.append(temp_df[mask])

        if not all_replicas:
            return False

        self.atoms = pd.concat(all_replicas, ignore_index=True)

        self.atoms.index = range(1, len(self.atoms) + 1)
        if "id" in self.atoms.columns:
            self.atoms["id"] = range(1, len(self.atoms) + 1)

        self.set_box(tuple(new_h_matrix))

    def unskew(self) -> None:
        """
        Unskew the box if the non-diagonal parameters of the box matrix
        :math:`|H_{i,j}| > H_{j,j}` for :math:`i \\neq j` by adding
        :math:`\\pm H_{j,:}` to :math:`H_{i,:}` while :math:`|H_{i,j}| > H_{j,j}`.
        """

        boxm = np.array(lattice2vectors(self.box))

        for i in range(3):
            for j in range(3):
                if i != j:
                    while boxm[i, j] > 0.5 * boxm[i, i]:
                        boxm[i] -= boxm[j]
                    while boxm[i, j] < -0.5 * boxm[i, i]:
                        boxm[i] += boxm[j]

        self.set_box((boxm[0], boxm[1], boxm[2]))
        self.wrap()

    def center_on_com(self, atom_types: list[str] | None = None) -> None:
        """Translates the system so its center of mass is at the center of the box."""

        u = self.to_mda()

        if atom_types is not None:
            type_str = " ".join(atom_types)
            sel = u.select_atoms(f"type {type_str}")
        else:
            sel = u.atoms

        com = sel.center_of_mass(pbc=True)  # Bai & Breen method, PBC-aware, triclinic-compatible

        center = u.dimensions[:3] / 2
        shift = center - com

        self.atoms[["x", "y", "z"]] += shift

        self.wrap()
