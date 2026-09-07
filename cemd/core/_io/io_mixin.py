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

from pathlib import Path
from typing import Any, Self

import numpy as np


class IOMixin:
    """Mixin for file I/O operations."""

    @classmethod
    def from_file(cls, path: str | Path, **kwargs) -> Self:
        """
        Load a system from a file.

        Parameters
        ----------
        path : str
            Path to the file (supports .cif, .data, .pdb, etc.)
        **kwargs
            Additional arguments passed to the specific reader:

            - For CIF files: ``primitive`` (bool, default ``False``) to
              use the primitive cell; ``refine`` (bool, default ``True``)
              to refine the structure using SpacegroupAnalyzer.
            - For LAMMPS data files: ``atom_style`` (str, default
              ``"full"``).
        """

        file_path = Path(path)
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        ext = file_path.suffix.lower()

        readers = {
            ".data": cls._from_lammps_data,
            ".lmp": cls._from_lammps_data,
            ".cif": cls._from_cif,
            ".pdb": cls._from_pdb,
            ".lt": cls._from_lammpstemplate,
            ".sdf": cls._from_sdf,
        }

        if ext not in readers:
            raise ValueError(f"Unsupported file format: {ext}")

        topology = readers[ext](path, **kwargs)
        return cls(topology)

    @classmethod
    def from_mda(cls, obj) -> Self:
        """Create from MDAnalysis Universe or AtomGroup."""
        topology = cls._from_mda(obj)
        return cls(topology)

    @classmethod
    def from_pmg(cls, struct, refine=True) -> Self:
        """
        Create from Pymatgen Structure.

        Parameters
        ----------
        structure : Structure
            Pymatgen Structure object
        refine : bool
            If True, refine the structure using SpacegroupAnalyzer.
        """
        topology = cls._from_pmg(struct, refine)
        return cls(topology)

    @classmethod
    def from_smiles(cls, smiles: str) -> Self:
        """Create from SMILES string."""
        topology = cls._from_smiles(smiles)
        return cls(topology)

    # ========================================================================
    # Database exploration methods
    # ========================================================================

    @classmethod
    def from_cod(cls, cod_id: int | str | None = None) -> Self:
        """
        Load a structure from the Crystallography Open Database.

        Parameters
        ----------
        cod_id : int or str, optional
            COD entry to download. Given one, the structure is fetched
            straight away; omitted, an interactive explorer opens instead.
            Pass the identifier in a script: the explorer cannot be
            replayed, and a script that names its structure says where its
            model came from.

        Returns
        -------
        AtomicSystem or None
            The structure, or None if the explorer was cancelled.

        Examples
        --------
        >>> calcite = AtomicSystem.from_cod(9016705)   # reproducible
        >>> system = AtomicSystem.from_cod()           # browse instead
        """
        from .sources.cod import explore_cod, get_structure_by_cod_id

        if cod_id is not None:
            return get_structure_by_cod_id(int(cod_id))
        return explore_cod()

    @classmethod
    def from_pubchem(cls, cid: int | str | None = None) -> Self:
        """
        Load a molecule from PubChem.

        Parameters
        ----------
        cid : int or str, optional
            PubChem compound identifier to download. Given one, the
            molecule is fetched straight away; omitted, an interactive
            explorer opens instead.

        Returns
        -------
        AtomicSystem or None
            The molecule, or None if it could not be fetched or the
            explorer was cancelled.

        Examples
        --------
        >>> caffeine = AtomicSystem.from_pubchem(2519)  # reproducible
        >>> molecule = AtomicSystem.from_pubchem()      # browse instead

        Notes
        -----
        Only coordinates and bonds are returned. Call
        :meth:`~cemd.core.atomic_system.AtomicSystem.guess_connections` for
        the angles and dihedrals.
        """
        from .sources.pubchem import explore_pubchem, get_structure

        if cid is not None:
            return get_structure(int(cid))
        return explore_pubchem()

    # ========================================================================
    # Write methods
    # ========================================================================

    def write(
        self, path: str | Path, atom_style: str = "full", oldstyle: bool = False
    ) -> None:
        """
        Write the system to a file.

        Parameters
        ----------
        path : str or Path
            Path to the output file (supports .data, .pdb, etc.)
        atom_style : str, default='full'
            LAMMPS atom style.
        oldstyle : bool, default=False
            Use old style formatting if applicable.
        """
        file_path = Path(path)

        ext = file_path.suffix.lower()

        writers = {
            ".data": self._write_lammps_data,
            ".pdb": self._write_pdb,
        }

        if ext not in writers:
            raise ValueError(f"Unsupported output format: {ext}")

        writers[ext](str(file_path), atom_style=atom_style, oldstyle=oldstyle)

    # ========================================================================
    # Converters
    # ========================================================================

    def to_mda(self) -> Any:
        """
        Convert this AtomicSystem instance to an MDAnalysis Universe.

        Maps per-type properties (masses, elements) to individual atoms
        as required by MDAnalysis topology attributes.

        Returns
        -------
        mda.Universe
            An MDAnalysis Universe populated with atom types, names,
            charges, coordinates, masses, and elements.
        """
        import MDAnalysis as mda

        universe = mda.Universe.empty(self.num_atoms, trajectory=True)

        atom_types = self.atoms["type"].to_numpy()

        universe.add_TopologyAttr("type", atom_types)
        universe.add_TopologyAttr("name", atom_types)
        universe.add_TopologyAttr("charge", self.atoms["charge"].to_numpy())
        universe.add_TopologyAttr("ids", self.atoms.index.to_numpy())

        elem_map = self.elements  # Dictionnaire {type_id: symbole}
        atom_elements = np.array(
            [elem_map.get(t, elem_map.get(str(t), "X")) for t in atom_types],
            dtype=str,
        )
        universe.add_TopologyAttr("elements", atom_elements)

        mass_map = self.masses
        atom_masses = np.array(
            [mass_map.get(t, mass_map.get(str(t), 1.0)) for t in atom_types],
            dtype=float,
        )
        universe.add_TopologyAttr("masses", atom_masses)

        universe.atoms.positions = self.atoms[["x", "y", "z"]].to_numpy()
        universe.dimensions = self.box
        self._add_connectivity_to_mda(universe)

        return universe

    def to_pmg(self) -> Any:
        """Convert to Pymatgen Structure."""
        from pymatgen.core import Lattice, Structure

        lattice = Lattice.from_parameters(*self.box)
        structure = Structure(
            lattice=lattice,
            species=self.atoms["type"].tolist(),
            coords=self.atoms[["x", "y", "z"]].to_numpy(),
            coords_are_cartesian=True,
        )
        return structure

    # ========================================================================
    # Private readers
    # ========================================================================

    @staticmethod
    def _from_lammps_data(path: str) -> dict:
        """Read LAMMPS data file."""
        from .formats.lammps import LAMMPSReader

        return LAMMPSReader.read(path)

    @staticmethod
    def _from_cif(path: str, primitive: bool = False, refine: bool = True) -> dict:
        """Read CIF file."""
        from .formats.cif import CIFReader

        return CIFReader.read(path, primitive, refine)

    @staticmethod
    def _from_pdb(path: str) -> dict:
        """Read PDB file."""
        from .formats.pdb import PDBReader

        return PDBReader.read(path)

    @staticmethod
    def _from_lammpstemplate(path: str) -> dict:
        """Read LT file."""
        from .formats.lt import LTReader

        return LTReader.read(path)

    @staticmethod
    def _from_sdf(path: str) -> dict:
        """Read SDF file."""
        from .formats.sdf import SDFReader

        return SDFReader.read(path)

    @staticmethod
    def _from_mda(obj) -> dict:
        """Read from MDAnalysis."""
        from .formats.mda import MDAReader

        return MDAReader.read(obj)

    @staticmethod
    def _from_pmg(struct, refine: bool = True) -> dict:
        """Read from Pymatgen."""
        from .formats.pmg import PMGReader

        return PMGReader.read(struct, refine)

    @staticmethod
    def _from_smiles(smiles: str) -> dict:
        """Read from SMILES."""
        from .formats.smiles import SMILESReader

        return SMILESReader.read(smiles)

    # ========================================================================
    # Private writers
    # ========================================================================

    def _write_lammps_data(
        self, path: str, atom_style: str = "full", oldstyle: bool = False
    ) -> None:
        """Write to LAMMPS data file."""
        from .formats.lammps import LAMMPSWriter

        LAMMPSWriter.write(self, path, atom_style=atom_style, oldstyle=oldstyle)

    def _write_pdb(self, path: str, **kwargs) -> None:
        """Write to PDB file."""
        from .formats.pdb import PDBWriter

        PDBWriter.write(self, path)

    def _add_connectivity_to_mda(self, universe) -> None:
        """Add bonds, angles, dihedrals, impropers to MDAnalysis Universe."""
        connectivity = [
            ("bonds", 2),
            ("angles", 3),
            ("dihedrals", 4),
            ("impropers", 4),
        ]
        for name, n_cols in connectivity:
            df = getattr(self, name)
            if df is not None and not df.empty:
                indices = df[[f"atom_{i}" for i in range(1, n_cols + 1)]].to_numpy() - 1
                universe.add_TopologyAttr(name, [tuple(row) for row in indices])
