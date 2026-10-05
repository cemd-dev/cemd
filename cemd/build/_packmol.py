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

import shutil
import subprocess
import tempfile
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np

from .. import AtomicSystem
from .._constants import MASSES_DICT
from ..core._format import canonical_ff_type, get_ff_key
from ._structures import STRUCTURES_DIR


def write_mda_pdb(obj, path: str | Path) -> None:
    """Write an MDAnalysis Universe/AtomGroup to a scratch PDB, quietly.

    `AtomicSystem` only carries the topology attributes it needs, but
    MDAnalysis's PDB writer wants several more (altLocs, resnames, icodes,
    segids, chainIDs, resids, occupancies, tempfactors, record_types,
    formalcharges) and warns once per missing one. These are Packmol
    scratch files read straight back into an `AtomicSystem` afterwards, so
    the values are irrelevant -- silence the warning here instead of
    inventing attributes nothing uses, or leaving every caller script to
    filter it out for itself.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")
        obj.write(str(path))


@dataclass
class PackmolStructure:
    """Molecular structure definition for Packmol."""

    structure: object
    number: int

    # Structure constraints
    fixed: tuple[float, float, float, float, float, float] | None = None
    center: bool = False

    inside_cube: tuple[float, float, float, float] | None = None
    outside_cube: tuple[float, float, float, float] | None = None
    inside_box: tuple[float, float, float, float, float, float] | None = None
    outside_box: tuple[float, float, float, float, float, float] | None = None
    inside_sphere: tuple[float, float, float, float] | None = None
    outside_sphere: tuple[float, float, float, float] | None = None
    inside_ellipsoid: tuple[float, float, float, float, float, float, float] | None = (
        None
    )
    outside_ellipsoid: tuple[float, float, float, float, float, float, float] | None = (
        None
    )
    above_plane: tuple[float, float, float, float] | None = None
    below_plane: tuple[float, float, float, float] | None = None
    inside_cylinder: (
        tuple[float, float, float, float, float, float, float, float] | None
    ) = None
    outside_cylinder: (
        tuple[float, float, float, float, float, float, float, float] | None
    ) = None
    over_xygauss: tuple[float, float, float, float, float, float] | None = None
    below_xygauss: tuple[float, float, float, float, float, float] | None = None

    # Structure options
    radius: float | None = None
    resnumbers: Literal[0, 1, 2, 3] | None = None
    chain: str | None = None
    changechains: bool = False
    segid: str | None = None
    connect: bool | None = None

    restart_to: str | None = None
    restart_from: str | None = None

    maxmove: int | None = None
    nloop: int | None = None

    # Rotation constraints
    constrain_rotation: dict[
        Literal["x", "y", "z"],
        tuple[float, float],
    ] = field(default_factory=dict)

    # Penalty function
    fscale: float | None = None
    use_short_tol: bool = False
    short_tol_dist: float | None = None
    short_tol_scale: float | None = None

    # Escape hatch for Packmol directives with no dedicated field above
    # (e.g. the "atoms N ... end atoms" sub-block to fix specific atoms).
    # Each entry is written verbatim as one line before "end structure".
    extra_instructions: list[str] = field(default_factory=list)

    def to_input(self, structure_path: str) -> str:
        """Convert the structure to Packmol syntax."""
        output = f"structure {structure_path}\n"
        output += f"  number {self.number}\n"

        if self.center:
            output += "  center\n"

        if self.fixed is not None:
            output += f"  fixed {_format_values(self.fixed)}\n"

        constraints = {
            "inside cube": self.inside_cube,
            "outside cube": self.outside_cube,
            "inside box": self.inside_box,
            "outside box": self.outside_box,
            "inside sphere": self.inside_sphere,
            "outside sphere": self.outside_sphere,
            "inside ellipsoid": self.inside_ellipsoid,
            "outside ellipsoid": self.outside_ellipsoid,
            "above plane": self.above_plane,
            "below plane": self.below_plane,
            "inside cylinder": self.inside_cylinder,
            "outside cylinder": self.outside_cylinder,
            "over xygauss": self.over_xygauss,
            "below xygauss": self.below_xygauss,
        }

        for keyword, values in constraints.items():
            if values is not None:
                output += f"  {keyword} {_format_values(values)}\n"

        if self.radius is not None:
            output += f"  radius {self.radius}\n"

        if self.resnumbers is not None:
            output += f"  resnumbers {self.resnumbers}\n"

        if self.chain is not None:
            output += f"  chain {self.chain}\n"

        if self.changechains:
            output += "  changechains\n"

        if self.segid is not None:
            output += f"  segid {self.segid}\n"

        if self.connect is not None:
            output += f"  connect {'yes' if self.connect else 'no'}\n"

        if self.restart_to is not None:
            output += f"  restart_to {self.restart_to}\n"

        if self.restart_from is not None:
            output += f"  restart_from {self.restart_from}\n"

        if self.maxmove is not None:
            output += f"  maxmove {self.maxmove}\n"

        if self.nloop is not None:
            output += f"  nloop {self.nloop}\n"

        for axis, (angle, tolerance) in self.constrain_rotation.items():
            output += f"  constrain_rotation {axis} {angle} {tolerance}\n"

        if self.fscale is not None:
            output += f"  fscale {self.fscale}\n"

        if self.use_short_tol:
            output += "  use_short_tol\n"

        if self.short_tol_dist is not None:
            output += f"  short_tol_dist {self.short_tol_dist}\n"

        if self.short_tol_scale is not None:
            output += f"  short_tol_scale {self.short_tol_scale}\n"

        for instruction in self.extra_instructions:
            output += f"  {instruction}\n"

        output += "end structure\n"

        return output


@dataclass
class PackmolInput:
    """Complete Packmol input."""

    # Mandatory/global options
    tolerance: float = 2.0
    output: str = "output.pdb"
    filetype: Literal["pdb", "tinker", "xyz"] = "pdb"

    # Periodic boundary conditions
    pbc: tuple[float, ...] | None = None

    # Global optimization options
    discale: float | None = None
    maxit: int | None = None
    movebadrandom: bool = False
    movefrac: float | None = None
    disable_movebad: bool = False

    # PDB / output options
    ignore_conect: bool = False
    non_standard_conect: bool = False
    writecrd: str | None = None
    add_amber_ter: bool = False
    amber_ter_preserve: bool = False
    add_box_sides: bool = False

    # System size / randomization
    sidemax: float | None = None
    seed: int | None = None
    randominitialpoint: bool = False
    avoid_overlap: bool | None = None

    # Precision / output
    precision: float | None = None
    writeout: int | None = None
    writebad: bool = False
    iprint1: int | None = None
    iprint2: int | None = None

    # Technical
    fbins: float | None = None
    hexadecimal_indices: bool = False
    chkgrad: bool = False

    # Restart
    restart_to: str | None = None
    restart_from: str | None = None

    # Packing structures
    structures: list[PackmolStructure] = field(default_factory=list)

    def to_input(self, structure_paths: list[str]) -> str:
        """Generate the Packmol input file."""
        if len(structure_paths) != len(self.structures):
            raise ValueError(
                "The number of structure paths must match the number of structures."
            )

        output = f"tolerance {self.tolerance}\n"
        output += f"output {self.output}\n"
        output += f"filetype {self.filetype}\n"

        if self.pbc is not None:
            output += f"pbc {_format_values(self.pbc)}\n"

        if self.discale is not None:
            output += f"discale {self.discale}\n"

        if self.maxit is not None:
            output += f"maxit {self.maxit}\n"

        if self.movebadrandom:
            output += "movebadrandom\n"

        if self.movefrac is not None:
            output += f"movefrac {self.movefrac}\n"

        if self.disable_movebad:
            output += "disable_movebad\n"

        if self.ignore_conect:
            output += "ignore_conect\n"

        if self.non_standard_conect:
            output += "non_standard_conect\n"

        if self.writecrd is not None:
            output += f"writecrd {self.writecrd}\n"

        if self.add_amber_ter:
            output += "add_amber_ter\n"

        if self.amber_ter_preserve:
            output += "amber_ter_preserve\n"

        if self.add_box_sides:
            output += "add_box_sides\n"

        if self.sidemax is not None:
            output += f"sidemax {self.sidemax}\n"

        if self.seed is not None:
            output += f"seed {self.seed}\n"

        if self.randominitialpoint:
            output += "randominitialpoint\n"

        if self.avoid_overlap is not None:
            output += f"avoid_overlap {'yes' if self.avoid_overlap else 'no'}\n"

        if self.precision is not None:
            output += f"precision {self.precision}\n"

        if self.writeout is not None:
            output += f"writeout {self.writeout}\n"

        if self.writebad:
            output += "writebad\n"

        if self.iprint1 is not None:
            output += f"iprint1 {self.iprint1}\n"

        if self.iprint2 is not None:
            output += f"iprint2 {self.iprint2}\n"

        if self.fbins is not None:
            output += f"fbins {self.fbins}\n"

        if self.hexadecimal_indices:
            output += "hexadecimal_indices\n"

        if self.chkgrad:
            output += "chkgrad\n"

        if self.restart_to is not None:
            output += f"restart_to {self.restart_to}\n"

        if self.restart_from is not None:
            output += f"restart_from {self.restart_from}\n"

        output += "\n"

        for structure, path in zip(
            self.structures,
            structure_paths,
            strict=True,
        ):
            output += structure.to_input(path)
            output += "\n"

        return output


def _format_values(values: tuple[float, ...]) -> str:
    """Format numerical values for Packmol."""
    return " ".join(str(value) for value in values)


def _rebuild_topology_from_templates(
    solution: AtomicSystem,
    structures: list[PackmolStructure],
    structure_paths: list[str],
) -> AtomicSystem:
    """Rebuild topology and restore force-field assignments."""

    u = solution.to_mda()

    topology_attrs = {
        "bonds": 2,
        "angles": 3,
        "dihedrals": 4,
        "impropers": 4,
    }

    dict_names = {
        "bonds": "bond",
        "angles": "angle",
        "dihedrals": "dihedral",
        "impropers": "improper",
    }

    new_topology = {attr: [] for attr in topology_attrs}
    atom_ff_keys = {}
    atom_charges = {}
    interaction_ff_keys = {kind: {} for kind in dict_names.values()}
    atom_type_mapping = {}
    used_global_atom_types = set()
    packmol_type_mapping = {}
    atom_offset = 0

    for struct_idx, (structure, structure_path) in enumerate(
        zip(structures, structure_paths, strict=True)
    ):
        if isinstance(structure.structure, AtomicSystem):
            template = structure.structure
        else:
            template = AtomicSystem.from_file(structure_path)

        n_copies = structure.number
        n_atoms = template.num_atoms
        template_atom_types = template.atoms["type"]

        # ==============================================================
        # Build local atom type -> global atom type mapping
        # ==============================================================

        local_to_global = {}

        # Iterate over every atom type actually present in the template,
        # not just the ones with an assigned ff_key: a structure loaded
        # straight from a file (e.g. a CIF) may have no ff_keys at all,
        # and its atom types still need a global-type mapping below.
        for local_type in template_atom_types.unique():
            local_type = str(local_type)
            ff_key = template.ff_keys.atom.get(local_type)
            identity = (local_type, ff_key)

            if identity in atom_type_mapping:
                global_type = atom_type_mapping[identity]
            elif local_type not in used_global_atom_types:
                global_type = local_type
                atom_type_mapping[identity] = global_type
                used_global_atom_types.add(global_type)
                if ff_key is not None:
                    atom_ff_keys[global_type] = ff_key
            else:
                index = 2
                while f"{local_type}_{index}" in used_global_atom_types:
                    index += 1
                global_type = f"{local_type}_{index}"
                atom_type_mapping[identity] = global_type
                used_global_atom_types.add(global_type)
                if ff_key is not None:
                    atom_ff_keys[global_type] = ff_key

            local_to_global[local_type] = global_type

            # A PDB round-trip (Packmol's own output included) can't carry
            # charge -- the standard format has no such field, and cemd's
            # own writer/reader pair used to silently corrupt it to 1.0
            # (see `PDBReader`). Carry the template's real per-type charge
            # through explicitly instead of trusting the reloaded PDB.
            atom_charges[global_type] = template.charges.get(local_type, 0.0)

        # ==============================================================
        # Build Packmol type -> global type mapping
        # ==============================================================

        new_types = np.array(u.atoms.types, dtype=object)
        for copy_idx in range(n_copies):
            offset = atom_offset + copy_idx * n_atoms
            for local_index, local_type in enumerate(template_atom_types):
                global_type = local_to_global[str(local_type)]
                tmp = f"__packmol_{offset + local_index}"
                new_types[offset + local_index] = tmp
                packmol_type_mapping[tmp] = global_type
        u.atoms.types = new_types

        # ==============================================================
        # Rebuild topology
        # ==============================================================

        for attr, n_per_interaction in topology_attrs.items():
            interactions = getattr(template, attr, None)

            if interactions is None or interactions.empty:
                continue

            kind = dict_names[attr]
            columns = [f"atom_{i}" for i in range(1, n_per_interaction + 1)]
            template_ff_keys = getattr(template.ff_keys, kind)

            for atom_indices in interactions[columns].to_numpy(dtype=int) - 1:
                local_atom_types = tuple(
                    str(template_atom_types.iloc[i]) for i in atom_indices
                )

                # Resolve FF key
                ff_key = get_ff_key(template_ff_keys, local_atom_types, kind)

                # Convert local types -> global types
                global_atom_types = tuple(
                    local_to_global[atom_type] for atom_type in local_atom_types
                )

                # For dihedrals, keep the original order
                if kind == "dihedral":
                    global_interaction_type = "-".join(global_atom_types)
                else:
                    global_interaction_type = canonical_ff_type(global_atom_types, kind)

                # A template can carry real connectivity while having no
                # force-field assignment yet (a molecule read from an SDF or
                # a PDB, say). That is chemistry, not bookkeeping: keep the
                # interaction and leave the ff key to be assigned later,
                # instead of silently dropping the bond.
                if ff_key is not None:
                    existing_ff_key = interaction_ff_keys[kind].get(
                        global_interaction_type
                    )
                    if existing_ff_key is not None and existing_ff_key != ff_key:
                        raise ValueError(
                            f"{kind.capitalize()} type {global_interaction_type!r} "
                            f"has multiple force-field keys: {existing_ff_key!r} "
                            f"and {ff_key!r}"
                        )

                    interaction_ff_keys[kind][global_interaction_type] = ff_key

                # Add interaction for every copy
                for copy_idx in range(n_copies):
                    offset = atom_offset + copy_idx * n_atoms
                    new_indices = tuple(index + offset for index in atom_indices)
                    new_topology[attr].append(new_indices)

        atom_offset += n_copies * n_atoms

    # ==================================================================
    # Apply topology to MDAnalysis
    # ==================================================================

    for attr, interactions in new_topology.items():
        if interactions:
            u.add_TopologyAttr(attr, interactions)

    # ==================================================================
    # Convert back to AtomicSystem
    # ==================================================================

    result = AtomicSystem.from_mda(u)

    # ==================================================================
    # Restore atom FF keys
    # ==================================================================

    result._ff_keys.atom.update(atom_ff_keys)

    # ==================================================================
    # Restore interaction FF keys
    # ==================================================================

    for kind, keys in interaction_ff_keys.items():
        target = getattr(result._ff_keys, kind)
        target.update(keys)

    # ==================================================================
    # Set atom types (after restoring FF keys)
    # ==================================================================

    result.set_types(packmol_type_mapping)

    # ==================================================================
    # Restore per-type charges from the original templates (see the
    # `atom_charges` comment above for why this can't be trusted to a
    # PDB round-trip); `atom_charges` is already keyed by global type,
    # matching what `set_types` just renamed atoms to.
    # ==================================================================

    result.set_charges(atom_charges)

    # ==================================================================
    # Restore force-field parameters
    # ==================================================================

    for structure in structures:
        if not isinstance(structure.structure, AtomicSystem):
            continue
        template = structure.structure
        result._ff_params.bond.update(template._ff_params.bond)
        result._ff_params.angle.update(template._ff_params.angle)
        result._ff_params.dihedral.update(template._ff_params.dihedral)
        result._ff_params.improper.update(template._ff_params.improper)
        result._ff_params.pair.update(template._ff_params.pair)

    return result


def get_structure_path(
    name: str,
    temp_dir: str | Path,
) -> Path:
    """Return the path to a structure file as a Path object.

    The function first looks for an LT file, then a PDB file, then an SDF file
    in the CEMD structure directory. If neither exists and the
    species is monoatomic, a temporary PDB file is generated.
    """
    name_lower = name.lower()
    temp_path = Path(temp_dir)

    # Check for a .lt file
    lt_path = STRUCTURES_DIR / f"{name_lower}.lt"
    if lt_path.exists():
        return lt_path

    # Check for a .pdb file
    pdb_path = STRUCTURES_DIR / f"{name_lower}.pdb"
    if pdb_path.exists():
        return pdb_path

    # Check for a .sdf file
    sdf_path = STRUCTURES_DIR / f"{name_lower}.sdf"
    if sdf_path.exists():
        return sdf_path

    # Generate a temporary file if the atom is monoatomic. MASSES_DICT is
    # keyed by proper element symbols (e.g. "Ca"), so normalize the case
    # regardless of how the caller spelled it (e.g. "ca", "CA").
    element = name.capitalize()
    if element in MASSES_DICT:
        temp_pdb = temp_path / f"{name_lower}.pdb"

        temp_pdb.write_text(
            f"HETATM    1 {element:>2s}  {element:>3s} A   1"
            "       0.000   0.000   0.000"
            "  1.00  0.00\n"
            "END\n",
            encoding="ascii",
        )
        return temp_pdb

    # Raise an error if no structure is found
    raise FileNotFoundError(
        f"Structure for '{name}' not found as ATB (MOLTEMPLATE), PDB or SDF and cannot be generated."
    )


def run_packmol(
    packmol: PackmolInput,
    output_path: str | Path | None = None,
) -> AtomicSystem | Path | str:
    """
    Run Packmol and return the resulting system.

    Parameters
    ----------
    packmol
        Packmol input definition.

    output_path
        If provided, copy the PDB generated by Packmol to this path
        and return the path. Otherwise, return an AtomicSystem with
        its topology reconstructed from the input templates.
    """
    with tempfile.TemporaryDirectory(dir=".") as tmp:
        tmp_path = Path(tmp)
        pdb_output_file = tmp_path / "tmp_out.pdb"
        packmol_input_file = tmp_path / "tmp.inp"

        structure_paths = []

        for index, structure in enumerate(packmol.structures):
            template = structure.structure

            if isinstance(template, AtomicSystem):
                structure_path = tmp_path / f"structure_{index}.pdb"

                template.write(str(structure_path))

            elif Path(template).is_file():
                structure_path = Path(template)

            else:
                structure_path = get_structure_path(
                    template,
                    tmp_path,
                )

            # Packmol only reads PDB here, but the bundled library also
            # holds SDF and moltemplate files (e.g. `ho.sdf`), which it
            # rejects with "Could not read any atom from file". Species
            # named in a blueprint rather than passed as an AtomicSystem
            # went through unconverted, so "HO", usable via
            # `structures={...}`, failed when used by name alone.
            if structure_path.suffix.lower() != ".pdb":
                converted = tmp_path / f"structure_{index}.pdb"
                AtomicSystem.from_file(str(structure_path)).write(str(converted))
                structure_path = converted

            structure_paths.append(structure_path)

        input_content = packmol.to_input(
            [str(p) for p in structure_paths],
        )

        input_content = input_content.replace(
            f"output {packmol.output}",
            f"output {pdb_output_file}",
            1,
        )

        packmol_input_file.write_text(input_content, encoding="utf-8")

        result = subprocess.run(
            f"packmol < {packmol_input_file}",
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # Packmol reports a dense-but-usable packing ("ended without perfect
        # packing") with a non-zero status, while still writing the best
        # solution it found -- which it describes as a reasonable starting
        # configuration. Tight interlayers hit this routinely, so judge the
        # run by whether an output was produced, not by the status alone.
        if not pdb_output_file.exists():
            raise RuntimeError(
                f"Packmol failed with return code {result.returncode} "
                "and produced no output."
            )

        if result.returncode != 0:
            warnings.warn(
                f"Packmol ended without a perfect packing (return code "
                f"{result.returncode}); using the best solution it found. "
                "Relax the structure before production, or give the "
                "molecules more room.",
                UserWarning,
                stacklevel=2,
            )

        if output_path is not None:
            out_path = Path(output_path)
            shutil.copy(pdb_output_file, out_path)
            return output_path

        solution = AtomicSystem.from_file(
            pdb_output_file,
        )

        return _rebuild_topology_from_templates(
            solution,
            packmol.structures,
            structure_paths,
        )
