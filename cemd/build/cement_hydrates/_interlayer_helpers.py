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

import MDAnalysis as mda
import numpy as np

from ...core._format import lattice2vectors
from ...core.atomic_system import AtomicSystem
from .._packmol import (
    PackmolInput,
    PackmolStructure,
    get_structure_path,
    run_packmol,
    write_mda_pdb,
)
from ._silicate_helpers import grouped_average

# Number density of bulk liquid water (1 g/cm^3, 18 g/mol), in molecules/A^3.
# Below this density, Packmol's default tolerance (2.0 A) packs quickly. Above
# it -- i.e. once an interlayer is asked to hold water denser than the bulk
# liquid, as happens with nanoconfined water at high ws_ratio -- the same
# tolerance describes a target that is barely (or no longer) geometrically
# achievable in a thin slab, and Packmol's GENCAN optimizer can spend minutes
# hunting for a packing that doesn't exist before giving up with the "ended
# without a perfect packing" warning.
_BULK_WATER_DENSITY = 0.03346  # molecules / A^3

# Never relax below this: still far above typical bond lengths, so Packmol's
# random start stays a safe input for the ReaxFF relaxation that follows.
_MIN_PACKING_TOLERANCE = 1.5


def _estimate_packing_tolerance(
    n_molecules: int,
    volume: float,
    max_tolerance: float = 2.0,
) -> float:
    """Pick a Packmol tolerance from the target packing density.

    Scales down linearly, past the bulk-water density, from
    ``max_tolerance`` at bulk density to ``_MIN_PACKING_TOLERANCE`` once the
    target is ~1.7x bulk density -- calibrated against CSH interlayer
    packings, where a small reduction (2.0 -> 1.8 A) turns a multi-minute,
    warning-prone packing into one that converges in seconds.
    """
    if n_molecules == 0 or volume <= 0:
        return max_tolerance

    density = n_molecules / volume
    excess = density / _BULK_WATER_DENSITY - 1.0

    if excess <= 0:
        return max_tolerance

    tolerance = max_tolerance - 1.0 * excess
    return max(_MIN_PACKING_TOLERANCE, min(max_tolerance, tolerance))


def _plane(vecp: np.ndarray, vecq: np.ndarray, vecr: np.ndarray) -> str:
    """Given three 3D-vectors, return the parameters a, b, c, d of the equation of the plane of coordinates (x, y, z) that goes through them, where:
    .. math:: ax + by + cz = d
    """

    vecu = vecq - vecp
    vecv = vecr - vecp
    vecw = np.cross(vecu, vecv)

    return f"{vecw[0]} {vecw[1]} {vecw[2]} {np.dot(vecw, vecp)}"


def _get_packmol_bounding_planes(box: np.ndarray) -> tuple:
    """Computes the 4 bounding planes for Packmol from box vectors."""
    vec0 = np.zeros(3)
    vec1, vec2, vec3 = lattice2vectors(box)
    vec1, vec2, vec3 = vec1 * 0.95, vec2 * 0.95, vec3 * 0.95

    return (
        _plane(vec0, vec2, vec3),  # Left
        _plane(vec1, vec1 + vec2, vec1 + vec3),  # Right
        _plane(vec0, vec3, vec1),  # Bottom
        _plane(vec2, vec2 + vec3, vec1 + vec2),  # Top
    )


def distribute_species_in_layers(
    num_water: int,
    nca_to_add: int,
    si_posz: np.ndarray,
    supercell_factor: int,
    box_z: float,
):  # -> tuple[list, list, list, int]:
    """
    Calculates layer positions and distributes H2O and Ca per layer.
    """
    # Find layers and keep only the main silicate planes (max atoms)
    si_layers_pos, _, list_num = grouped_average(si_posz, 8 * supercell_factor)
    si_layers_pos = si_layers_pos[np.array(list_num) == max(list_num)]

    nplanes = len(si_layers_pos)
    interlayers_bounds = []

    for i in range(nplanes):
        z_start = si_layers_pos[i]
        z_end = si_layers_pos[(i + 1) % nplanes]

        # Management of periodicity at the edges of the box (PBC)
        dist = (z_end - z_start) % box_z

        # An intersheet space is typically > 5.0 Å (while an intrasheet space is ~3.5 Å)
        if dist > 5.0:
            interlayers_bounds.append((z_start, z_end, dist))

    nlayers = len(interlayers_bounds)

    # Helper to distribute a total number into n parts as evenly as possible
    def _distribute_evenly(total, n):
        if n == 0:
            return []
        base, extra = divmod(total, n)
        return [base + (1 if i < extra else 0) for i in range(n)]

    nw_layers = _distribute_evenly(num_water, nlayers)
    nca_layers = _distribute_evenly(nca_to_add, nlayers)

    return interlayers_bounds, nw_layers, nca_layers, nlayers


def fill_csh_interlayers(
    univ: mda.Universe,
    box: np.ndarray,
    interlayers_bounds: list[tuple[float, float, float]],
    nw_layers: list,
    nca_layers: list,
    nca_to_add: int,
    nlayers: int,
    tmp: str | Path,
    progress_callback=None,
) -> Path:
    """Runs Packmol iteratively to fill each CSH interlayer with water and Ca."""
    tmp_path = Path(tmp)
    left, right, bottom, top = _get_packmol_bounding_planes(box)

    # get_structure_path may resolve "h2o" to a .lt (moltemplate) file,
    # which Packmol cannot read directly; load it as an AtomicSystem so
    # run_packmol rewrites it to PDB first.
    h2o = AtomicSystem.from_file(get_structure_path("h2o", tmp_path))
    ca_pdb = get_structure_path("ca", tmp_path) if nca_to_add != 0 else None

    current_pdb = tmp_path / "tmp0.pdb"
    write_mda_pdb(univ, current_pdb)

    for i in range(nlayers):
        z_start, z_end, dist = interlayers_bounds[i]
        z_upper = z_start + dist if z_end < z_start else z_end

        # Restrict the fill to this interlayer and to the box's transverse
        # bounds. Packmol only accepts one plane per "above/below plane"
        # line, so each simultaneous half-space constraint needs its own
        # instruction (PackmolStructure.above_plane/below_plane only model
        # a single plane each).
        plane_instructions = [
            f"above plane 0.0 0.0 1.0 {z_start + 1.5}",
            f"above plane {bottom}",
            f"above plane {left}",
            f"below plane 0.0 0.0 1.0 {z_upper - 1.5}",
            f"below plane {top}",
            f"below plane {right}",
        ]

        structures = [
            PackmolStructure(
                structure=current_pdb,
                number=1,
                inside_box=(0.0, 0.0, 0.0, float(box[0]), float(box[1]), float(box[2])),
                fixed=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            ),
            PackmolStructure(
                structure=h2o,
                number=int(nw_layers[i]),
                extra_instructions=list(plane_instructions),
            ),
        ]

        if nca_to_add != 0:
            structures.append(
                PackmolStructure(
                    structure=ca_pdb,
                    number=int(nca_layers[i]),
                    extra_instructions=list(plane_instructions),
                )
            )

        next_pdb = tmp_path / f"tmp{i + 1}.pdb"

        # Pick the tolerance from how dense this particular interlayer is
        # asked to be. Nanoconfined layers (high ws_ratio, thin gap) can
        # exceed bulk water density, at which point Packmol's default
        # tolerance describes a packing that barely exists (or doesn't),
        # and spends minutes searching before giving up -- see
        # `_estimate_packing_tolerance`.
        layer_volume = float(box[0]) * float(box[1]) * (dist - 3.0)
        n_layer_molecules = int(nw_layers[i]) + int(nca_layers[i])
        tolerance = _estimate_packing_tolerance(n_layer_molecules, layer_volume)

        # Construction of the Packmol configuration object
        packmol_input = PackmolInput(
            tolerance=tolerance,
            output=str(next_pdb),
            structures=structures,
        )

        run_packmol(packmol_input, next_pdb)
        current_pdb = next_pdb

        msg = f"Adding {nw_layers[i]} H2O" + (
            f" and {nca_layers[i]} Ca2+" if nca_to_add != 0 else ""
        )
        print(f"{msg} in layer {i + 1}")
        if progress_callback:
            p = 20 + int((i / nlayers) * 70)
            progress_callback(p, f"{msg} in interlayer {i + 1}/{nlayers}...")

    return current_pdb
