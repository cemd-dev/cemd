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

import random
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ...analysis.util import require_box
from .._packmol import write_mda_pdb

if TYPE_CHECKING:
    import MDAnalysis as mda


def grouped_average(
    array: np.ndarray, expected_size: int
) -> tuple[np.ndarray, list[np.ndarray], list[int]]:
    """
    Groups coordinates that are close to each other and returns their averages.
    Optimized to avoid large memory allocation and fix boolean ambiguity errors.
    """
    if array.size == 0:
        return np.array([]), [], []

    # Sort the array to make grouping much faster (O(N log N))
    sorted_arr = np.sort(array)
    tol = 0.2

    def get_groups(data, tolerance):
        # Find gaps between points larger than tolerance
        gaps = np.diff(data) > tolerance
        # Create group IDs (0, 0, 1, 1, 1, 2...)
        group_ids = np.concatenate(([0], np.cumsum(gaps)))
        return group_ids

    # Iterate tolerance until we reach the expected number of layers
    group_ids = get_groups(sorted_arr, tol)
    num_groups = group_ids[-1] + 1

    while num_groups > expected_size:
        tol += 0.2
        group_ids = get_groups(sorted_arr, tol)
        num_groups = group_ids[-1] + 1

    # Calculate averages and collect indices
    averaged_values = []
    indices_list = []
    counts = []

    for g_id in range(num_groups):
        # Get indices where the sorted array belongs to this group
        mask = group_ids == g_id
        group_data = sorted_arr[mask]

        averaged_values.append(np.mean(group_data))
        counts.append(len(group_data))

        # To get indices relative to the ORIGINAL (unsorted) array:
        # We find where original values match the values in this group
        original_indices = np.where(np.isin(array, group_data))[0]
        indices_list.append(original_indices)

    return np.array(averaged_values), indices_list, counts


def _detect_z_factor(universe: mda.Universe) -> int:
    """
    Detects the supercell multiplication factor along the z-axis.
    Assumes a unit cell height of ~22 Å (Tobermorite 11) or ~28 Å (Tobermorite 14).
    """
    # Get the z-dimension of the simulation box
    box_z = require_box(universe, "Silicate analysis")[2]

    # Determine the reference unit cell height based on the total box size
    # If the box is a multiple of ~28 Å, we use 28, otherwise we default to 22.
    if abs(box_z % 28.0) < abs(box_z % 22.0):
        unit_cell_z = 28.0
    else:
        unit_cell_z = 22.0

    # Calculate the multiplier factor (minimum of 1)
    multiplier = max(1, int(round(box_z / unit_cell_z)))
    return multiplier


def _get_bridging_silicates_indices(universe: mda.Universe) -> np.ndarray:
    """Returns indices of bridging silicates using _get_bridging_zplanes."""
    si_sel = universe.select_atoms("type Si")
    si_indices = si_sel.indices

    (_, si_indices_per_plane, si_count_per_plane) = _get_silicate_planes(universe)

    max_num = max(si_count_per_plane)
    bridging_groups = [
        g for g, n in zip(si_indices_per_plane, si_count_per_plane) if n != max_num
    ]
    return si_indices[np.concatenate(bridging_groups).flatten()]


def _get_silicate_planes(
    universe: mda.Universe,
) -> tuple[np.ndarray, list[np.ndarray], list[int]]:
    """Groups the Si by plane z and sorts them, keeping the indices
    origin of each plane (used by _find_symmetric_bridging_pairs)."""
    # Perform the selection internally to simplify the function call
    si_sel = universe.select_atoms("type Si")
    si_pos_z = si_sel.positions[:, 2]

    supercell_factor = _detect_z_factor(universe)

    (si_plane_positions, si_indices_per_plane, si_count_per_plane) = grouped_average(
        si_pos_z, 8 * supercell_factor
    )

    order = np.argsort(si_plane_positions)
    si_plane_positions = si_plane_positions[order]
    si_indices_per_plane = [si_indices_per_plane[i] for i in order]
    si_count_per_plane = [si_count_per_plane[i] for i in order]
    return si_plane_positions, si_indices_per_plane, si_count_per_plane


def find_dense_layer_seam_offset(universe: mda.Universe) -> float | None:
    """
    Find a z-offset that puts the periodic seam (z=0 / box_z) inside a
    dense, intra-sheet silicate gap instead of inside an interlayer.

    Packmol has no notion of periodicity: when a C-S-H interlayer gap
    straddles the box's z boundary (the last silicate plane wrapping
    around to the first one), water/Ca get packed by
    :func:`~._interlayer_helpers.fill_csh_interlayers` at raw
    z-coordinates beyond ``box_z``, without ever being checked against the
    real framework atoms sitting at the periodic-equivalent position near
    z=0 -- producing severe atom overlaps right at the seam. Moving the
    seam into solid material instead sidesteps the problem: every
    interlayer then lies fully inside ``[0, box_z)``, so the existing
    plane-bounded Packmol packing (already correct in x/y) just works.

    Returns
    -------
    float or None
        The offset to subtract from every z-coordinate (then wrap into
        ``[0, box_z)``), or ``None`` if fewer than two silicate planes are
        found, or none of the gaps between them is a dense (<=5 A)
        intra-sheet gap to move the seam into.
    """
    box_z = require_box(universe, "Seam alignment")[2]
    layers, _, _ = _get_silicate_planes(universe)
    nplanes = len(layers)
    if nplanes < 2:
        return None

    for i in range(nplanes):
        z_start = layers[i]
        z_end = layers[(i + 1) % nplanes]
        dist = (z_end - z_start) % box_z
        if dist <= 5.0:
            return (z_start + dist / 2.0) % box_z

    return None


def _find_symmetric_bridging_pairs(
    universe: mda.Universe, planar_distance_cutoff: float | None = 5.0
) -> list[list[tuple[int, int]]]:
    """
    Identifies bridging silicates and associates them in symmetrical pairs
    across each interlayer pore, grouped by pore.

    Parameters
    ----------
    universe : mda.Universe
        The MDAnalysis universe containing the structure.
    planar_distance_cutoff : float or None, default=5.0
        Maximum planar distance (in XY plane) for two Si atoms to be considered
        a symmetric pair. If None, planar distance is not used and atoms are
        paired in order of appearance.

    Returns
    -------
    list[list[tuple[int, int]]]
        List of pores, each containing a list of (idx1, idx2) pairs of bridging
        silicates that are symmetrically positioned.
    """
    box = require_box(universe, "Bridging-silicate pairing")
    box_x, box_y, box_z = box[0], box[1], box[2]

    si_sel = universe.select_atoms("type Si")
    si_indices = si_sel.indices
    si_pos = si_sel.positions

    (si_plane_positions, si_indices_per_plane, si_count_per_plane) = (
        _get_silicate_planes(universe)
    )

    max_num = max(si_count_per_plane)
    nplanes = len(si_plane_positions)

    def _pbc_dist_xy(p, q) -> float:
        dx = p[0] - q[0]
        dy = p[1] - q[1]
        dx -= box_x * round(dx / box_x)
        dy -= box_y * round(dy / box_y)
        return (dx * dx + dy * dy) ** 0.5

    bridging_plane_indices = [
        i for i in range(nplanes) if si_count_per_plane[i] != max_num
    ]

    n_bridging_planes = len(bridging_plane_indices)

    bridging_pairs_per_pore = []
    used = set()

    for idx_b in range(n_bridging_planes):
        i = bridging_plane_indices[idx_b]
        j = bridging_plane_indices[(idx_b + 1) % n_bridging_planes]

        gap = (si_plane_positions[j] - si_plane_positions[i]) % box_z

        if planar_distance_cutoff is not None and gap > planar_distance_cutoff:
            continue

        atoms_i = si_indices[np.array(si_indices_per_plane[i])]
        atoms_j = si_indices[np.array(si_indices_per_plane[j])]
        xy_i = si_pos[np.array(si_indices_per_plane[i])][:, :2]
        xy_j = si_pos[np.array(si_indices_per_plane[j])][:, :2]

        current_pore_pairs = []

        for k, idx1 in enumerate(atoms_i):
            if idx1 in used:
                continue

            best_l, min_dist = None, float("inf")

            for l_, idx2 in enumerate(atoms_j):
                if idx2 in used:
                    continue

                dist = _pbc_dist_xy(xy_i[k], xy_j[l_])

                if dist < min_dist:
                    min_dist, best_l = dist, l_

            if best_l is not None and min_dist < 5.0:
                idx2 = atoms_j[best_l]
                current_pore_pairs.append((idx1, idx2))
                used.add(idx1)
                used.add(idx2)

        if current_pore_pairs:
            bridging_pairs_per_pore.append(current_pore_pairs)

    return bridging_pairs_per_pore


def calculate_csh_modifiers(
    nsi: int, nca: int, target_cs_ratio: float, min_mcl: float
) -> tuple[int, int, float]:
    """Calculate the vacancy fraction and check that it is higher than 1 /(1 + MCL) to keep a minimum of bridging silicates. Calculate the number of bridging SiO2 and Ca2+ to add to reach the target C/S."""
    nsi_to_remove = nsi - round(nca / target_cs_ratio)
    vacancy_fraction = nsi_to_remove / nsi
    min_vacancy_fraction = 1 / (1 + min_mcl)

    nca_to_add = 0
    if vacancy_fraction > min_vacancy_fraction:
        while vacancy_fraction > min_vacancy_fraction:
            nsi_to_remove -= 1
            vacancy_fraction = nsi_to_remove / nsi

        # The product is rounded to 9 decimals before flooring: in binary,
        # 1.4 * 180 is 251.99999999999997, so a bare floor() dropped one
        # calcium and the built system came out at Ca/Si = 1.3944 instead
        # of the requested 1.4. The rounding only absorbs the
        # representation error -- a genuinely fractional target still
        # floors down, keeping the ratio from overshooting.
        n_ca_target = np.floor(round(target_cs_ratio * (nsi - nsi_to_remove), 9))
        nca_to_add = int(max(0, n_ca_target - nca))

    return nsi_to_remove, nca_to_add, vacancy_fraction


def remove_bridging_silicates(
    univ: mda.Universe, nsi_to_remove: int, symmetry: bool = True
) -> mda.Universe:
    if nsi_to_remove <= 0:
        return univ

    if symmetry:
        pores = _find_symmetric_bridging_pairs(univ)
        for pore in pores:
            np.random.shuffle(pore)
        n_pores = len(pores)
        to_remove = []
        count = 0
        pore_idx = 0

        seconds_to_process = [[] for _ in range(n_pores)]

        while count < nsi_to_remove:
            if not any(pores):
                break

            current_pore = pores[pore_idx % n_pores]

            if current_pore:
                item = current_pore.pop(0)

                if isinstance(item, (list, tuple)) and len(item) == 2:
                    idx1, idx2 = item
                    options = [idx1, idx2] if idx2 is not None else [idx1]
                    idx2remove = np.random.choice(options)
                    options.remove(idx2remove)
                    to_remove.append(idx2remove)
                    seconds_to_process[pore_idx % n_pores].append(options[0])
                else:
                    to_remove.append(item)

                count += 1
            pore_idx += 1

        pore_idx = 0
        while count < nsi_to_remove:
            if not any(seconds_to_process):
                break

            current_seconds = seconds_to_process[pore_idx % n_pores]
            if current_seconds:
                to_remove.append(current_seconds.pop(0))
                count += 1
            pore_idx += 1

    else:
        bridging_indices = _get_bridging_silicates_indices(univ)
        to_remove = np.random.choice(bridging_indices, nsi_to_remove, replace=False)

    # Update the universe
    idx_str = " ".join(map(str, to_remove))
    univ = univ.select_atoms(f"not index {idx_str}")
    univ = univ.select_atoms("not (name O and not around 2.1 type Si)")

    return univ


def substitute_si_by_al(
    univ: mda.Universe, nal: int, symmetry: bool = True
) -> tuple[mda.Universe, list[int]]:
    """
    Substitutes bridging Si with Al according to a ratio.

    Parameters
    ----------
    univ : mda.Universe
        The MDAnalysis Universe containing the structure.
    nal : int
        Number of Al atoms to substitute.
    symmetry : bool, optional
        If True, substitutes Si by Al in symmetric pairs across interlayers.
        If False, substitutes randomly. Default is True.
    """
    if nal <= 0:
        return univ, []

    to_sub = []

    if symmetry:
        pores = _find_symmetric_bridging_pairs(univ, planar_distance_cutoff=None)
        random.shuffle(pores)

        n_pores = len(pores)
        count = 0
        pore_idx = 0

        while count < nal:
            if not any(pores):
                break

            current_pore = pores[pore_idx % n_pores]

            if current_pore:
                idx1, idx2 = current_pore.pop(0)

                to_sub.append(idx1)
                count += 1

                if idx2 is not None and count < nal:
                    to_sub.append(idx2)
                    count += 1

            pore_idx += 1
    else:
        bridging_indices = _get_bridging_silicates_indices(univ)
        to_sub = list(np.random.choice(bridging_indices, nal, replace=False))

    for idx in to_sub:
        atom = univ.atoms[idx]
        atom.type = "Al"
        atom.name = "Al"

    substituted_ids = univ.atoms[to_sub].ids.tolist()

    return univ, substituted_ids


def neutralize_csh_charge(pdb_path: str | Path) -> float:
    """Remove hydrogens to neutralize the C-S-H system"""

    import MDAnalysis as mda

    univ = mda.Universe(pdb_path)

    si_sel = univ.select_atoms("name Si SI")
    ca_sel = univ.select_atoms("name Ca CA")
    # Water atoms come from h2o.lt and are named O1/H1/H2, so O and H are
    # selected by element; a "name H" selection silently matched nothing.
    o_sel = univ.select_atoms("element O")
    h_sel = univ.select_atoms("element H")
    total_charge = 2 * len(ca_sel) + len(h_sel) + 4 * len(si_sel) - 2 * len(o_sel)
    h_candidates = list(h_sel.indices[::2])
    if len(h_candidates) < total_charge:
        total_charge = len(h_candidates)
    to_remove = random.sample(h_candidates, total_charge)

    # Filter universe and overwrite PDB
    remaining_indices = np.setdiff1d(univ.atoms.indices, to_remove)
    write_mda_pdb(univ.atoms[remaining_indices], pdb_path)

    return total_charge / 2


# def get_interlayer_ca_indices(universe: mda.Universe) -> np.ndarray:
#     """Get IDs of interlayer calcium in a C-S-H structure."""
#     box = universe.dimensions

#     si_plane_positions, _, si_count_per_plane = _get_silicate_planes(universe)
#     max_num = max(si_count_per_plane)
#     n_planes = len(si_plane_positions)

#     si_pairing_plane_indices = [
#         i for i in range(n_planes) if si_count_per_plane[i] == max_num
#     ]
#     si_pairing_plane_positions = si_plane_positions[si_pairing_plane_indices]

#     interlayer_ca_indices = []

#     # Internal layers
#     for i in range(len(si_pairing_plane_positions) - 1):
#         gap = si_pairing_plane_positions[i + 1] - si_pairing_plane_positions[i]
#         if gap > 5:
#             sel = universe.select_atoms(
#                 f"(type Ca) and "
#                 f"(prop z > {si_pairing_plane_positions[i]}) and "
#                 f"(prop z < {si_pairing_plane_positions[i + 1]})"
#             )
#             interlayer_ca_indices.append(sel.ids)

#     # Periodic layer (box edge)
#     gap_pbc = abs(
#         si_pairing_plane_positions[-1] - si_pairing_plane_positions[0] - box[2]
#     )
#     if gap_pbc > 5:
#         sel = universe.select_atoms(
#             f"(type Ca) and "
#             f"((prop z > {si_pairing_plane_positions[-1]}) or "
#             f"(prop z < {si_pairing_plane_positions[0]}))"
#         )
#         interlayer_ca_indices.append(sel.ids)

#     return (
#         np.concatenate(interlayer_ca_indices) if interlayer_ca_indices else np.array([])
#     )
