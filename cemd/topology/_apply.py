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
# from __future__ import annotations

from __future__ import annotations

import re
from itertools import combinations
from typing import TYPE_CHECKING

import numpy as np

from .rules import DihedralRule, TopologyRule

if TYPE_CHECKING:
    import MDAnalysis as mda


# Updating the MDAnalysis topology
def _update_topo(univ, attr, new_entries):
    if not new_entries:
        return
    existing = list(getattr(univ, attr).indices) if hasattr(univ, attr) else []
    combined = list({tuple(e) for e in existing} | {tuple(e) for e in new_entries})
    if hasattr(univ, attr):
        univ.del_TopologyAttr(attr)
    univ.add_TopologyAttr(attr, combined)


#: Selection keywords whose result for one atom depends on other atoms.
#: A neighbor selection using one of them goes through the per-center
#: path, which evaluates it against the whole universe every time.
_NON_LOCAL = re.compile(
    r"\b(around|sphzone|sphlayer|cyzone|cylayer|isolayer|point|bonded|byres|"
    r"same|global|group|fullgroup)\b"
)


def apply_single_rule_to_universe(
    universe: mda.Universe, rule: TopologyRule
) -> mda.Universe:
    """
    Applies a single geometry-based connectivity rule in-place on a MDAnalysis Universe.

    Parameters
    ----------
    universe : mda.Universe
        The MDAnalysis universe to modify.
    rule : TopologyRule
        The topology rule to apply.

    Returns
    -------
    mda.Universe
        The modified universe.
    """
    if any(_NON_LOCAL.search(c.selection) for c in rule.neighbors):
        return _apply_rule_per_center(universe, rule)
    return _apply_rule_with_neighbor_lists(universe, rule)


def _type_names(selection: str) -> set[str] | None:
    """The types of a plain "type A B ..." selection, else None."""
    match = re.fullmatch(r"\s*type\s+(.+?)\s*", selection)
    if match is None:
        return None
    names = match.group(1).split()
    if any(n in ("and", "or", "not") or set(n) & set("*?[]()") for n in names):
        return None
    return set(names)


def _periodic_box(universe: mda.Universe):
    """The box `around` would use: the dimensions, if they are all set."""
    box = universe.dimensions
    if box is None or not np.all(box > 0):
        return None
    return box


def _apply_rule_with_neighbor_lists(
    universe: mda.Universe, rule: TopologyRule
) -> mda.Universe:
    """The per-center rule, with every neighbor search done up front.

    Selecting "around <cutoff> index <c>" for each center searched the
    whole system once per center: 70 s for the ClayFF rules on 24 000
    atoms, four times longer at each doubling. Here one `capped_distance`
    per criterion -- the search `around` itself runs, same cutoff, same
    periodic box -- lists every center's neighbors at once.

    Rules rename atoms as they go, and a later center sees the new names:
    a hydrogen already renamed "Hw" no longer counts as "type H". Each
    criterion therefore keeps a membership mask over all atoms, updated
    for the atoms each match renames. That is exact for selections where
    an atom's membership depends on that atom alone, which is why the
    non-local keywords go through `_apply_rule_per_center` instead.
    """
    from MDAnalysis.lib.distances import capped_distance

    center_atoms = universe.select_atoms(rule.center)
    atoms = universe.atoms
    centers = center_atoms.indices
    box = _periodic_box(universe)

    criteria = []
    for criterion in rule.neighbors:
        member = np.zeros(len(atoms), dtype=bool)
        member[universe.select_atoms(criterion.selection).indices] = True

        if len(centers):
            pairs = capped_distance(
                center_atoms.positions,
                atoms.positions,
                max_cutoff=criterion.cutoff,
                box=box,
                return_distances=False,
            )
        else:
            pairs = np.empty((0, 2), dtype=np.int64)
        # Grouped by center, each group in index order as a selection is.
        pairs = pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))]
        bounds = np.searchsorted(pairs[:, 0], np.arange(len(centers) + 1))
        criteria.append(
            (criterion, member, pairs[:, 1], bounds, _type_names(criterion.selection))
        )

    # Writing a type through MDAnalysis rewrites the whole type array: one
    # rename per match made the loop quadratic again (11 s of 15 on 60 000
    # atoms). When every criterion only tests types, the renames go into a
    # local array, written back once. A selection on anything else has to
    # see the universe up to date, so it keeps renaming in place.
    batch = all(type_names is not None for *_, type_names in criteria)
    types = atoms.types.copy() if batch else None
    renamed_any = False

    new_bonds = []
    new_angles = []
    new_impropers = []

    for row, c in enumerate(centers.tolist()):
        matched_neighbors = []
        neighbor_type_tasks = []
        is_valid = True

        for criterion, member, neighbors, bounds, _ in criteria:
            candidates = neighbors[bounds[row] : bounds[row + 1]]
            found = candidates[(candidates != c) & member[candidates]]

            if criterion.exact_match:
                if len(found) != criterion.count:
                    is_valid = False
                    break
            else:
                if len(found) < criterion.count:
                    is_valid = False
                    break

            matched_neighbors.extend(found.tolist())

            if criterion.new_type is not None:
                neighbor_type_tasks.append((found, criterion.new_type))

        if not is_valid:
            continue

        renamed = []
        if rule.new_type is not None:
            if batch:
                types[c] = rule.new_type
            else:
                atoms[c].type = rule.new_type
            renamed.append(c)
        for found, new_type in neighbor_type_tasks:
            if batch:
                types[found] = new_type
            else:
                atoms[found].types = new_type
            renamed.extend(found.tolist())

        if renamed:
            renamed_any = True
            renamed = np.unique(renamed)
            if batch:
                current_types = types[renamed]
            else:
                current = atoms[renamed]
                current_types = current.types
            for criterion, member, _, _, type_names in criteria:
                if type_names is not None:
                    member[renamed] = np.isin(current_types, list(type_names))
                else:
                    member[renamed] = False
                    member[current.select_atoms(criterion.selection).indices] = True

        if rule.bonds:
            for n in matched_neighbors:
                new_bonds.append(tuple(sorted((c, n))))

        if rule.angles and len(matched_neighbors) >= 2:
            for n1, n2 in combinations(matched_neighbors, 2):
                new_angles.append((n1, c, n2))

        if rule.impropers and len(matched_neighbors) >= 3:
            for n1, n2, n3 in combinations(matched_neighbors, 3):
                new_impropers.append((c, n1, n2, n3))

    if batch and renamed_any:
        atoms.types = types

    _update_topo(universe, "bonds", new_bonds)
    _update_topo(universe, "angles", new_angles)
    _update_topo(universe, "impropers", new_impropers)

    return universe


def _apply_rule_per_center(universe: mda.Universe, rule: TopologyRule) -> mda.Universe:
    """The rule evaluated with one selection per center: slow, but exact
    for any selection string."""
    center_atoms = universe.select_atoms(rule.center)

    new_bonds = []
    new_angles = []
    new_impropers = []

    for c in center_atoms:
        matched_neighbors = []
        neighbor_type_tasks = []
        is_valid = True

        for n_criterion in rule.neighbors:
            found = universe.select_atoms(
                f"({n_criterion.selection}) and "
                f"(around {n_criterion.cutoff} index {c.index}) and "
                f"not index {c.index}"
            )

            if n_criterion.exact_match:
                if len(found) != n_criterion.count:
                    is_valid = False
                    break
            else:
                if len(found) < n_criterion.count:
                    is_valid = False
                    break

            matched_neighbors.extend(list(found))

            if n_criterion.new_type is not None:
                neighbor_type_tasks.append((found, n_criterion.new_type))

        if not is_valid:
            continue

        # Apply new type to center
        if rule.new_type is not None:
            c.type = rule.new_type

        # Apply new types to neighbors
        for atoms, new_type in neighbor_type_tasks:
            atoms.types = new_type

        # Create bonds
        if rule.bonds:
            for n in matched_neighbors:
                new_bonds.append(tuple(sorted((c.index, n.index))))

        # Create angles
        if rule.angles and len(matched_neighbors) >= 2:
            for n1, n2 in combinations(matched_neighbors, 2):
                new_angles.append((n1.index, c.index, n2.index))

        # Create impropers
        if rule.impropers and len(matched_neighbors) >= 3:
            for n1, n2, n3 in combinations(matched_neighbors, 3):
                new_impropers.append((c.index, n1.index, n2.index, n3.index))

    _update_topo(universe, "bonds", new_bonds)
    _update_topo(universe, "angles", new_angles)
    _update_topo(universe, "impropers", new_impropers)

    return universe


def apply_single_dihedral_rule_to_universe(
    universe: mda.Universe,
    rule: DihedralRule,
    default_cutoff: float = 2.0,
) -> mda.Universe:
    """Generate proper dihedrals for a single DihedralRule, in-place.

    Parameters
    ----------
    universe : mda.Universe
        MDAnalysis universe containing the system.
    rule : DihedralRule
        The dihedral rule defining the i-j-k-l selection pattern.
    default_cutoff : float
        Default cutoff if not specified in rule.

    Returns
    -------
    mda.Universe
        The modified universe.
    """

    from MDAnalysis.lib.distances import capped_distance

    dihedrals = []
    seen_dihedrals = set()

    sel_i = universe.select_atoms(rule.i)
    sel_j = universe.select_atoms(rule.j)
    sel_k = universe.select_atoms(rule.k)
    sel_l = universe.select_atoms(rule.l_)

    cutoffs = rule.cutoffs
    if len(cutoffs) < 3:
        cutoffs = list(cutoffs) + [default_cutoff] * (3 - len(cutoffs))

    if len(sel_i) == 0 or len(sel_j) == 0 or len(sel_k) == 0 or len(sel_l) == 0:
        return universe

    pairs_ij, _ = capped_distance(
        sel_i.positions,
        sel_j.positions,
        max_cutoff=cutoffs[0],
        box=universe.dimensions,
    )

    pairs_jk, _ = capped_distance(
        sel_j.positions,
        sel_k.positions,
        max_cutoff=cutoffs[1],
        box=universe.dimensions,
    )

    pairs_kl, _ = capped_distance(
        sel_k.positions,
        sel_l.positions,
        max_cutoff=cutoffs[2],
        box=universe.dimensions,
    )

    idx_ij = {(sel_i.indices[i], sel_j.indices[j]) for i, j in pairs_ij}
    idx_jk = {(sel_j.indices[j], sel_k.indices[k]) for j, k in pairs_jk}
    idx_kl = {(sel_k.indices[k], sel_l.indices[l_]) for k, l_ in pairs_kl}

    # Indexed by the shared atom: matching every i-j pair against every
    # j-k pair and every k-l pair was cubic in the number of pairs.
    k_of_j = {}
    for j, k in idx_jk:
        k_of_j.setdefault(j, []).append(k)
    l_of_k = {}
    for k, l_ in idx_kl:
        l_of_k.setdefault(k, []).append(l_)

    for i, j in idx_ij:
        for k in k_of_j.get(j, ()):
            for l_ in l_of_k.get(k, ()):
                if i == l_:
                    continue

                dih = (i, j, k, l_)

                dih_rev = (l_, k, j, i)
                if dih not in seen_dihedrals and dih_rev not in seen_dihedrals:
                    seen_dihedrals.add(dih)
                    dihedrals.append(dih)

    _update_topo(universe, "dihedrals", dihedrals)

    return universe


def _generate_dihedrals_from_bonds(universe: mda.Universe) -> list[tuple]:
    """Generate all proper dihedrals from existing bonds.

    A dihedral i-j-k-l is defined by a central bond j-k,
    where i is a neighbor of j and l is a neighbor of k.
    """
    if not hasattr(universe, "bonds") or len(universe.bonds) == 0:
        return []

    neighbors = {i: set() for i in range(len(universe.atoms))}
    for bond in universe.bonds.indices:
        i, j = bond
        neighbors[i].add(j)
        neighbors[j].add(i)

    dihedrals = set()

    for j, k in universe.bonds.indices:
        for i in neighbors[j] - {k}:
            for l_ in neighbors[k] - {j}:
                if i == l_:
                    continue
                dih = (i, j, k, l_)
                dih_rev = (l_, k, j, i)
                if dih_rev not in dihedrals:
                    dihedrals.add(dih)

    return list(dihedrals)


# def apply_clayff_rules(universe: mda.Universe) -> tuple[mda.Universe, dict]:
#     """Applies ClayFF to the universe."""
#     for rule in CLAYFF_RULES:
#         universe = apply_single_rule_to_universe(universe, rule)
#     return universe, {}


# def apply_cshff_rules(universe: mda.Universe) -> tuple[mda.Universe, dict]:
#     """Applies CSHFF to the universe and returns the calcium indices to modify."""
#     from ..build.cement_hydrates._silicate_helpers import get_interlayer_ca_indices

#     universe, _ = apply_clayff_rules(universe)
#     list_ids_cw = get_interlayer_ca_indices(universe)
#     actions = {"rename_atoms": (list_ids_cw, "Cw")}
#     return universe, actions
