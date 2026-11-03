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

from functools import singledispatch

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

from ..core.atomic_system import AtomicSystem
from .silicates import TYPES_PRESET


def _coordination_numbers(center, neighbor, cutoff, box) -> np.ndarray:
    """Number of `neighbor` atoms within `cutoff` of each `center` atom, in this frame."""
    pairs = capped_distance(
        center.positions, neighbor.positions, cutoff, box=box, return_distances=False
    )
    # An atom selected in both groups would count itself
    same = center.indices[pairs[:, 0]] == neighbor.indices[pairs[:, 1]]
    return np.bincount(pairs[~same, 0], minlength=len(center))


def _process_coordination(
    universe, center_types, neighbor_types, cutoff, start, stop, step, max_cn
):
    center = universe.select_atoms(f"type {center_types}")
    neighbor = universe.select_atoms(f"type {neighbor_types}")
    if len(center) == 0:
        raise ValueError(f"No atom of type(s) {center_types!r} in the system.")
    if len(neighbor) == 0:
        raise ValueError(f"No atom of type(s) {neighbor_types!r} in the system.")

    rows, frames = [], []
    for ts in universe.trajectory[start:stop:step]:
        cn = _coordination_numbers(center, neighbor, cutoff, ts.dimensions)
        counts = np.bincount(np.minimum(cn, max_cn), minlength=max_cn + 1)
        rows.append(np.append(100.0 * counts / len(center), cn.mean()))
        frames.append(ts.frame)

    # The last column gathers every coordination of max_cn and above
    columns = [*range(max_cn), f"{max_cn}+", "mean"]
    return pd.DataFrame(rows, index=pd.Index(frames, name="frame"), columns=columns)


@singledispatch
def analyze_coordination(
    source,
    center_types: str = "Al",
    neighbor_types: str = TYPES_PRESET["o_types"],
    cutoff: float = 2.4,
    start: int = 0,
    stop: int | None = None,
    step: int = 1,
    max_cn: int = 8,
) -> pd.DataFrame:
    """Distribution of the coordination number of a cation, frame by frame.

    The fractions of 4-, 5- and 6-fold aluminium are what 27Al MAS NMR
    measures in a calcined or ground clay: Al(VI) in the octahedral sheet of
    the raw mineral, Al(IV) and Al(V) as it is disordered. Counting them in
    a simulation gives a number that compares directly.

    Parameters
    ----------
    source : AtomicSystem or mda.Universe
        The structure, or the trajectory, to analyze.
    center_types : str, optional
        Space-separated atom types whose coordination is counted.
    neighbor_types : str, optional
        Space-separated atom types counted as neighbours. Default: the
        oxygen types of the bundled force fields, plus plain ``O``.
    cutoff : float, optional
        Neighbour distance cutoff in Angstroms. Take it at the first minimum
        of the center-neighbour RDF (`compute_rdf`); 2.4 A suits Al-O and
        Si-O in aluminosilicates.
    start, stop, step : int, optional
        Frames analyzed (only for mda.Universe). ``stop=None`` runs to the
        last frame.
    max_cn : int, optional
        Coordination numbers of ``max_cn`` and above are gathered in the
        last column.

    Returns
    -------
    pandas.DataFrame
        One row per frame (indexed by frame number; a single row 0 for an
        AtomicSystem). Columns ``0`` to ``max_cn - 1`` and ``"{max_cn}+"``
        hold the percentage of center atoms with that coordination number;
        ``"mean"`` holds the mean coordination number.
    """
    raise TypeError(
        f"analyze_coordination() does not support source of type {type(source).__name__}; "
        "expected an AtomicSystem or an MDAnalysis Universe."
    )


# The defaults are repeated on each implementation: `singledispatch` calls it
# directly, without going through the generic function above (see the note
# in silicates.py).
@analyze_coordination.register(AtomicSystem)
def _(
    source: AtomicSystem,
    center_types: str = "Al",
    neighbor_types: str = TYPES_PRESET["o_types"],
    cutoff: float = 2.4,
    start: int = 0,
    stop: int | None = None,
    step: int = 1,
    max_cn: int = 8,
) -> pd.DataFrame:
    return _process_coordination(
        source.to_mda(), center_types, neighbor_types, cutoff, 0, None, 1, max_cn
    )


@analyze_coordination.register(mda.Universe)
def _(
    source: mda.Universe,
    center_types: str = "Al",
    neighbor_types: str = TYPES_PRESET["o_types"],
    cutoff: float = 2.4,
    start: int = 0,
    stop: int | None = None,
    step: int = 1,
    max_cn: int = 8,
) -> pd.DataFrame:
    return _process_coordination(
        source, center_types, neighbor_types, cutoff, start, stop, step, max_cn
    )
