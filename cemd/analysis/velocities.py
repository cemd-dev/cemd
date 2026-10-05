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

import MDAnalysis as mda
import numpy as np
import pandas as pd
from tqdm import tqdm


def velocity_profile(
    universe: mda.Universe, atom_types: str | list[str], bin_size: float = 0.1
) -> pd.Series:
    """Mean speed of the selected atoms as a function of their z position.

    Parameters
    ----------
    universe : mda.Universe
        Universe whose trajectory carries velocities (a LAMMPS dump or a
        TRR, not a DCD).
    atom_types : str or list of str
        Atom type(s) to select.
    bin_size : float, default=0.1
        Width of the z bins, in Angstroms.

    Returns
    -------
    pd.Series
        Mean speed ``|v|`` in each bin, in the trajectory's velocity unit
        (Angstrom/ps for MDAnalysis readers), indexed by the bin centre.
        A bin no atom visited holds NaN.
    """
    if isinstance(atom_types, str):
        atom_types = [atom_types]
    sel = universe.select_atoms("type {}".format(" ".join(map(str, atom_types))))

    positions = []
    speeds = []

    for ts in tqdm(universe.trajectory):
        positions.append(sel.positions[:, 2])
        # One speed per atom, so that it lines up with its z position.
        speeds.append(np.linalg.norm(sel.velocities, axis=1))

    positions_array = np.concatenate(positions)
    speeds_array = np.concatenate(speeds)

    edges = np.arange(positions_array.min(), positions_array.max(), bin_size)
    centres = (edges[:-1] + edges[1:]) / 2

    counts, _ = np.histogram(positions_array, bins=edges)
    sums, _ = np.histogram(positions_array, bins=edges, weights=speeds_array)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_speed = sums / counts

    return pd.Series(mean_speed, centres)
