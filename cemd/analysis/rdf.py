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
from tqdm import tqdm

from ..core.atomic_system import AtomicSystem

# __all__ = ['compute_rdf']


@singledispatch
def compute_rdf(
    source, type1: str, type2: str, cutoff: float = 10.0, dr: float = 0.1, **kwargs
) -> tuple[pd.DataFrame, float]:
    """Calculate the Radial Distribution Function (RDF) g(r), and its derivatives G(r) and n(r).

    This function uses singledispatch to handle either a static AtomicSystem or
    an MDAnalysis Universe (for trajectory analysis).

    Parameters
    ----------
    source : AtomicSystem or mda.Universe
        The source structure or trajectory to analyze.
    type1 : str
        Atom type for the first selection ('all' or specific type).
    type2 : str
        Atom type for the second selection ('all' or specific type).
    cutoff : float, optional
        Maximum distance for the RDF calculation (in Angstroms).
    dr : float, optional
        Bin width for the histogram (in Angstroms).
    skip : int, optional
        Number of frames to skip during trajectory analysis (only for mda.Universe).

    Returns
    -------
    pd.DataFrame
        Indexed by distance ``r``, with columns ``g_r``, ``G_r`` and ``n_r``.
    float
        The number density of the second selection, in atoms per cubic
        Angstrom, that ``g_r`` was normalised by. Returned because the
        normalisation is what makes an RDF comparable between systems, and
        checking it is the quickest way to tell a real first peak from one
        produced by the wrong particle count.

    Raises
    ------
    TypeError
        If the provided source type is not supported.
    """
    raise TypeError(f"Le type {type(source)} n'est pas supporté par compute_rdf.")


@compute_rdf.register(AtomicSystem)
def _(
    source: AtomicSystem, type1: str, type2: str, cutoff: float = 10.0, dr: float = 0.1
) -> tuple[pd.DataFrame, float]:
    # The AtomicSystem knows how to convert to a single frame MDA Universe
    u = source.to_mda()
    return compute_rdf(u, type1=type1, type2=type2, cutoff=cutoff, dr=dr)


@compute_rdf.register(mda.Universe)
def _(
    source: mda.Universe,
    type1: str,
    type2: str,
    cutoff: float = 10.0,
    dr: float = 0.1,
    skip: int = 1,
) -> tuple[pd.DataFrame, float]:
    bins = np.arange(0, cutoff + dr, dr)

    sel1 = source.select_atoms("all" if type1 == "all" else f"type {type1}")
    sel2 = source.select_atoms("all" if type2 == "all" else f"type {type2}")

    total_hist = np.zeros(len(bins) - 1, dtype=np.int64)
    nframes = 0
    volumes = []

    for ts in tqdm(source.trajectory[::skip], desc="Calcul de la RDF"):
        nframes += 1
        volumes.append(ts.volume)

        _, dists = capped_distance(
            sel1.positions, sel2.positions, max_cutoff=cutoff, box=ts.dimensions
        )
        dists = dists[dists > 0.01]

        hist, edges = np.histogram(dists, bins=bins)
        total_hist += hist

    # Normalization and physical processing (identical for one or a thousand frames)
    mean_volume = np.mean(volumes)
    mean_hist = total_hist / nframes

    r = (edges[:-1] + edges[1:]) / 2
    rho_target = len(sel2) / mean_volume
    shell_vols = 4 * np.pi * r**2 * dr
    g_r = mean_hist / (len(sel1) * rho_target * shell_vols)

    rho0 = source.atoms.n_atoms / mean_volume
    G_r = 4 * np.pi * r * rho0 * (g_r - 1)
    n_r = np.cumsum(4 * np.pi * r**2 * rho_target * g_r * dr)

    df_res = pd.DataFrame({"g_r": g_r, "G_r": G_r, "n_r": n_r}, index=r)

    df_res.index.name = "r"

    return df_res, rho_target
