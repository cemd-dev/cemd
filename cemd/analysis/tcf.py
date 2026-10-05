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
from MDAnalysis.lib.distances import calc_bonds, capped_distance
from scipy.optimize import leastsq
from tqdm import tqdm

from .util import require_box

__all__ = ["bondcorr", "lifetime"]


def bondcorr(
    universe: mda.Universe,
    atom_types_a: str | list,
    atom_types_b: str | list,
    distance: float,
    dt: float = 100,
    nblocks: int = None,
    corrlength: int = 5000,
    gaplength: int = None,
) -> pd.DataFrame:
    """Calculate the time correlation function (TCF) of a bond.

    Tracks the survival probability of chemical bonds between two atom types
    over a block-averaged trajectory. A bond is considered broken as soon as
    its length exceeds the defined distance threshold.

    Parameters
    ----------
    universe : mda.Universe
        The input MDAnalysis Universe to analyze.
    atom_types_a : str or list
        The atom type(s) for the first element forming the bond.
    atom_types_b : str or list
        The atom type(s) for the second element forming the bond.
    distance : float
        The maximum cutoff distance to define a bond (in Angstroms).
    dt : float, optional
        The simulation timestep of the trajectory (in femtoseconds).
    nblocks : int, optional
        Number of blocks used for averaging. If None, as many as fit in the
        trajectory.
    corrlength : int, optional
        Correlation length (number of frames per block).
    gaplength : int, optional
        Number of frames between the start of consecutive blocks. Defaults
        to a fiftieth of the trajectory.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing "t [ps]" (correlation time) and "TCF"
        (survival probability).

    Raises
    ------
    ValueError
        If `corrlength`, `gaplength`, or `nblocks` are too large for
        the available trajectory length.
    """

    def _make_type_selection(univ, atom_types):
        if isinstance(atom_types, str):
            atom_types = [atom_types]
        return univ.select_atoms("type {}".format(" ".join(atom_types)))

    nframes = len(universe.trajectory)

    if gaplength is None:
        gaplength = max(1, nframes // 50)

    # Blocks start at 0, gaplength, 2*gaplength... and each needs corrlength
    # frames: the last one must start at or before nframes - corrlength.
    max_blocks = (nframes - corrlength) // gaplength + 1 if corrlength <= nframes else 0
    if nblocks is None:
        # One short of max_blocks when the last block would end exactly on
        # the last frame: the count this function has always used.
        nblocks = (nframes - corrlength) // gaplength

    if nblocks < 1 or nblocks > max_blocks:
        raise ValueError(
            f"{nblocks} block(s) of {corrlength} frames, {gaplength} frames "
            f"apart, do not fit in a trajectory of {nframes} frames."
        )

    sel_a = _make_type_selection(universe, atom_types_a)
    sel_b = _make_type_selection(universe, atom_types_b)

    starts = np.arange(0, nblocks * gaplength, gaplength)
    ends = starts + corrlength

    master_results = np.zeros_like(np.arange(starts[0], ends[0], 1), dtype=np.float32)
    nobs = 0

    for i, (start, end) in enumerate(tqdm(zip(starts, ends), total=len(starts)), 1):
        # Initialize the bond recording
        ts = universe.trajectory[start]

        pairs = capped_distance(
            sel_a.positions,
            sel_b.positions,
            max_cutoff=distance,
            box=require_box(ts, "bondcorr"),
            return_distances=False,
        )
        idx1, idx2 = pairs[:, 0], pairs[:, 1]

        # An atom of both selections (like-species pair, Ow-Ow) is found at
        # zero distance from itself: a "bond" that never breaks.
        distinct = sel_a.indices[idx1] != sel_b.indices[idx2]
        idx1, idx2 = idx1[distinct], idx2[distinct]

        nbonds = len(idx1)

        if nbonds != 0:
            results = np.zeros_like(np.arange(start, end, 1), dtype=np.float32)

            # Check if the initial bond still exist during the correlation time
            for j, ts in enumerate(universe.trajectory[start:end]):
                b = calc_bonds(
                    sel_a.positions[idx1],
                    sel_b.positions[idx2],
                    box=require_box(ts, "bondcorr"),
                )

                winners = b < distance
                results[j] = winners.sum()

                idx1 = idx1[winners]
                idx2 = idx2[winners]

                if len(idx1) == 0:  # Once everyone has lost, the fun stops
                    break

            if np.any(results):
                results /= nbonds
            master_results += results

            nobs += 1

    master_results /= nobs

    time = np.arange(starts[0], ends[0], 1) * dt / 1000  # time in ps
    odf = pd.DataFrame(np.c_[time, master_results], columns=["t [ps]", "TCF"])

    return odf


def lifetime(odf: pd.DataFrame, corrtime: float) -> tuple[float, np.ndarray]:
    """Calculate the bond lifetime by fitting the time correlation function.

    Fits the TCF with a double-exponential decay model using a bounded
    least-squares approach. The average lifetime is computed as the
    weighted sum of the characteristic decay times.

    Parameters
    ----------
    odf : pd.DataFrame
        Input DataFrame containing "t [ps]" and "TCF" columns.
    corrtime : float
        A characteristic time scale used for the initial optimization guess.

    Returns
    -------
    tau : float
        The calculated average bond lifetime.
    p : np.ndarray
        The optimized fit parameters: [A1, tau1, tau2].

    Notes
    -----
    The fitting function used is:
    $$f(x) = A_1 \\cdot \\exp(-x / \\tau_1) + (1 - A_1) \\cdot \\exp(-x / \\tau_2)$$
    """

    def double_exp(x, A1, tau1, tau2) -> np.ndarray:
        """Sum of two exponential functions"""
        A2 = 1 - A1
        return A1 * np.exp(-x / tau1) + A2 * np.exp(-x / tau2)

    def within_bounds(p) -> bool:
        """Returns True/False if boundary conditions are met or not.
        Uses length of p to detect whether it's handling continuous /
        intermittent

        Boundary conditions are:
            0 < A_x < 1
            sum(A_x) < 1
            0 < tau_x
        """
        A1, tau1, tau2 = p
        return 0.0 < A1 < 1.0 and tau1 > 0.0 and tau2 > 0.0
        # if len(p) == 3:
        #     A1, tau1, tau2 = p
        #     return (A1 > 0.0) & (A1 < 1.0) & \
        #             (tau1 > 0.0) & (tau2 > 0.0)
        # elif len(p) == 5:
        #     A1, A2, tau1, tau2, tau3 = p
        #     return (A1 > 0.0) & (A1 < 1.0) & (A2 > 0.0) & \
        #             (A2 < 1.0) & ((A1 + A2) < 1.0) & \
        #             (tau1 > 0.0) & (tau2 > 0.0) & (tau3 > 0.0)

    def err(p, x, y):
        """Custom residual function, returns real residual if all
        boundaries are met, else returns a large number to trick the
        leastsq algorithm
        """
        if within_bounds(p):
            return y - double_exp(x, *p)
        else:
            return np.full_like(y, 100000)

    p_guess = (0.5, 10 * corrtime, corrtime)

    t = odf["t [ps]"]
    tcf = odf["TCF"]

    p, _, _, _, _ = leastsq(err, p_guess, args=(t, tcf), full_output=True)

    A1, tau1, tau2 = p
    A2 = 1 - A1
    tau = A1 * tau1 + A2 * tau2

    return tau, p
