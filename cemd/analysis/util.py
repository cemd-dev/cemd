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

import os
from collections.abc import Sequence

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import minimize_vectors, transform_RtoS, transform_StoR
from tqdm import tqdm


def require_box(source, caller: str) -> np.ndarray:
    """Return the periodic box of a universe or timestep, or say it is missing.

    MDAnalysis reports ``dimensions`` as ``None`` for a trajectory that
    carries no cell (some readers report zeros instead). Analyses that wrap
    positions or divide by the cell volume cannot proceed without it, and
    would otherwise fail on ``TypeError: 'NoneType' object is not
    subscriptable`` -- or, for the RDF, quietly return ``inf``.

    Parameters
    ----------
    source : mda.Universe or Timestep
        Anything exposing ``dimensions``.
    caller : str
        Name of the analysis, quoted in the error message.

    Returns
    -------
    np.ndarray
        ``[lx, ly, lz, alpha, beta, gamma]``.

    Raises
    ------
    ValueError
        If the box is missing, non-finite or has a non-positive length.
    """
    dimensions = source.dimensions
    if (
        dimensions is None
        or not np.all(np.isfinite(dimensions))
        or np.any(np.asarray(dimensions)[:3] <= 0)
    ):
        raise ValueError(
            f"{caller} needs the periodic box, but this trajectory carries "
            "none. Set it with `universe.dimensions = [lx, ly, lz, 90, 90, 90]` "
            "(for every frame of an in-memory trajectory) or load a file that "
            "stores the cell."
        )
    return np.asarray(dimensions)


def write_dcd(
    universe: mda.Universe,
    output_path: str,
    selection: str = "all",
    start: int = 0,
    end: int = None,
    step: int = 1,
) -> None:
    """Write trajectory frames from an MDAnalysis Universe to a new DCD file.

    Parameters
    ----------
    universe : mda.Universe
        MDAnalysis Universe object containing topology and trajectory data.
    output_path : str
        Path where the new DCD file will be saved.
    selection : str, optional
        Selection string (MDAnalysis syntax) specifying which atoms to write.
    start : int, optional
        Starting frame index (0-indexed).
    end : int, optional
        Ending frame index. If None, writes up to the last available frame.
    step : int, optional
        Step size (stride) for writing frames.
    """
    # check on file extension
    _, ext = os.path.splitext(output_path)
    if ext.lower() != ".dcd":
        raise TypeError("The output file extension must be '.dcd'.")

    # apply the atom selection
    selected_atoms = universe.select_atoms(selection)
    if len(selected_atoms) == 0:
        raise ValueError(f"Selection string '{selection}' matched 0 atoms.")

    print(f"Preparing to write DCD file: {output_path}")
    print(f"Selected {len(selected_atoms)} out of {len(universe.atoms)} atoms.")

    # handle the trajectory slice definition
    trajectory_slice = universe.trajectory[start:end:step]
    total_frames_to_write = len(trajectory_slice)

    print(
        f"Writing {total_frames_to_write} frames (start={start}, end={end}, step={step})..."
    )

    # open the writer and stream the frames
    # the Writer requires the output filename and the exact number of atoms being written
    with mda.Writer(output_path, selected_atoms.n_atoms) as W:
        for ts in trajectory_slice:
            W.write(selected_atoms)

    print("DCD trajectory written successfully!")


def _whole_center_of_mass(group: mda.AtomGroup, box: np.ndarray) -> np.ndarray:
    """Centre of mass of a group made whole across the periodic boundaries.

    Along each cell vector, the group is cut through the widest band it
    leaves empty -- the pore, for a slab -- rather than at the box edge, so
    a group straddling the boundary is averaged in one piece.
    """
    frac = transform_RtoS(group.positions, box) % 1.0
    for axis in range(3):
        values = frac[:, axis]
        ordered = np.sort(values)
        gaps = np.diff(ordered)
        # The band across the boundary, from the last atom back to the first.
        wrap_gap = ordered[0] + 1.0 - ordered[-1]
        if len(gaps) and gaps.max() > wrap_gap:
            cut = ordered[gaps.argmax()]
            values[values <= cut] += 1.0
    positions = transform_StoR(frac.astype(np.float32), box)
    return np.average(positions, axis=0, weights=group.masses)


def _com_trajectory(group: mda.AtomGroup, unwrap: bool, caller: str):
    """Yield the centre of mass of `group` at every frame of its trajectory.

    With `unwrap`, each centre is the periodic image closest to the
    previous one, so that it never jumps by a box length between frames.
    """
    previous = None
    for ts in tqdm(group.universe.trajectory):
        if not unwrap:
            yield group.center_of_mass()
            continue
        box = require_box(ts, caller)
        current = _whole_center_of_mass(group, box)
        if previous is not None:
            current = previous + minimize_vectors(
                (current - previous)[None].astype(np.float32), box
            )[0]
        previous = current
        yield current


def shift2com(
    universe: mda.Universe,
    atom_types: list[str | int],
    output_trajectory: str = "recentered_traj.dcd",
    unwrap: bool = False,
) -> None:
    """Recenter all atoms relative to the center of mass (COM) of a selection.

    Parameters
    ----------
    universe : mda.Universe
        MDAnalysis universe object.
    atom_types : list
        List of atom types used to define the reference center of mass.
    output_trajectory : str, optional
        Path to the output recentered DCD trajectory file.
    unwrap : bool, default=False
        Make the reference group whole across the periodic boundaries
        before taking its COM. Set it for a compact group the boundary may
        cut through, such as a slab: without it, a slab split by the
        boundary has its COM in the middle of the pore. Leave it off for a
        group spread through the whole box (ions in solution), whose COM is
        then the plain average of the wrapped positions.
    """

    selection_string = f"type {' '.join(atom_types)}"
    ref_atoms = universe.select_atoms(selection_string)

    if len(ref_atoms) == 0:
        raise ValueError(f"No atoms found for the type(s): {atom_types}")

    all_atoms = universe.atoms

    with mda.Writer(output_trajectory, all_atoms.n_atoms) as W:
        for com in _com_trajectory(ref_atoms, unwrap, "shift2com"):
            all_atoms.positions -= com

            W.write(all_atoms)


# def shift2com(itraj, psf, atypes, otraj='recentered_traj.dcd'):
#     """Recenter all atoms w.r.t to the COM of a selection based on atom types

#     Parameters
#     ----------
#         itraj: str
#             An input DCD trajectory
#         psf: str
#             A PSF file
#         atypes: list of str
#             List of the atom types in the reference selection
#         otraj: str
#             An output DCD trajectory

#     """

#     sel_str = " ".join(map(str, atypes))

#     subprocess.run(['vmd', '-dispdev', 'text', '-e', shift2com_tcl, '-args', itraj, psf, sel_str, otraj], stdout=subprocess.DEVNULL, check=True)


def minmax_position(
    universe: mda.Universe,
    atom_types: list[str | int],
    axis: str = "z",
    bounds: Sequence[float] = None,
    start: int = 0,
    end: int = -1,
) -> tuple[float, float]:
    """Calculate the mean of the minimum and maximum coordinates along an axis.

    Parameters
    ----------
    universe : mda.Universe
        MDAnalysis Universe object containing topology and trajectory.
    atom_types : list
        List of atom type strings to select.
    axis : str, optional
        Axis along which to calculate the coordinates ('x', 'y', or 'z').
    bounds : Sequence[float], optional
        Optional spatial limits [min, max] in Angstroms to filter atoms.
    start : int, optional
        Starting frame index for the trajectory slice.
    end : int, optional
        Ending frame index for the trajectory slice.

    Returns
    -------
    tuple
        A tuple containing (mean_min, mean_max) coordinates.
    """

    typestr = " ".join(map(str, atom_types))

    if bounds is None:
        selstr = f"type {typestr}"
        print(f"Compute the min and max coordinate along {axis} for {typestr} atoms...")
    else:
        selstr = f"type {typestr} and prop {axis} >= {bounds[0]} and prop {axis} < {bounds[1]}"
        print(
            f"Compute the min and max coordinate along {axis} for {typestr} atoms between {axis}={bounds[0]} angströms and {axis}={bounds[1]}..."
        )

    sel = universe.select_atoms(selstr)

    axes = {"x": 0, "y": 1, "z": 2}
    if axis not in axes:
        raise ValueError(f"axis must be 'x', 'y' or 'z', not {axis!r}.")
    axid = axes[axis]

    mins, maxs = [], []

    for ts in tqdm(universe.trajectory[start:end]):
        mins.append(sel.positions[:, axid].min())
        maxs.append(sel.positions[:, axid].max())

    return np.mean(np.array(mins)), np.mean(np.array(maxs))


def mean_pos(universe) -> pd.DataFrame:
    """Return a DataFrame containing the mean position of each atom over the trajectory.

    Parameters
    ----------
    universe : mda.Universe
        MDAnalysis universe object.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns 'type', 'x', 'y', and 'z'.
    """

    ag = universe.atoms

    # Molecules are made whole frame by frame, without touching the
    # universe: an unwrap *transformation* stayed attached to it, and a
    # second call failed on "Can't add transformations again".
    mean_pos = np.zeros((len(ag), 3))
    for ts in tqdm(universe.trajectory):
        mean_pos += ag.unwrap(compound="fragments", inplace=False)

    mean_pos /= len(universe.trajectory)

    df = pd.DataFrame(mean_pos, columns=["x", "y", "z"])
    df.insert(0, "type", ag.types)

    return df


def com(
    universe: mda.Universe, atom_types: list[str | int], unwrap: bool = False
) -> np.ndarray:
    """Calculate the mean position of the center of mass of a selection.

    Parameters
    ----------
    universe : mda.Universe
        MDAnalysis universe object.
    atom_types : list
        List of atom type strings to select.
    unwrap : bool, default=False
        Make the group whole across the periodic boundaries first; see
        :func:`shift2com`. The result is then brought back into the box.

    Returns
    -------
    np.ndarray
        The mean position of the center of mass over the trajectory.
    """

    sel = universe.select_atoms("type {}".format(" ".join(atom_types)))

    mean = np.mean(list(_com_trajectory(sel, unwrap, "com")), axis=0)

    if unwrap:
        box = require_box(universe, "com")
        frac = transform_RtoS(mean[None].astype(np.float32), box) % 1.0
        mean = transform_StoR(frac, box)[0]

    return mean
