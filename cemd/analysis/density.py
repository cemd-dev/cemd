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

import dask
import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy import integrate
from tqdm import tqdm

_AXIS_MAP = {
    "x": {"axid": 0, "axida": 1, "axidb": 2},
    "y": {"axid": 1, "axida": 0, "axidb": 2},
    "z": {"axid": 2, "axida": 0, "axidb": 1},
}


def _get_axis_ids(axis: str) -> dict:
    if axis not in _AXIS_MAP:
        raise ValueError(f"axis must be 'x', 'y' or 'z', got '{axis}'")
    return _AXIS_MAP[axis]


def density_profile(
    universe: mda.Universe,
    atom_types: list[str | int],
    axis: str = "z",
    start: int = 0,
    end: int = -1,
    bin_size: float = 0.1,
) -> pd.DataFrame:
    """Create a DataFrame with the average density of atoms of given types along a given axis.

    Parameters
    ----------
    universe : mda.Universe
        The input MDAnalysis Universe to analyze.
    atom_types : list
        Atom types to compute the density for.
    axis : str
        Axis along which to calculate the density ('x', 'y', or 'z').
    start : int
        Starting frame index.
    end : int
        Ending frame index.
    bin_size : float
        Size of bins for the histogram (in Angstroms).

    Returns
    -------
    pd.DataFrame
        Average density profile for the specified atom types.
    """

    box = universe.dimensions

    ids = _get_axis_ids(axis)
    axid = ids["axid"]
    axida, axidb = ids["axida"], ids["axidb"]
    slice_vol = bin_size * box[axida] * box[axidb]

    bins = np.arange(0, box[axid], bin_size)
    pos = (bins[1:] + bins[:-1]) / 2
    density_total = []
    columns = []

    def count_pframe(frame_index, sel):

        sel.universe.trajectory[frame_index]

        posi = sel.positions[:, axid]

        posi = posi % box[axid]

        count = np.histogram(posi, bins=bins, range=[0, box[axid]])[0]

        return count

    if atom_types == "all":
        atom_types = np.unique(universe.atoms.types)

    for t in atom_types:
        print(f"Compute 1D atomic density of {t} atoms...")

        sel = universe.select_atoms(f"type {t}")

        # nframes = len(universe.trajectory[start:end])

        # job_list = []
        # for frame_index in tqdm( range(nframes) ):
        #     job_list.append(dask.delayed(count_pframe)(frame_index, sel))

        frames = (
            range(len(universe.trajectory))[start:end]
            if end != -1
            else range(len(universe.trajectory))[start:]
        )
        nframes = len(frames)

        if nframes == 0:
            raise ValueError(
                "Le slice de la trajectoire [start:end] ne contient aucune frame."
            )

        job_list = []
        for frame_index in tqdm(frames):
            job_list.append(dask.delayed(count_pframe)(frame_index, sel))

        result = dask.compute(job_list)
        atom_count = np.sum(result[0], axis=0)

        density = atom_count / slice_vol / nframes * 1000

        density_total.append(density)

        columns.append(f"{t}")

    return pd.DataFrame(np.array(density_total).T, columns=columns, index=pos)


def density_map(
    univ: mda.Universe,
    atom_types: str | list[str | int],
    interface_coordinate: float,
    axis: str = "z",
    eps: float = 3.0,
    start: int = 0,
    end: int = -1,
    bin_size: float = 0.1,
) -> pd.DataFrame:
    """Create a 2D density map of atoms within a specified distance of an interface.

    Parameters
    ----------
    univ : mda.Universe
        The input MDAnalysis Universe to analyze.
    atom_types : str or list
        Atom types to consider.
    interface_coordinate : float
        Interface coordinate along the given axis.
    axis : str
        Axis parallel to which to calculate the density.
    eps : float
        Half-thickness of the slab: atoms between ``interface_coordinate -
        eps`` and ``interface_coordinate + eps`` are counted.
    start : int
        Starting trajectory frame.
    end : int
        Ending trajectory frame.
    bin_size : float
        Size of bins for the histogram in Angstroms.

    Returns
    -------
    pd.DataFrame
        2D average density map.
    """

    box = univ.dimensions

    type_str = " ".join(atom_types) if isinstance(atom_types, list) else atom_types

    # A slab of half-thickness `eps` centred on the interface, which is what
    # "within a distance of an interface" means and what this function is
    # for -- the adsorbed layer, not everything on one side of the plane.
    # The selection used to be `prop axis < interface + eps`, i.e. the whole
    # half-cell below: on a solid-liquid system with the solid underneath,
    # that swept the entire solid into a map meant to show the liquid.
    lower = interface_coordinate - eps
    upper = interface_coordinate + eps
    sel = univ.select_atoms(
        f"type {type_str} and prop {axis} > {lower} and prop {axis} < {upper}",
        updating=True,
    )

    ids = _get_axis_ids(axis)
    axid = ids["axid"]
    axida, axidb = ids["axida"], ids["axidb"]
    bins_a = np.arange(0, box[axida], bin_size)
    bins_b = np.arange(0, box[axidb], bin_size)

    # A bin of the map is a column through the selected slab: bin_size by
    # bin_size in the plane, and as deep as the selection reaches along the
    # perpendicular axis. The 1D slab volume used before
    # (bin_size * box[a] * box[b]) is the volume of a whole slice, so the
    # reported densities were low by box[a] * box[b] / (bin_size * depth) --
    # a factor of 40 on a 20 A box at bin_size 0.5, and one that changes
    # with the box and the binning rather than being a constant offset.
    # It went unnoticed because a density map is read for its contrast.
    slab_depth = min(upper, float(box[axid])) - max(lower, 0.0)
    column_vol = bin_size * bin_size * slab_depth

    nframes = len(univ.trajectory[start:end])

    print(f"Compute 2D atomic density of {type_str} atoms...")

    pos_a_list, pos_b_list = [], []
    for ts in tqdm(univ.trajectory[start:end]):
        posi, posj = sel.positions[:, axida], sel.positions[:, axidb]

        posi = sel.positions[:, axida] % box[axida]
        posj = sel.positions[:, axidb] % box[axidb]

        pos_a_list.append(posi)
        pos_b_list.append(posj)

    pos_a = np.concatenate(pos_a_list)
    pos_b = np.concatenate(pos_b_list)

    hist, edges_a, edges_b = np.histogram2d(pos_a, pos_b, bins=(bins_a, bins_b))

    ra = (edges_a[1:] + edges_a[:-1]) / 2
    rb = (edges_b[1:] + edges_b[:-1]) / 2

    density = hist / column_vol / nframes * 1000

    return pd.DataFrame(density, index=ra, columns=rb)


def find_interfaces_coordinates(
    input_df: pd.DataFrame, solid_types: list[str | int], liquid_types: list[str | int]
) -> tuple[float, float, float, float]:
    """Calculate the solid/liquid interface coordinates for an interfacial system.

    Parameters
    ----------
    input_df : pd.DataFrame
        Input DataFrame of 1D atomic density profiles.
    solid_types : list
        List of atom types composing the solid phase.
    liquid_types : list
        List of atom types composing the liquid phase.

    Returns
    -------
    tuple
        Left liquid interface, left solid interface, right solid interface,
        and right liquid interface coordinates.
    """

    solid_series = pd.Series(input_df[solid_types].max(axis=1), index=input_df.index)
    liquid_series = pd.Series(input_df[liquid_types].max(axis=1), index=input_df.index)

    # On the left
    solid_left = (solid_series == 0).idxmax()
    solution_left = (liquid_series != 0).idxmax()

    # On the right
    df_solide_r = solid_series[::-1]
    df_solution = liquid_series.loc[solution_left:]
    solid_right = (df_solide_r == 0).idxmax()
    solution_right = (df_solution == 0).idxmax()

    return solution_left, solid_left, solid_right, solution_right


def shift_profile(
    input_df: pd.DataFrame, shift: float, csv_output: str = "new_density_profile.csv"
):
    """Shift the coordinates of a 1D density profile.

    Parameters
    ----------
    input_df : pd.DataFrame
        Input DataFrame of 1D atomic density profiles.
    shift : float
        Positive or negative shift in Angstroms.
    csv_output : str
        Output CSV file path.

    Returns
    -------
    pd.DataFrame
        Shifted density profile data.
    """

    # to avoid pylint warning
    # pylint: disable=E1101
    columns = input_df.columns

    r = input_df.index.values
    dr = r[1] - r[0]
    densities = input_df.values

    idx = int(shift / dr)
    new_densities = np.roll(densities, idx, axis=0)

    new_data = np.hstack((r.reshape(r.size, 1), new_densities))

    output_df = pd.DataFrame(new_data, columns=columns)

    if csv_output is not None:
        output_df.to_csv(csv_output)

    return output_df


def electrostatic_potential(
    input_df: pd.DataFrame, list_charges: list[float]
) -> tuple[pd.Series, pd.Series, pd.Series]:
    """Calculate charge distribution, electric field, and electrostatic potential.

    Parameters
    ----------
    input_df : pd.DataFrame
        DataFrame of atomic densities.
    list_charges : list
        List of charges for each atom type in the CSV file.

    Returns
    -------
    charge_density : pd.Series
        Calculated charge density profile.
    efield : pd.Series
        Calculated electric field profile.
    potential : pd.Series
        Calculated electrostatic potential profile.
    """

    def shift_constant(array: np.ndarray, threshold: float = 0.1) -> np.ndarray:
        """Shifts the array by subtracting the mean of its constant (flat) regions.

        This function automatically detects stationary or flat regions within a 1D signal by analyzing relative consecutive differences. It computes the average intensity of the signal across these stable baselines and subtracts it from the entire array, effectively centering the baseline/background around zero.

        Parameters
        ----------
            array
            The 1D input array containing the signal to be corrected.
            threshold
            The maximum relative difference between consecutive elements to consider a region as 'constant' or flat. Evaluated after scaling the array by its maximum value.

        Returns:
            np.ndarray: The baseline-corrected array where the flat regions
                are centered around zero.

        Example:
            >>> signal = np.array([10.0, 10.1, 10.0, 25.0, 50.0, 10.2, 10.0])
            >>> shift_constant(signal, threshold=0.02)
            array([-0.025,  0.075, -0.025, 14.975, 39.975,  0.175, -0.025])
        """

        narr = array / array.max()
        # Calculate the absolute differences between consecutive elements
        differences = np.abs(np.diff(narr))

        # Identify indices where the difference is within the threshold
        constant_indices = np.where(differences <= threshold)[0]

        # Adjust indices to include the next element in the array
        constant_indices = np.append(constant_indices, constant_indices + 1)

        # Get unique indices and sort them
        constant_indices = np.unique(constant_indices)

        result_array = array - np.mean(array[constant_indices])

        return result_array

    charge = np.array(list_charges)

    r = input_df.index
    dr = (r[1] - r[0]) * 1e-10

    charge_density = np.sum(input_df.values * charge.T, axis=1)
    charge_density_series = pd.Series(charge_density, index=r)

    efield = (
        integrate.cumulative_trapezoid(charge_density * 1.602e-19 * 1e27, dx=dr)
        / 8.854e-12
    )
    r = (r[1:] + r[:-1]) / 2
    efield_series = pd.Series(shift_constant(efield), index=r)

    potential = -integrate.cumulative_trapezoid(efield, dx=dr)
    r = (r[1:] + r[:-1]) / 2
    potential_series = pd.Series(shift_constant(potential), index=r)

    return charge_density_series, efield_series, potential_series
