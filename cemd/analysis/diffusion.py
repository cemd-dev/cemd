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

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
import scipy.stats as stats
from tqdm import tqdm


def _linear_fit(
    t: np.ndarray, msd: np.ndarray, sigma_msd: np.ndarray = None
) -> tuple[float, float]:

    if sigma_msd is not None:
        w = 1 / sigma_msd

        fit, cov_matrix = np.polyfit(t, msd, 1, w=w, cov="unscaled")

    else:
        fit, cov_matrix = np.polyfit(t, msd, 1, cov=True)

    m, _ = fit

    # Extract standard errors from the scaled covariance matrix
    se_m_cov = np.sqrt(cov_matrix[0, 0])

    # Calculate the confidence intervals (e.g., 95% confidence)
    alpha = 0.05
    t_value = stats.t.ppf(1 - alpha / 2, df=len(t) - 2)
    ci_m = t_value * se_m_cov

    return m, ci_m


def msd(
    univ: mda.Universe,
    atom_type: str,
    dt: float = 100.0,
    nblocks: int = None,
    corrlength: int = None,
    gaplength: int = None,
) -> pd.DataFrame:
    """Calculate the diffusion coefficient in a bulk solution from the mean squared displacement (MSD).

    Uses Einstein's equation:
    .. math:: D = \\frac{1}{6} \\frac{\\langle [\\mathbf{r}(t_0 + \\tau) - \\mathbf{r}(t_0)]^2 \\rangle}{\\tau}

    Parameters
    ----------
    univ : MDAnalysis.Universe
        The input MDAnalysis Universe to analyze.
    atom_type : str
        The atom type to compute the diffusion coefficient.
    dt : float
        Time between two saved frames, in femtoseconds -- the MD timestep
        multiplied by the dump interval, **not** the MD timestep itself.

        .. warning::

           Do not take this from ``universe.trajectory.dt``. A DCD written
           by LAMMPS stores a time field that MDAnalysis reinterprets, and
           the result need not be the real interval: on the reference
           trajectory it reads 4.888821 ps where the frames are 100 fs
           apart, a factor of 49. The mean-squared displacement stays
           perfectly linear either way, so the diffusion coefficient comes
           out wrong by that factor with nothing to suggest anything is
           amiss. Take the value from the input deck that produced the
           trajectory.

    nblocks : int, optional
        Number of blocks.
    corrlength : int, optional
        Number of frames in a block.
    gaplength : int, optional
        Number of frames between blocks.

    Returns
    -------
    pd.DataFrame
        DataFrame containing MSD values over time.
    """

    sel = univ.select_atoms(f"type {atom_type}")
    box = univ.dimensions

    if nblocks is None or corrlength is None or gaplength is None:
        corrlength = int(len(univ.trajectory) / 2)
        gaplength = int(len(univ.trajectory) / 50)
        if gaplength == 0:
            gaplength = 1
            print(
                "Warning: The trajectory is too short so the gap between blocks was set to 1."
            )
        nblocks = int((len(univ.trajectory) - corrlength) / gaplength)

    if (gaplength * nblocks + corrlength) > len(univ.trajectory):
        raise ValueError(
            "Gap between correlation block, correlation length, or number of block too large."
        )

    print(
        f"Compute the profile of the mean-squared displacement on {len(sel)} {atom_type} atoms..."
    )
    print(f"Total number of frames: {len(univ.trajectory)}")
    print(f"Timestep between frames: {dt} fs")
    print(f"Number of blocks for calculation: {nblocks}")
    print(f"Length of blocks: {corrlength} frames / {corrlength * dt * 1e-3:.2f} ps")
    print(
        f"Length between origins: {gaplength} frames / {gaplength * dt * 1e-3:.2f} ps\n"
    )

    if nblocks > 1:
        starts = np.arange(0, nblocks * gaplength, gaplength)
        ends = starts + corrlength
    else:
        starts = np.array([0])
        ends = np.array([corrlength])

    t = np.arange(0, corrlength, 1) * dt / 1000  # time in picosecond

    msd_dic = {
        "xx": [],
        "yy": [],
        "zz": [],
        "xy": [],
        "xz": [],
        "yz": [],
    }

    for i, (start, end) in enumerate(tqdm(zip(starts, ends), total=len(starts)), 1):
        # initialize the trajectory
        univ.trajectory[start]

        # initialize positions w.r.t the center of mass to remove drift
        pos = sel.positions
        rr = pos - np.mean(pos, axis=0)

        # initialize total dr and msd
        d = np.zeros(pos.shape)
        msd_values = {key: [] for key in msd_dic.keys()}

        for ts in univ.trajectory[start:end]:
            # collect positions w.r.t the center of mass to remove the drift
            pos = sel.positions
            r = pos - np.mean(pos, axis=0)

            # calculate the differential displacement w.r.t last frame positions
            dd = r - rr

            # Taking PBC into account
            dd -= box[:3].T * (dd / box[:3].T).round()

            # Add the the differential displacement to the total displacement
            d += dd

            msd_values["xx"].append(np.mean(d[:, 0] * d[:, 0]))
            msd_values["yy"].append(np.mean(d[:, 1] * d[:, 1]))
            msd_values["zz"].append(np.mean(d[:, 2] * d[:, 2]))
            msd_values["xy"].append(np.mean(d[:, 0] * d[:, 1]))
            msd_values["xz"].append(np.mean(d[:, 0] * d[:, 2]))
            msd_values["yz"].append(np.mean(d[:, 1] * d[:, 2]))

            # save last frame position w.r.t the center of mass to remove the drift
            rr = r

        for key, values in msd_values.items():
            msd_dic[key].append(np.array(values))

    data_dic = {key: np.mean(values, axis=0) for key, values in msd_dic.items()}

    dfo = pd.DataFrame(data_dic, index=t)
    dfo.index.name = "time"

    return dfo


def msd_profile(
    univ: mda.Universe,
    atom_type: str,
    dt: int = 100,
    axis: str = "z",
    nblocks: int = None,
    corrlength: int = 100,
    gaplength: int = 100,
    delta: float = 1.0,
) -> tuple[pd.MultiIndex, pd.MultiIndex]:
    """Calculate a profile of diffusion coefficients at a liquid/solid interface.

    Parameters
    ----------
    univ : MDAnalysis.Universe
        The input MDAnalysis Universe to analyze.
    atom_type : str
        The atom type to compute the diffusion coefficient.
    dt : int
        Timestep between each frame of the trajectory (fs).
    axis : str
        Axis along which the profile will be calculated ('x', 'y', or 'z').
    nblocks : int, optional
        Number of blocks.
    corrlength : int, optional
        Number of frames in a block.
    gaplength : int, optional
        Number of frames between blocks.
    delta : float
        Size of the spatial binning (Angstrom).

    Returns
    -------
    pd.MultiIndex
        MSD values per bin.
    pd.MultiIndex
        Standard deviation of MSD per bin.
    dict
        Coordinate variation during measurement.
    """

    sel = univ.select_atoms(f"type {atom_type}")
    box = univ.dimensions

    if nblocks is None:
        nblocks = int((len(univ.trajectory) - corrlength) / gaplength)

    if (gaplength * nblocks + corrlength) > len(univ.trajectory):
        raise Exception(
            "Gap between correlation block, correlation length, or number of block too large."
        )

    if nblocks > 1:
        starts = np.arange(0, nblocks * gaplength, gaplength)
        ends = starts + corrlength
    else:
        starts = np.array([0])
        ends = np.array([corrlength])

    print(
        f"Compute the profile of the mean-squared displacement on {len(sel)} {atom_type} atoms along the {axis} axis..."
    )
    print(f"Total number of frames: {len(univ.trajectory)}")
    print(f"Timestep between frames: {dt} fs")
    print(f"Number of blocks: {nblocks}")
    print(f"Length of blocks: {corrlength} frames / {corrlength * dt * 1e-3} ps")
    print(f"Gap between origins: {gaplength} frames / {gaplength * dt * 1e-3} ps\n")

    if axis == "x":
        axid = 0
    if axis == "y":
        axid = 1
    if axis == "z":
        axid = 2

    # if dmin == None:
    #     dmin = 0
    # if dmax == None:
    #     dmax = box[axid]

    t = np.arange(0, corrlength, 1) * dt / 1000

    # positions_array = np.array([], dtype=np.float64).reshape(0, corrlength)
    # avg_pos = np.array([], dtype=np.float64)
    # std_pos = np.array([], dtype=np.float64)

    msd_per_block = {
        "xx": [],
        "yy": [],
        "zz": [],
        "xy": [],
        "xz": [],
        "yz": [],
    }

    msd_df_dic = {
        "xx": [],
        "yy": [],
        "zz": [],
        "xy": [],
        "xz": [],
        "yz": [],
    }

    std_df_dic = {
        "xx": [],
        "yy": [],
        "zz": [],
        "xy": [],
        "xz": [],
        "yz": [],
    }

    error_df_dic = {
        "xx": [],
        "yy": [],
        "zz": [],
        "xy": [],
        "xz": [],
        "yz": [],
    }

    for i, (start, end) in enumerate(tqdm(zip(starts, ends), total=len(starts)), 1):
        # initialize the trajectory
        univ.trajectory[start]

        # create the selection of atoms
        sel = univ.select_atoms(f"type {atom_type}")

        # get the reference positions at the start of the block
        rpos = sel.positions
        rposbyframe = sel.positions  # to avoid problem due to PBC

        # initialize total dr and msd
        d = np.zeros(rpos.shape)
        msd_values = {key: [] for key in msd_per_block.keys()}
        positions_intra_block = []
        positions_per_block = []

        for ts in univ.trajectory[start:end]:
            # collect positions
            pos = sel.positions

            # calculate the differential displacement w.r.t last frame positions
            dd = pos - rpos

            # for i in range(3):
            #     idx = dd[:,i] > 0.5 * box[i]
            #     dd[idx] -= boxv[i]
            #     idx = -dd[:,i] > 0.5 * box[i]
            #     dd[idx] += boxv[i]

            # take PBC into account
            dd -= box[:3].T * (dd / box[:3].T).round()

            # Add the the differential displacement to the total displacement
            d += dd

            msd_values["xx"].append(d[:, 0] * d[:, 0])
            msd_values["yy"].append(d[:, 1] * d[:, 1])
            msd_values["zz"].append(d[:, 2] * d[:, 2])
            msd_values["xy"].append(d[:, 0] * d[:, 1])
            msd_values["xz"].append(d[:, 0] * d[:, 2])
            msd_values["yz"].append(d[:, 1] * d[:, 2])

            positions_intra_block.append((rposbyframe + dd)[:, axid])

            # save last frame position
            rpos = pos
            rposbyframe += dd

        for key, values in msd_values.items():
            msd_per_block[key].append(np.array(values).T)
            # msd_dic[key].extend(np.array(values).T)
        # positions_array = np.append(positions_array, np.array(positions_intra_block).T, axis=0)
        positions_per_block.append(np.array(positions_intra_block).T)

    msd_dic = {
        key: np.concatenate(blocks, axis=0)  # shape (n_atoms_total, corrlength)
        for key, blocks in msd_per_block.items()
    }

    positions_array = np.concatenate(positions_per_block, axis=0)
    avg_pos = np.mean(positions_array, axis=1)
    drange = np.arange(np.min(positions_array), np.max(positions_array), delta)
    binned_pos = (drange[:-1] + drange[1:]) / 2

    # Dictionnary with variation in position during the measurement of MSD
    column_groups = {i: [] for i in range(len(drange))}
    for col_index, avg in enumerate(avg_pos):
        bin_index = np.digitize(avg, drange) - 1
        if bin_index == len(drange):
            bin_index -= 1
        column_groups[bin_index].append(col_index)

    coords_var_dict = {
        j: positions_array[column_groups[i], :].flatten()
        - np.median(positions_array[column_groups[i], :].flatten())
        for i, j in enumerate(np.round(binned_pos, 2))
    }

    def remove_empty_arrays(input_dict):
        keys_to_remove = [key for key, value in input_dict.items() if value.size == 0]
        for key in keys_to_remove:
            del input_dict[key]
        return input_dict

    coords_var_dict = remove_empty_arrays(coords_var_dict)

    # create the histogram for MSD
    for key, values in msd_dic.items():
        binned_msd = np.array(
            [
                np.mean(values[np.where((avg_pos > low) & (avg_pos <= high))], axis=0)
                for low, high in zip(drange[:-1], drange[1:])
            ]
        )
        binned_std_msd = np.array(
            [
                np.std(values[np.where((avg_pos > low) & (avg_pos <= high))], axis=0)
                for low, high in zip(drange[:-1], drange[1:])
            ]
        )

        binned_pos_final = binned_pos[~np.isnan(binned_msd).any(axis=1)]
        binned_msd = binned_msd[~np.isnan(binned_msd).any(axis=1)]
        binned_std_msd = binned_std_msd[~np.isnan(binned_std_msd).any(axis=1)]

        msd_df_dic[key] = pd.DataFrame(
            binned_msd.T, index=t, columns=np.round(binned_pos_final, 2)
        )
        std_df_dic[key] = pd.DataFrame(
            binned_std_msd.T, index=t, columns=np.round(binned_pos_final, 2)
        )

        msd_df_dic[key].index.name = "time"
        std_df_dic[key].index.name = "time"

    return pd.concat(msd_df_dic, axis=1), pd.concat(std_df_dic, axis=1), coords_var_dict


def diffusion_coefficient(df_msd: pd.DataFrame, start: float = 0.0) -> pd.DataFrame:
    (
        """Calculate the diffusion coefficients for a bulk solution from an MSD DataFrame.

    Parameters
    ----------
    df_msd : pd.DataFrame
        DataFrame containing the MSD data.
    start : float, optional
        Time at which the linear regression starts (ps).

    Returns
    -------
    pd.DataFrame
        Calculated diffusion coefficients and errors in m²/s.
    """
        """Calculate the diffusion coefficients for a bulk solution."""
    )
    df_msd = df_msd.loc[start:]
    t = df_msd.index

    dc_array = []
    error_array = []

    for c in df_msd.columns:
        msd = df_msd[c]

        m, ci_m = _linear_fit(t, msd)
        dc_array.append(m)
        error_array.append(ci_m)

    res = pd.DataFrame(
        [dc_array, error_array],
        index=["DC (m2/s)", "Error (m2/s)"],
        columns=df_msd.columns,
    )
    res *= 1e-8 / 2  # Conversion from A2/ps to m2/s

    dc_iso = res[["xx", "yy", "zz"]].iloc[0].mean()
    sd_iso = np.sqrt((res[["xx", "yy", "zz"]].iloc[1] ** 2).sum()) / 3

    res["3d"] = [dc_iso, sd_iso]

    return res


def diffusion_coefficient_profile(
    df_msd: pd.MultiIndex, df_std: pd.MultiIndex = None, start: float = 1.0
) -> pd.DataFrame:
    """Calculate the spatial diffusion coefficient profile from an MSD MultiIndex DataFrame.

    Parameters
    ----------
    df_msd : pd.MultiIndex
        DataFrame containing the MSD profiles.
    df_std : pd.MultiIndex, optional
        DataFrame containing the standard deviations.
    start : float, optional
        Time at which the linear regression starts (ps).

    Returns
    -------
    pd.DataFrame
        Calculated diffusion coefficients and errors per spatial bin.
    """
    df_msd = df_msd.loc[start:]
    t = df_msd.index

    if df_std is not None:
        df_std = df_std.loc[start:]

    dc_array = []
    error_array = []

    components = ["xx", "yy", "zz", "xy", "xz", "yz"]
    for comp in components:
        msd_comp = df_msd[comp]
        std_comp = df_std[comp] if df_std is not None else None

        for c in msd_comp.columns:
            msd = msd_comp[c]
            std = std_comp[c] if std_comp is not None else None

            m, ci_m = _linear_fit(t, msd, std)
            dc_array.append(m)
            error_array.append(ci_m)

    res = pd.DataFrame(
        [dc_array, error_array],
        index=["DC (m2/s)", "Error (m2/s)"],
        columns=df_msd.columns,
    )
    res *= 1e-8 / 2

    # Isotropic calculation on the MultiIndex
    dc_iso = res["xx"].iloc[0] + res["yy"].iloc[0] + res["zz"].iloc[0]
    dc_iso /= 3
    sd_iso = (
        np.sqrt(
            res["xx"].iloc[1] ** 2 + res["yy"].iloc[1] ** 2 + res["zz"].iloc[1] ** 2
        )
        / 3
    )

    iso_values = pd.DataFrame(
        [dc_iso, sd_iso], index=res.index, columns=res["xx"].columns
    )
    multi_index = pd.MultiIndex.from_product(
        [["iso"], iso_values.columns], names=df_msd.columns.names
    )
    iso_values.columns = multi_index

    res = pd.concat([res, iso_values], axis=1)
    return res


def plot_msd(df_msd: pd.DataFrame) -> None:
    """Plot the MSD values from a pandas DataFrame.

    Parameters
    ----------
    df_msd : pd.DataFrame
        DataFrame containing the MSD data to plot.
    """

    # Plot using Seaborn
    half_size = int(len(df_msd.columns) / 2)

    # Get the colormaps
    cmap = plt.get_cmap("coolwarm")
    cmap_r = plt.get_cmap("coolwarm_r")

    # Generate arrays of values from 0.0 to 1.0 to sample colors
    colors_forward = cmap(np.linspace(0, 1, half_size))

    if len(df_msd.columns) % 2 == 0:
        colors_reverse = cmap_r(np.linspace(0, 1, half_size))
    else:
        colors_reverse = cmap_r(np.linspace(0, 1, int(half_size + 1)))

    # Combine the lists (Matplotlib handles RGBA arrays natively, no need for hex string conversion)
    palettef = np.vstack([colors_forward, colors_reverse])

    # Dashes
    dashes = []
    for i in range(len(df_msd.columns)):
        if i == 0 or i == len(df_msd.columns) - 1:
            dashes.append("--")
        else:
            dashes.append("-")

    fig, ax = plt.subplots()
    for i, (msd, col) in enumerate(zip(df_msd.values.T, df_msd.columns)):
        ax.plot(
            df_msd.index,
            msd,
            color=palettef[i],
            linestyle=dashes[i],
            label=np.round(col, 2),
        )

    # Set labels and title
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"MSD ($\AA^2$)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid()

    # Legend labels
    num_cols = int(np.ceil(len(df_msd.columns) / 20))
    ax.legend(
        title=r"z ($\AA$)", bbox_to_anchor=(1.04, 0.5), loc="center left", ncols=num_cols
    )

    # Adjust the plot layout to accommodate the legend
    fig.subplots_adjust(right=1 - num_cols * 0.20)

    # Show the plot
    fig.tight_layout()
    fig.show()


def plot_diffusion_profile(input_df: pd.DataFrame) -> None:
    """Plot the spatial diffusion coefficient profile from a pandas DataFrame.

    Parameters
    ----------
    df_dc : pd.DataFrame
        DataFrame containing the diffusion coefficient profile.
    """

    if input_df.shape[1] > input_df.shape[0]:
        input_df = input_df.T

    t = input_df.index
    dc = input_df["DC (m2/s)"]
    error = input_df["Error (m2/s)"]

    plt.scatter(t, dc, color="royalblue")
    plt.fill_between(t, dc - error, dc + error, alpha=0.2, color="royalblue")

    imax = dc.argmax()

    # Set labels and title
    plt.xlabel(r"Distance ($\AA$)")
    plt.ylabel("Diffusion coefficient ($m^2/s$)")

    plt.ylim(dc.min(), dc.max() + dc.iloc[imax])
    plt.grid()
    plt.show()
