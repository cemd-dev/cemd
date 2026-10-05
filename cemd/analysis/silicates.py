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
from typing import Any

import MDAnalysis as mda
import numpy as np

from ..core.atomic_system import AtomicSystem

# Centralized Defaults
TYPES_PRESET = {
    "si_types": "Si",
    "o_types": "O Ob Osi Osih Oa Oah Oh Ow Oc Os Obs",
    "al_types": "Al",
    "ca_types": "Ca Cw",
}


def _check_elements(source: AtomicSystem) -> None:
    """
    Check if the system contains Silicon and Oxygen.

    Raises
    ------
    ValueError
        If the system does not contain both Silicon and Oxygen.
    """
    elements = source.elements

    has_si = any("Si" in elem for elem in elements)
    has_o = any("O" in elem for elem in elements)

    if not has_si:
        raise ValueError(
            f"System does not contain Silicon (Si). Found elements: {elements}"
        )

    if not has_o:
        raise ValueError(
            f"System does not contain Oxygen (O). Found elements: {elements}"
        )


def _process_analyze(source, si_types, o_types, al_types, ca_types, cutoff):
    """Perform the core silicate network analysis.

    Calculates chemical stoichiometric ratios, polymerization (Q^n) distribution,
    and the Average Chain Length (MCL).

    Parameters
    ----------
    source : mda.Universe
        The MDAnalysis universe to analyze.
    si_types : str
        String of atom types defined as Silicon.
    o_types : str
        String of atom types defined as Oxygen.
    al_types : str
        String of atom types defined as Aluminum.
    ca_types : str
        String of atom types defined as Calcium.
    cutoff : float
        The distance cutoff for neighbor selection.

    Returns
    -------
    dict
        Calculated network metrics including 'Ca/(Si+Al)', 'Al/Si', 'H2O/(Si+Al)',
        'MCL', and 'Qn_distribution'.
    """

    def safe_select(types_str):
        if not types_str.strip():
            return None
        return source.select_atoms(f"type {types_str}")

    sel_si = safe_select(si_types)
    sel_o = safe_select(o_types)
    sel_al = safe_select(al_types)
    sel_ca = safe_select(ca_types)

    # "H*h*" only ever matched types carrying a second, lowercase "h"
    # (i.e. "Hh"), so the water hydrogens -- "H" straight from a build,
    # "Hw" after the ClayFF/CSHFF rules -- were skipped and the H2O
    # ratios below came out as 0. Count every hydrogen-named type.
    sel_h = source.select_atoms("type H*")

    if sel_si is None or sel_o is None:
        raise ValueError(
            "The atom types for Silicon (Si) and Oxygen (O) must be defined and valid."
        )

    n_si = len(sel_si)
    n_al = len(sel_al) if sel_al is not None else 0
    n_ca = len(sel_ca) if sel_ca is not None else 0
    n_h = len(sel_h)

    # Calculation of chemical stoichiometric ratios
    denom = n_si + n_al
    ca_si_ratio = n_ca / n_si if n_si > 0 else 0.0
    h2o_si_ratio = (n_h / 2) / n_si if n_si > 0 else 0.0
    ca_m_ratio = n_ca / denom if denom > 0 else 0.0
    al_si_ratio = n_al / n_si if n_si > 0 else 0.0
    h2o_m_ratio = (n_h / 2) / denom if denom > 0 else 0.0

    # Polymerization Analysis (Q^n Units)
    cqsi = np.zeros(5)

    # Definition of network formers (Si + Al)
    al_text = f" {al_types}" if al_types.strip() else ""
    network_formers = source.select_atoms(f"type {si_types}{al_text}")
    nf_indices = set(network_formers.indices)
    o_indices = set(sel_o.indices)

    for s_idx in sel_si.indices:
        neighbors_o = source.select_atoms(f"around {cutoff} index {s_idx}")
        n_bridging = 0
        for o_idx in neighbors_o.indices:
            if o_idx in o_indices:
                potential_bridges = source.select_atoms(
                    f"around {cutoff} index {o_idx}"
                )
                # If oxygen touches ANOTHER network trainer, it's an oxygen bridge
                bridges = [
                    idx
                    for idx in potential_bridges.indices
                    if idx in nf_indices and idx != s_idx
                ]
                if len(bridges) > 0:
                    n_bridging += 1
        if n_bridging <= 4:
            cqsi[n_bridging] += 1

    pqsi = (cqsi / n_si) * 100

    # Calculation of Average Chain Length (MCL)
    mcl = 2 * (1 + pqsi[2] / pqsi[1]) if pqsi[1] > 1e-3 else float("inf")

    return {
        "Ca/(Si+Al)": ca_m_ratio,
        "Ca/Si": ca_si_ratio,
        "Al/Si": al_si_ratio,
        "H2O/(Si+Al)": h2o_m_ratio,
        "H2O/Si": h2o_si_ratio,
        "MCL": float(mcl),
        "Qn_distribution": pqsi,
    }


@singledispatch
def analyze_silicates(
    source, types_map: dict = None, cutoff: float = 1.85
) -> dict[str, Any]:
    """Entry point for analyzing silicate network structures.

    Merges provided custom atom types with presets and delegates to the appropriate
    analysis method based on input type (AtomicSystem or mda.Universe).

    Parameters
    ----------
    source : AtomicSystem or mda.Universe
        The source structure or trajectory to analyze.
    types_map : dict, optional
        Custom mapping of atom types for Si, O, Al, and Ca.
    cutoff : float, optional
        The distance cutoff for defining atomic neighbors (in Angstroms).

    Returns
    -------
    dict
        A dictionary containing chemical ratios, mean chain length (MCL),
        and the Q^n distribution.
    """
    raise TypeError(
        f"analyze_silicates() does not support source of type {type(source).__name__}; "
        "expected an AtomicSystem or an MDAnalysis Universe."
    )


# NB: `functools.singledispatch` calls the registered implementation below
# directly with whatever arguments the caller passed (it does not route
# through the generic function body above once a type is registered), so
# the preset-merging has to happen in each implementation, and the
# parameter must be named `types_map` here too -- a caller passing
# `types_map=...` by keyword (as documented) would otherwise hit
# `TypeError: got an unexpected keyword argument 'types_map'`.
@analyze_silicates.register(AtomicSystem)
def _(source: AtomicSystem, types_map: dict = None, cutoff: float = 1.85) -> dict[str, Any]:
    _check_elements(source)
    config = {**TYPES_PRESET, **(types_map or {})}
    u = source.to_mda()
    return _process_analyze(
        u,
        si_types=config.get("si_types", "Si"),
        o_types=config.get("o_types", "O"),
        al_types=config.get("al_types", ""),
        ca_types=config.get("ca_types", ""),
        cutoff=cutoff,
    )


@analyze_silicates.register(mda.Universe)
def _(source: mda.Universe, types_map: dict = None, cutoff: float = 1.85) -> dict[str, Any]:
    config = {**TYPES_PRESET, **(types_map or {})}
    return _process_analyze(
        source,
        si_types=config.get("si_types", "Si"),
        o_types=config.get("o_types", "O"),
        al_types=config.get("al_types", ""),
        ca_types=config.get("ca_types", ""),
        cutoff=cutoff,
    )
