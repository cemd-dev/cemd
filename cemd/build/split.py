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

import warnings
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from .._constants import MASSES_DICT
from ..core._format import lattice2vectors, vectors2lattice
from .base import BaseBuilder

if TYPE_CHECKING:
    from ..core.atomic_system import AtomicSystem
    from .solution import SolutionBuilder


# Interatomic distances (Å) below which a contact is considered a real,
# structure-bearing bond that a cut would sever. Deliberately restricted to
# short, essentially covalent contacts: ionic coordination shells (Ca...Ow,
# Na...Cl) are *not* listed, since a cut passing through them is physically
# fine and must not be counted as damage. Keys can be atom types (checked
# first) or elements (fallback) -- see `Splitter.find_broken_bonds`.
# Si-O and H-O reach 2.0 and 1.2 A, not 1.8 and 1.1, so that the
# stretched bonds of disordered or ReaxFF-relaxed structures are not missed:
# nothing unbonded comes that close (Si...O second shell ~3 A, O-H...O
# hydrogen bonds >= 1.5 A).
DEFAULT_BOND_CUTOFFS: dict[tuple[str, str], float] = {
    ("Si", "O"): 2.0,
    ("Al", "O"): 1.8,
    ("C", "O"): 1.6,
    ("S", "O"): 1.7,
    ("H", "O"): 1.2,
}

# Cation-anion distances (Å) of the first coordination shell. Ionic
# contacts are not bonds: a surface always cuts some of them. They are only
# used to put isolated cations on the side where they keep most of their
# shell (``keep_molecules=True``), to rank cut positions, and reported. The
# first member of each key is the cation, i.e. the atom that gets placed.
DEFAULT_IONIC_CUTOFFS: dict[tuple[str, str], float] = {
    ("Ca", "O"): 3.0,
    ("Mg", "O"): 2.6,
    ("Na", "O"): 3.0,
    ("K", "O"): 3.3,
}


@dataclass
class Splitter(BaseBuilder):
    """
    Split a system along an axis and optionally add a solution.

    A cut can sever real interatomic bonds. Bonds that exist explicitly in
    the system's topology are dropped automatically. Contacts that are
    structure-bearing but carry no explicit bond -- the Si-O framework
    under ClayFF/CSHFF, for instance -- are found geometrically instead,
    from ``bonds_dict``; use :meth:`count_broken_bonds` or
    :meth:`scan_broken_bonds` to choose where to cut, and ``repair=True``
    to rebuild what the cut breaks.

    With ``repair=True``, a hydrogen always goes with its covalent partner,
    and every other severed bond is hydrolysed: of an X-O bond (Si-O,
    Al-O, C-O...), X gets an OH back where its O was and the O gets an H
    -- a bridging Si-O-Si gives two silanols, a Si-O...Ca a silanol and a
    Ca-OH. The caps take the type and charge of an atom of the system in
    the same covalent environment (Osih/Hsi, Oh/Hh under CSHFF), and an
    explicit bond or angle wherever the system already has one between the
    same types. A severed bond without an O (C-C, Si-Si...) has no such
    repair: the split refuses it, see ``keep_molecules``.

    In a disordered solid there is often no plane that cuts no bond at
    all. With ``keep_molecules=True`` the cut follows the molecules
    instead of a plane: every covalent unit (explicit bonds plus the
    contacts of ``bonds_dict``: a silicate chain with its O and H, a water
    molecule, a hydroxide...) moves as a whole, on the side holding most of
    its mass, and every isolated cation (Ca, Na...) goes to the side where
    it keeps most of its coordination shell (``ionic_dict``). No covalent
    bond is then broken, whatever the cut position; ``coordinate=None``
    picks the position that cuts the fewest ionic contacts
    (:meth:`best_coordinate`).

    After :meth:`split`, ``split_report`` gives the contacts measured as
    broken on the actual geometry, before vs. after the split.

    Parameters
    ----------
    system : AtomicSystem
        The system to split.
    coordinate : float, optional
        Position along ``axis`` at which to cut, as a projection on the
        axis unit vector. ``None`` uses :meth:`best_coordinate`.
    axis : int or str, default=2
        Axis along which to split ('x'/'y'/'z' or 0/1/2).
    gap_size : float, default=20.0
        Width of the gap to open, in Å.
    bonds_dict : dict, optional
        ``{(type_or_element, type_or_element): cutoff}`` used to detect
        severed contacts. Defaults to :data:`DEFAULT_BOND_CUTOFFS`.
    keep_molecules : bool, default=False
        Move whole covalent units rather than a half-space (see above).
    repair : bool, default=False
        Rebuild the severed bonds (see above). Set here, it also applies
        to :meth:`find_broken_bonds` and the scans, which then count the
        bonds the repair will rebuild.
    ionic_dict : dict, optional
        ``{(cation, anion): cutoff}`` coordination shells, types or
        elements. Defaults to :data:`DEFAULT_IONIC_CUTOFFS`.

    Examples
    --------
    >>> # Simple split: cut at z=15, opening a 20 A gap there
    >>> result = Splitter(system, coordinate=15.0, axis='z', gap_size=20.0).split()

    >>> # Split with solution (fluent API)
    >>> result = (Splitter(system, coordinate=15.0, axis='z', gap_size=30.0)
    ...           .add_solution(blueprint, padding=2.0)
    ...           .split())

    >>> # Look for the least damaging cut, then split and cap what breaks
    >>> splitter = Splitter(system, coordinate=0.0, axis='z', gap_size=20.0)
    >>> scan = splitter.scan_broken_bonds(step=0.25)
    >>> best = scan.loc[scan["n_broken"].idxmin(), "coordinate"]
    >>> result = splitter.with_coordinate(best).split(repair=True)
    >>> splitter.repair_report["charge_change"]  # to neutralize afterwards

    >>> # Pore in a disordered solid: no covalent bond broken, cut placed
    >>> # where it costs the fewest ionic contacts, liquid in the pore
    >>> result = (Splitter(system, axis='z', gap_size=30.0, keep_molecules=True)
    ...           .add_solution(blueprint)
    ...           .split())
    """

    system: AtomicSystem
    coordinate: float | None = None
    axis: int | str = 2
    gap_size: float = 20.0
    bonds_dict: dict[tuple[str, str], float] | None = None
    keep_molecules: bool = False
    ionic_dict: dict[tuple[str, str], float] | None = None
    repair: bool = False

    # Solution options (default: no solution)
    _solution_blueprint: SolutionBuilder | None = field(default=None, repr=False)
    _padding: float = field(default=2.0, repr=False)
    _vacuum: float = field(default=0.0, repr=False)
    _has_solution: bool = field(default=False, repr=False)
    _solution_thickness: float | None = field(default=None, repr=False, init=False)

    # Internal state (non-serializable)
    _unit_norm: np.ndarray | None = field(default=None, repr=False, init=False)
    _vec_list: list | None = field(default=None, repr=False, init=False)
    _axis_name: str | None = field(default=None, repr=False, init=False)

    _contact_cache: dict = field(default_factory=dict, repr=False, init=False)
    _molecule_cache: tuple | None = field(default=None, repr=False, init=False)

    # Filled in by `split(repair=True)`
    repair_report: dict[str, int] | None = field(default=None, repr=False, init=False)
    # Filled in by `split()`
    split_report: dict | None = field(default=None, repr=False, init=False)

    def __post_init__(self):
        """Initialize derived attributes."""
        # Convert axis if string
        if isinstance(self.axis, str):
            axis_map = {"x": 0, "y": 1, "z": 2}
            self.axis = axis_map[self.axis.lower()]

        self._axis_name = ["x", "y", "z"][self.axis]

        # Compute axis information
        vec_list = list(lattice2vectors(self.system.box))
        target_vec = vec_list[self.axis]
        self._unit_norm = target_vec / np.linalg.norm(target_vec)
        self._vec_list = vec_list

    def add_solution(
        self, blueprint: SolutionBuilder, padding: float = 2.0, vacuum: float = 0.0
    ) -> Splitter:
        """
        Add a liquid solution in the gap.

        Parameters
        ----------
        blueprint : SolutionBuilder
            The solution blueprint defining composition.
        padding : float, default=2.0
            Empty space between solution and surfaces (Å).
        vacuum : float, default=0.0
            Additional vacuum above the solution (Å).

        Returns
        -------
        Splitter
            Self for method chaining.
        """
        self._solution_blueprint = blueprint
        self._padding = padding
        self._vacuum = vacuum
        self._has_solution = True
        return self

    def with_coordinate(self, coordinate: float) -> Splitter:
        """
        Set the cutting position.

        Parameters
        ----------
        coordinate : float
            Position along the axis where to cut.

        Returns
        -------
        Splitter
            Self for method chaining.
        """
        self.coordinate = coordinate
        return self

    # ------------------------------------------------------------------
    # Severed-contact detection
    # ------------------------------------------------------------------

    def _resolve_cutoffs(self, bonds_dict: dict | None) -> dict:
        """Return the cutoff table to use, falling back to the defaults."""
        if bonds_dict is not None:
            return bonds_dict
        if self.bonds_dict is not None:
            return self.bonds_dict
        return DEFAULT_BOND_CUTOFFS

    @staticmethod
    def _lookup_cutoff(
        cutoffs: dict, type_1: str, type_2: str, element_1: str, element_2: str
    ) -> float | None:
        """Find the cutoff for a pair: atom types first, elements as fallback.

        Types take priority so that a system already carrying force-field
        types can distinguish structure-bearing oxygens (``Ob``, ``Osih``)
        from ones that merely sit nearby (``Ow``), which elements alone
        cannot express.
        """
        for key in ((type_1, type_2), (type_2, type_1)):
            if key in cutoffs:
                return cutoffs[key]

        for key in ((element_1, element_2), (element_2, element_1)):
            if key in cutoffs:
                return cutoffs[key]

        return None

    def _contacts(self, cutoffs: dict) -> dict[str, np.ndarray]:
        """Pairs closer than their cutoff, on the current geometry.

        Returns positional indices ``i``, ``j``, the distance ``d`` and
        ``k``, the number of periodic images along the split axis crossed
        by the minimum-image vector going from ``i`` to ``j``. The split
        moves some atoms by ``gap_size`` and grows the box by as much, so a
        pair keeps its separation iff ``moved[j] - moved[i] == k``.

        Cached: positions do not change until the split itself.
        """
        from MDAnalysis.lib.distances import self_capped_distance

        key = tuple(sorted(cutoffs.items()))
        if key in self._contact_cache:
            return self._contact_cache[key]

        empty = np.array([], dtype=int)
        atoms = self.system.atoms
        result = {"i": empty, "j": empty, "d": np.array([]), "k": empty}
        if not cutoffs or len(atoms) < 2:
            self._contact_cache[key] = result
            return result

        positions = atoms[["x", "y", "z"]].to_numpy(dtype=np.float64)
        box = np.asarray(self.system.box, dtype=np.float64)
        types = atoms["type"].astype(str).to_numpy()
        elements = self.system.elements
        atom_elements = np.array([elements.get(t, t) for t in types])

        pairs, distances = self_capped_distance(
            positions, max_cutoff=max(cutoffs.values()), box=box
        )
        i, j = pairs[:, 0], pairs[:, 1]

        # One lookup per distinct pair of types, not per pair of atoms
        pair_cutoff = np.full(len(pairs), -1.0)
        type_pairs = pd.Series(list(zip(types[i], types[j])))
        for (type_1, type_2), index in type_pairs.groupby(type_pairs).groups.items():
            index = np.asarray(index)
            cutoff = self._lookup_cutoff(
                cutoffs,
                type_1,
                type_2,
                atom_elements[i[index[0]]],
                atom_elements[j[index[0]]],
            )
            if cutoff is not None:
                pair_cutoff[index] = cutoff
        keep = distances <= pair_cutoff
        i, j, distances = i[keep], j[keep], distances[keep]

        result = {"i": i, "j": j, "d": distances, "k": self._image_shift(i, j)}
        self._contact_cache[key] = result
        return result

    def _image_shift(self, i: np.ndarray, j: np.ndarray) -> np.ndarray:
        """Periodic images crossed along the split axis, from ``i`` to ``j``."""
        from MDAnalysis.lib.distances import minimize_vectors

        if len(i) == 0:
            return np.array([], dtype=int)

        positions = self.system.atoms[["x", "y", "z"]].to_numpy(dtype=np.float64)
        raw = positions[j] - positions[i]
        nearest = minimize_vectors(raw, np.asarray(self.system.box, dtype=np.float64))
        # raw - nearest is a lattice vector: express it in the lattice basis
        lattice = np.array(self._vec_list, dtype=np.float64)
        images = (raw - nearest) @ np.linalg.inv(lattice)
        return np.rint(images[:, self.axis]).astype(int)

    def _explicit_contacts(self) -> dict[str, np.ndarray]:
        """Explicit bonds, in the same form as :meth:`_contacts`."""
        key = "explicit"
        if key in self._contact_cache:
            return self._contact_cache[key]

        bonds = self.system.bonds
        empty = np.array([], dtype=int)
        if bonds is None or bonds.empty:
            result = {"i": empty, "j": empty, "k": empty}
        else:
            index = pd.Series(
                np.arange(len(self.system.atoms)), index=self.system.atoms.index
            )
            i = index.loc[bonds["atom_1"].astype(int)].to_numpy()
            j = index.loc[bonds["atom_2"].astype(int)].to_numpy()
            result = {"i": i, "j": j, "k": self._image_shift(i, j)}

        self._contact_cache[key] = result
        return result

    def _covalent_contacts(self, bonds_dict: dict | None = None) -> dict[str, np.ndarray]:
        """Explicit bonds and geometric contacts of ``bonds_dict``, deduplicated."""
        geometric = self._contacts(self._resolve_cutoffs(bonds_dict))
        explicit = self._explicit_contacts()

        i = np.concatenate([geometric["i"], explicit["i"]])
        j = np.concatenate([geometric["j"], explicit["j"]])
        k = np.concatenate([geometric["k"], explicit["k"]])

        # Orient every pair low -> high index, so duplicates coincide
        swap = i > j
        i, j = np.where(swap, j, i), np.where(swap, i, j)
        k = np.where(swap, -k, k)
        _, unique = np.unique(np.stack([i, j]), axis=1, return_index=True)
        return {"i": i[unique], "j": j[unique], "k": k[unique]}

    def _molecules(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Covalent units and the relative side of their atoms.

        Returns, per atom, the component label and ``rel``, the side label
        relative to the lowest one in the component: walking any bond
        ``i -> j`` adds ``k_ij`` (see :meth:`_contacts`), so that moving
        the atoms with ``root + rel == 1`` keeps every bond of the unit.
        Also returns, per component, the span of ``rel`` (0 or 1).
        """
        if self._molecule_cache is not None:
            return self._molecule_cache

        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import breadth_first_order, connected_components

        n_atoms = len(self.system.atoms)
        contacts = self._covalent_contacts()
        i, j, k = contacts["i"], contacts["j"], contacts["k"]

        graph = coo_matrix((np.ones(len(i)), (i, j)), shape=(n_atoms, n_atoms)).tocsr()
        n_components, labels = connected_components(graph, directed=False)

        # k along each directed edge, both directions
        shift = {}
        for a, b, kk in zip(i.tolist(), j.tolist(), k.tolist()):
            shift[(a, b)] = kk
            shift[(b, a)] = -kk

        rel = np.zeros(n_atoms, dtype=int)
        seen = np.zeros(n_components, dtype=bool)
        for start in range(n_atoms):
            if seen[labels[start]]:
                continue
            seen[labels[start]] = True
            order, predecessors = breadth_first_order(
                graph, start, directed=False, return_predecessors=True
            )
            for atom in order[1:]:
                parent = predecessors[atom]
                rel[atom] = rel[parent] + shift[(parent, atom)]

        # A unit closed on itself through the boundary along the axis (an
        # infinite chain, a bulk framework) cannot be kept whole by a cut.
        inconsistent = rel[j] - rel[i] != k
        if inconsistent.any():
            atom_id = self.system.atoms.index[i[inconsistent][0]]
            raise ValueError(
                f"The covalent unit of atom {atom_id} is periodic along "
                f"{self._axis_name}: no cut can keep it whole. Restrict "
                "`bonds_dict`, or split without `keep_molecules`."
            )

        rel_min = np.full(n_components, np.iinfo(int).max)
        np.minimum.at(rel_min, labels, rel)
        rel -= rel_min[labels]
        span = np.zeros(n_components, dtype=int)
        np.maximum.at(span, labels, rel)
        if span.max(initial=0) > 1:
            raise ValueError(
                "A covalent unit is longer than the cell along "
                f"{self._axis_name}: no cut can keep it whole."
            )

        self._molecule_cache = (labels, rel, span)
        return self._molecule_cache

    def side_mask(self, coordinate: float | None = None) -> np.ndarray:
        """Atoms moved by the split (by position in ``system.atoms``).

        A half-space for a plain split; whole covalent units and their
        cations with ``keep_molecules=True``. With ``repair``, every
        hydrogen goes on the side of its covalent partner.
        """
        coordinate = self._resolve_coordinate(coordinate)
        positions = self.system.atoms[["x", "y", "z"]].to_numpy(dtype=np.float64)
        moves = (positions @ self._unit_norm >= coordinate).astype(int)
        if not self.keep_molecules:
            if self.repair:
                moves = self._attach_hydrogens(moves)
            return moves.astype(bool)

        labels, rel, span = self._molecules()
        n_components = len(span)
        masses = self.system.atoms["type"].astype(str).map(self.system.masses)
        masses = masses.fillna(1.0).to_numpy(dtype=np.float64)

        # Root label of each unit: the one agreeing with the half-space for
        # most of its mass; forced to 0 when the unit spans both labels.
        agreement = np.zeros((n_components, 2))
        for root in (0, 1):
            np.add.at(agreement[:, root], labels, masses * (moves == root + rel))
        root = (agreement[:, 1] > agreement[:, 0]).astype(int)
        root[span == 1] = 0
        moves = root[labels] + rel

        # Isolated cations follow most of their coordination shell
        ionic = self._contacts(self._resolve_ionic())
        if len(ionic["i"]):
            sizes = np.bincount(labels, minlength=n_components)
            types = self.system.atoms["type"].astype(str).to_numpy()
            elements = self.system.elements
            cations = {key[0] for key in self._resolve_ionic()}
            is_cation = np.array(
                [t in cations or elements.get(t, t) in cations for t in types]
            ) & (sizes[labels] == 1)

            i, j, k = ionic["i"], ionic["j"], ionic["k"]
            cation_first = is_cation[i]
            cation = np.where(cation_first, i, j)
            anion = np.where(cation_first, j, i)
            shift = np.where(cation_first, k, -k)  # needed: m_anion - m_cation
            valid = is_cation[cation] & ~is_cation[anion]

            kept = np.zeros((len(moves), 2), dtype=int)
            for side in (0, 1):
                np.add.at(
                    kept[:, side],
                    cation[valid],
                    (moves[anion[valid]] - side) == shift[valid],
                )
            index = np.arange(len(moves))
            flip = is_cation & (kept[index, 1 - moves] > kept[index, moves])
            moves[flip] = 1 - moves[flip]

        if self.repair:
            moves = self._attach_hydrogens(moves)
        return moves.astype(bool)

    def _atom_elements(self) -> np.ndarray:
        """Element of each atom, its type where the mass names none."""
        types = self.system.atoms["type"].astype(str).to_numpy()
        elements = self.system.elements
        return np.array([elements.get(t, t) for t in types])

    def _attach_hydrogens(self, moves: np.ndarray) -> np.ndarray:
        """Put every hydrogen on the side of its covalent partner.

        Cutting an O-H leaves a bare proton on one side, which no repair
        can make sense of; the H goes with its partner instead.
        """
        contacts = self._covalent_contacts()
        i, j, k = contacts["i"], contacts["j"], contacts["k"]
        is_h = self._atom_elements() == "H"
        moves = moves.copy()
        # A pair keeps its separation iff moves[j] - moves[i] == k
        for h, partner, shift in ((j, i, k), (i, j, -k)):
            pick = is_h[h] & ~is_h[partner]
            target = moves[partner[pick]] + shift[pick]
            fits = (target == 0) | (target == 1)
            moves[h[pick][fits]] = target[fits]
        return moves

    def _resolve_ionic(self) -> dict:
        return DEFAULT_IONIC_CUTOFFS if self.ionic_dict is None else self.ionic_dict

    def _resolve_coordinate(self, coordinate: float | None) -> float:
        if coordinate is not None:
            return coordinate
        if self.coordinate is None:
            self.coordinate = self.best_coordinate()
        return self.coordinate

    @staticmethod
    def _broken(moves: np.ndarray, contacts: dict) -> np.ndarray:
        m = moves.astype(int)
        return m[contacts["j"]] - m[contacts["i"]] != contacts["k"]

    def find_broken_bonds(
        self,
        coordinate: float | None = None,
        bonds_dict: dict | None = None,
    ) -> list[tuple[int, int]]:
        """Find the bonds a cut would sever, from the current geometry.

        Looks at interatomic *distances* (``bonds_dict``), so it also finds
        structure-bearing contacts that carry no explicit bond (the Si-O
        framework under ClayFF, for instance), and at the explicit bonds
        of the topology.

        Parameters
        ----------
        coordinate : float, optional
            Cut position to test. Defaults to this splitter's own
            ``coordinate``.
        bonds_dict : dict, optional
            Cutoff table, see :data:`DEFAULT_BOND_CUTOFFS`.

        Returns
        -------
        list of tuple of int
            Pairs of atom ids whose contact the cut would break.

        Notes
        -----
        Call this *before* splitting: it reports on the geometry as it
        currently stands. Contacts that are neighbours through the
        periodic boundary survive when both atoms move with their image
        (the split shifts the box and the moving fragment by the same
        amount). With ``keep_molecules=True``, nothing covalent breaks.
        """
        moves = self.side_mask(coordinate)
        contacts = self._covalent_contacts(bonds_dict)
        broken = self._broken(moves, contacts)

        ids = self.system.atoms.index.to_numpy()
        return [
            (int(ids[a]), int(ids[b]))
            for a, b in zip(contacts["i"][broken], contacts["j"][broken])
        ]

    def count_broken_bonds(
        self,
        coordinate: float | None = None,
        bonds_dict: dict | None = None,
    ) -> int:
        """Count the bonds a cut at ``coordinate`` would sever.

        See :meth:`find_broken_bonds` for the detection rules.
        """
        return len(self.find_broken_bonds(coordinate, bonds_dict))

    def count_broken_ionic(self, coordinate: float | None = None) -> int:
        """Count the ionic contacts (``ionic_dict``) a cut would sever."""
        moves = self.side_mask(coordinate)
        return int(self._broken(moves, self._contacts(self._resolve_ionic())).sum())

    def describe_cut(
        self,
        coordinate: float | None = None,
        bonds_dict: dict | None = None,
    ) -> dict:
        """What a cut at ``coordinate`` moves and severs, in one pass.

        :meth:`find_broken_bonds` and :meth:`count_broken_ionic` each
        place the atoms again; with ``keep_molecules`` that placement is
        most of the cost, so a caller wanting both (a live preview) gets
        them here from a single :meth:`side_mask`.

        Returns
        -------
        dict
            ``moves``, the :meth:`side_mask`; ``broken``, the pairs of
            atom ids of the severed bonds, as :meth:`find_broken_bonds`;
            ``n_ionic``, the number of severed ionic contacts;
            ``unrepairable``, the pairs of ``broken`` that ``repair``
            refuses (no O, see :meth:`split`).
        """
        moves = self.side_mask(coordinate)
        contacts = self._covalent_contacts(bonds_dict)
        broken = self._broken(moves, contacts)
        ids = self.system.atoms.index.to_numpy()
        i, j = contacts["i"][broken], contacts["j"][broken]
        repairable = self._repairable(i, j)
        return {
            "moves": moves,
            "broken": [(int(ids[a]), int(ids[b])) for a, b in zip(i, j)],
            "unrepairable": [
                (int(ids[a]), int(ids[b]))
                for a, b in zip(i[~repairable], j[~repairable])
            ],
            "n_ionic": int(
                self._broken(moves, self._contacts(self._resolve_ionic())).sum()
            ),
        }

    def scan_broken_bonds(
        self,
        start: float | None = None,
        stop: float | None = None,
        step: float = 0.5,
        bonds_dict: dict | None = None,
    ) -> pd.DataFrame:
        """Count severed contacts across a range of cut positions.

        Useful to pick a cut plane that damages the structure as little as
        possible before committing to it.

        Parameters
        ----------
        start, stop : float, optional
            Range of cut coordinates to scan. Defaults to the span of the
            system along the split axis.
        step : float, default=0.5
            Scanning increment in Å.
        bonds_dict : dict, optional
            Cutoff table, see :data:`DEFAULT_BOND_CUTOFFS`.

        Returns
        -------
        pandas.DataFrame
            Columns ``coordinate``, ``n_broken`` (bonds) and ``n_ionic``
            (ionic contacts), one row per position.
        """
        if step <= 0:
            raise ValueError(f"step must be positive, got {step}")

        positions = self.system.atoms[["x", "y", "z"]].to_numpy(dtype=np.float64)
        projections = positions @ self._unit_norm

        start = float(projections.min()) if start is None else start
        stop = float(projections.max()) if stop is None else stop

        coordinates = np.arange(start, stop + step / 2, step)

        return pd.DataFrame(
            {
                "coordinate": coordinates,
                "n_broken": [
                    self.count_broken_bonds(c, bonds_dict) for c in coordinates
                ],
                "n_ionic": [self.count_broken_ionic(c) for c in coordinates],
            }
        )

    def best_coordinate(
        self,
        start: float | None = None,
        stop: float | None = None,
        step: float = 0.25,
    ) -> float:
        """Cut position breaking the fewest bonds, then the fewest ionic contacts.

        Remaining ties (with ``keep_molecules``, every position breaks no
        bond) go to the plane crossing the fewest bonds itself, i.e. the
        one dragging the fewest units across: the flattest surfaces.

        See :meth:`scan_broken_bonds` for the parameters.
        """
        return self.best_in_scan(self.scan_broken_bonds(start, stop, step))

    def best_in_scan(self, scan: pd.DataFrame) -> float:
        """The position :meth:`best_coordinate` picks, from a scan already made.

        ``scan`` is a :meth:`scan_broken_bonds` table of this splitter.
        """
        contacts = self._covalent_contacts()
        positions = self.system.atoms[["x", "y", "z"]].to_numpy(dtype=np.float64)
        projections = positions @ self._unit_norm
        # Not written into `scan`: it belongs to the caller.
        scan = scan.assign(
            n_crossing=[
                int(self._broken(projections >= c, contacts).sum())
                for c in scan["coordinate"]
            ]
        )

        best = scan.sort_values(
            ["n_broken", "n_ionic", "n_crossing"], kind="stable"
        ).iloc[0]
        return float(best["coordinate"])

    def _move_fragments(self, universe) -> None:
        """Move the second fragment to create the gap."""
        pos = universe.atoms.positions
        mask = self.side_mask()
        shift_vector = self._unit_norm * self.gap_size
        pos[mask] += shift_vector
        universe.atoms.positions = pos

    def _remove_crossing_bonds(self, original_positions, old_box) -> None:
        """Remove bonds stretched across the newly created gap.

        A bond whose two atoms ended up on opposite sides of the cut gets
        physically severed by the split: its length grows by roughly
        ``gap_size``. Compare each bond's length before/after the move to
        detect this, then drop the bond and any angle/dihedral/improper
        built on it (they'd otherwise reference a pair of atoms that are
        no longer meaningfully connected).

        Both lengths are computed under the minimum-image convention --
        ``old_box`` for the intact geometry, the (already grown)
        ``self.system.box`` for the post-move one. A bond is free to wrap
        through a periodic boundary that has nothing to do with the cut
        (most systems have at least one, simply from how the cell was
        built); measured as a raw Cartesian difference such a bond looks
        artificially long *before* the move already, and moving only one
        of its two atoms by ``gap_size`` then falsely reads as "stretched
        by the cut" even though neither atom is anywhere near the cut
        plane.
        """
        bonds = self.system.bonds
        if bonds is None or bonds.empty:
            return

        from MDAnalysis.lib.distances import calc_bonds

        current_positions = self.system.atoms[["x", "y", "z"]]

        a1_ids = bonds["atom_1"].astype(int).to_numpy()
        a2_ids = bonds["atom_2"].astype(int).to_numpy()

        orig1 = original_positions.loc[a1_ids].to_numpy(dtype=np.float32)
        orig2 = original_positions.loc[a2_ids].to_numpy(dtype=np.float32)
        cur1 = current_positions.loc[a1_ids].to_numpy(dtype=np.float32)
        cur2 = current_positions.loc[a2_ids].to_numpy(dtype=np.float32)

        d_before = calc_bonds(orig1, orig2, box=np.asarray(old_box, dtype=np.float32))
        d_after = calc_bonds(
            cur1, cur2, box=np.asarray(self.system.box, dtype=np.float32)
        )

        is_broken = d_after > d_before + self.gap_size / 2
        broken_pairs = {
            frozenset((int(a1), int(a2)))
            for a1, a2, broken in zip(a1_ids, a2_ids, is_broken)
            if broken
        }

        if not broken_pairs:
            return

        def _references_broken_pair(row, n_atoms) -> bool:
            atom_cols = [f"atom_{i}" for i in range(1, n_atoms + 1)]
            indices = [int(row[c]) for c in atom_cols]
            return any(
                frozenset((indices[i], indices[i + 1])) in broken_pairs
                for i in range(len(indices) - 1)
            )

        for name, n_atoms in [
            ("bonds", 2),
            ("angles", 3),
            ("dihedrals", 4),
            ("impropers", 4),
        ]:
            df = getattr(self.system, name)
            if df is None or df.empty:
                continue
            cleaned = df[~df.apply(_references_broken_pair, axis=1, args=(n_atoms,))]
            setattr(self.system, name, cleaned if not cleaned.empty else None)

    def _update_box(self) -> None:
        """Update the simulation box after splitting."""
        self._vec_list[self.axis] = (
            self._vec_list[self.axis] + self._unit_norm * self.gap_size
        )
        new_box = vectors2lattice(tuple(self._vec_list))
        self.system.set_box(new_box)

    def _pore_walls(self) -> tuple[float, float]:
        """Positions of the two walls of the opened gap, along the axis.

        The outermost atoms of each side, once the fragments have moved.
        For a plane cut these are the cut plane and the plane ``gap_size``
        above it; with ``keep_molecules`` the units that straddled the
        plane stick out into the gap, and the liquid must stay clear of
        them.
        """
        low = self.coordinate
        high = self.coordinate + self.gap_size
        positions = self.system.atoms[["x", "y", "z"]].to_numpy(dtype=np.float64)
        projections = positions @ self._unit_norm
        inside = projections[(projections > low) & (projections < high)]

        middle = (low + high) / 2
        lower, upper = inside[inside < middle], inside[inside >= middle]
        if len(lower):
            low = float(lower.max())
        if len(upper):
            high = float(upper.min())
        return low, high

    def _calculate_solution_thickness(self, width: float | None = None) -> float:
        """Calculate the optimal thickness for the liquid.

        ``width`` is the free space between the walls, ``gap_size`` when
        not given.
        """
        width = self.gap_size if width is None else width
        thickness = width - 2 * self._padding - self._vacuum
        if thickness <= 0:
            raise ValueError(
                f"Gap size ({width:.1f} Å free) too small for "
                f"padding ({self._padding:.1f} Å) and vacuum ({self._vacuum:.1f} Å). "
                f"Available space: {thickness:.1f} Å"
            )
        return thickness

    def _insert_solution(self) -> None:
        """Fill the gap opened by the split with the solution.

        The box is already the final one here (``_update_box`` has grown it
        by ``gap_size``), so the liquid is merged *into* it rather than
        stacked onto the system: ``_merge_structure`` -- which is what the
        surface builder uses to lay a film on top of a slab -- would have
        translated the liquid above the topmost atom and stretched the cell
        a second time, leaving the pore empty and the fluid outside the
        solid.
        """
        from .interface import _calculate_liquid_box, _merge_data

        low, high = self._pore_walls()
        liquid_thickness = self._calculate_solution_thickness(high - low)
        self._solution_thickness = liquid_thickness

        # Calculate the box for the liquid
        liquid_box = _calculate_liquid_box(
            self.system, liquid_thickness, self._axis_name
        )

        # Build the liquid
        liquid = self._solution_blueprint.build(liquid_box)

        # Packmol packs the liquid in its own cell starting at the origin;
        # the transverse dimensions match the host cell, so the liquid only
        # has to be moved onto the host's own origin.
        bounds = self.system._box_lmp
        for index, name in enumerate("xyz"):
            if index == self.axis:
                continue
            liquid.atoms[name] += float(bounds[index][0])

        # Center the liquid between the walls
        gap_center = (low + high) / 2.0
        liq_pos = liquid.atoms[self._axis_name].values
        liq_center = (liq_pos.max() + liq_pos.min()) / 2.0
        liquid.atoms[self._axis_name] += gap_center - liq_center

        # Merge everything, keeping the box set by `_update_box`
        result = _merge_data(self.system, liquid, self.system._box_lmp)

        # Replace the system with the merged result
        self.system._replace_internals(result)

    # ------------------------------------------------------------------
    # Repair
    # ------------------------------------------------------------------

    #: Closest a cap may come to an atom other than the one it caps, in Å.
    CAP_CLEARANCE = {"H": 1.4, "O": 2.0}

    #: Tilts (degrees) away from the broken bond tried, in order, when a
    #: cap's spot is taken; twelve directions around the bond each.
    CAP_TILTS = (15, 30, 45, 60, 75, 90, 105, 120)

    def _repairable(self, i: np.ndarray, j: np.ndarray) -> np.ndarray:
        """Which severed pairs a hydrolysis can rebuild: exactly one O."""
        elements = self._atom_elements()
        return (elements[i] == "O") != (elements[j] == "O")

    def _plan_repair(self, bonds_dict: dict | None) -> list[tuple]:
        """List the caps to place, from the intact geometry.

        Refuses the cut, before anything moves, if it severs a bond that
        cannot be hydrolysed.

        Returns
        -------
        list of (atom id, cap, unit vector toward the lost partner, length)
            ``cap`` is ``"H"`` for an O that lost its partner, ``"OH"``
            for the partner that lost its O. Translation leaves the
            vector valid on the split geometry.
        """
        from MDAnalysis.lib.distances import minimize_vectors

        moves = self.side_mask()
        contacts = self._covalent_contacts(bonds_dict)
        broken = self._broken(moves, contacts)
        i, j = contacts["i"][broken], contacts["j"][broken]

        elements = self._atom_elements()
        ids = self.system.atoms.index.to_numpy()
        repairable = self._repairable(i, j)
        if not repairable.all():
            pairs = pd.Series(
                [
                    "-".join(sorted((elements[a], elements[b])))
                    for a, b in zip(i[~repairable], j[~repairable])
                ]
            ).value_counts()
            raise ValueError(
                f"The cut severs {int((~repairable).sum())} bonds that "
                f"cannot be repaired ({pairs.to_dict()}): only a bond to an "
                "O can be hydrolysed. Move the cut, or use "
                "keep_molecules=True to keep these molecules whole."
            )

        positions = self.system.atoms[["x", "y", "z"]].to_numpy(dtype=np.float64)
        box = np.asarray(self.system.box, dtype=np.float64)
        bond = minimize_vectors(positions[j] - positions[i], box)
        length = np.linalg.norm(bond, axis=1)
        direction = bond / length[:, None]

        plan = []
        for a, b, d, r in zip(i, j, direction, length):
            # Seen from each atom: toward the partner it loses
            for atom, toward in ((a, d), (b, -d)):
                cap = "H" if elements[atom] == "O" else "OH"
                plan.append((int(ids[atom]), cap, toward, float(r)))
        return plan

    def _environment(self, bonds_dict: dict | None) -> tuple[list, list]:
        """Covalent environment of every atom, at two depths.

        The first depth is the element and its partners' elements; the
        second replaces each partner by its own first depth, which tells
        the H of a water (on an O-H-H) from the H of a silanol (on an
        O-Si-H).
        """
        contacts = self._covalent_contacts(bonds_dict)
        elements = self._atom_elements()
        partners: list[list[int]] = [[] for _ in range(len(elements))]
        for a, b in zip(contacts["i"], contacts["j"]):
            partners[a].append(b)
            partners[b].append(a)
        first = [
            (elements[a], tuple(sorted(elements[p] for p in partners[a])))
            for a in range(len(elements))
        ]
        second = [
            (elements[a], tuple(sorted(first[p] for p in partners[a])))
            for a in range(len(elements))
        ]
        return first, second

    def _templates(self, bonds_dict: dict | None) -> tuple[dict, dict, dict]:
        """Type and charge of each covalent environment of the system.

        Read from the intact system: its atoms are the models for the
        ones the repair adds or changes. Each environment keeps its most
        common type, with that type's mean charge there. The third table
        is each atom's own environment, by id, to tell later which ones
        the repair changed.
        """
        first, second = self._environment(bonds_dict)
        atoms = self.system.atoms
        tables = []
        for environment in (first, second):
            frame = pd.DataFrame(
                {
                    "env": environment,
                    "type": atoms["type"].to_numpy(),
                    "charge": atoms["charge"].to_numpy(dtype=np.float64),
                }
            )
            table = {}
            for env, group in frame.groupby("env", sort=False):
                atype = group["type"].value_counts().index[0]
                charge = group.loc[group["type"] == atype, "charge"].mean()
                table[env] = (atype, float(charge))
            tables.append(table)
        return tables[0], tables[1], dict(zip(atoms.index, second))

    def _place_caps(
        self, plan: list[tuple], oh_length: float
    ) -> tuple[list, list, int]:
        """Positions of the caps, clear of the atoms around them.

        A cap goes along the broken bond, at the bond length for an OH's O
        and ``oh_length`` for an H. When an atom (or a cap already placed)
        sits closer than :attr:`CAP_CLEARANCE`, the cap tilts away from
        the bond, as little as frees it; with no free direction, it takes
        the roomiest one and counts as crowded.

        Returns
        -------
        caps : list of (element, position, index of its partner)
            The partner is an atom id, or the position in ``caps`` of the
            cap it is bonded to (the O of an OH, for its H), as ``-1 - i``.
        bonds : list of (partner, cap) in the same notation
        crowded : int
        """
        from MDAnalysis.lib.distances import capped_distance, minimize_vectors

        if not plan:
            return [], [], 0

        atoms = self.system.atoms
        positions = atoms[["x", "y", "z"]].to_numpy(dtype=np.float64)
        box = np.asarray(self.system.box, dtype=np.float64)
        row = {atom_id: n for n, atom_id in enumerate(atoms.index)}

        anchors = np.array([row[step[0]] for step in plan], dtype=int)
        reach = (
            oh_length
            + 2 * max(step[3] for step in plan)
            + max(self.CAP_CLEARANCE.values())
        )
        near: list[list[int]] = [[] for _ in plan]
        pairs = capped_distance(
            positions[anchors],
            positions,
            max_cutoff=reach,
            box=box,
            return_distances=False,
        )
        for n, other in pairs:
            near[n].append(other)

        caps: list[tuple[str, np.ndarray, int]] = []
        bonds: list[tuple[int, int]] = []
        crowded = 0

        # The tilts and turns tried around a bond, the bond itself first
        tilt = np.radians(np.repeat([0, *self.CAP_TILTS], 12)[11:])
        turn = np.tile(
            np.linspace(0, 2 * np.pi, 12, endpoint=False), len(self.CAP_TILTS) + 1
        )[11:]

        def directions(axis: np.ndarray) -> np.ndarray:
            """`axis`, then cones of growing tilt around it."""
            helper = np.eye(3)[np.argmin(np.abs(axis))]
            u = np.cross(axis, helper)
            u /= np.linalg.norm(u)
            v = np.cross(axis, u)
            side = np.cos(turn)[:, None] * u + np.sin(turn)[:, None] * v
            return np.cos(tilt)[:, None] * axis + np.sin(tilt)[:, None] * side

        def place(element, origin, axis, length, others, own=None) -> np.ndarray:
            """Best spot for a cap; `own` is the cap it bonds to, if any."""
            nonlocal crowded
            candidates = origin + length * directions(axis)
            # Caps placed so far, the nearby ones only
            placed = np.array([c[1] for n, c in enumerate(caps) if n != own])
            if len(placed):
                offset = minimize_vectors(placed - origin, box)
                placed = placed[np.linalg.norm(offset, axis=1) < reach]
            reference = np.vstack([positions[others], placed.reshape(-1, 3)])
            if len(reference) == 0:
                return candidates[0]
            gaps = minimize_vectors(
                (reference[None, :, :] - candidates[:, None, :]).reshape(-1, 3),
                box,
            ).reshape(len(candidates), len(reference), 3)
            closest = np.linalg.norm(gaps, axis=2).min(axis=1)
            free = np.flatnonzero(closest >= self.CAP_CLEARANCE[element])
            if len(free):
                return candidates[free[0]]
            crowded += 1
            return candidates[np.argmax(closest)]

        for n, (atom_id, cap, toward, length) in enumerate(plan):
            anchor = anchors[n]
            others = [o for o in near[n] if o != anchor]
            origin = positions[anchor]
            if cap == "H":
                position = place("H", origin, toward, oh_length, others)
                caps.append(("H", position, atom_id))
                bonds.append((atom_id, -len(caps)))
                continue
            oxygen = place("O", origin, toward, length, others)
            caps.append(("O", oxygen, atom_id))
            bonds.append((atom_id, -len(caps)))
            oxygen_cap = -len(caps)
            # The H continues away from the anchor: Si-O-H is no line, but
            # the relaxation bends it.
            axis = (oxygen - origin) / np.linalg.norm(oxygen - origin)
            hydrogen = place("H", oxygen, axis, oh_length, others, own=len(caps) - 1)
            caps.append(("H", hydrogen, oxygen_cap))
            bonds.append((oxygen_cap, -len(caps)))

        return caps, bonds, crowded

    def _repair(
        self,
        plan: list[tuple],
        templates: tuple[dict, dict],
        explicit: tuple[set, set],
        bonds_dict: dict | None,
        oh_length: float,
    ) -> dict:
        """Add the caps, type them, and give them their explicit topology."""
        from ..core._format import canonical_ff_type as canonical

        system = self.system
        charge_before = float(system.atoms["charge"].sum())
        caps, cap_bonds, crowded = self._place_caps(plan, oh_length)

        first_id = int(system.atoms.index.max()) + 1
        cap_ids = [first_id + n for n in range(len(caps))]
        if caps:
            for element in {c[0] for c in caps}:
                if element not in system._masses:
                    system._masses[element] = MASSES_DICT[element]
            system.add_atoms([c[0] for c in caps], np.array([c[1] for c in caps]))

        def atom_id(ref: int) -> int:
            return cap_ids[-1 - ref] if ref < 0 else ref

        new_bonds = [(atom_id(a), atom_id(b)) for a, b in cap_bonds]

        # Type and charge the caps, and the atoms around them whose
        # environment the repair changed (an Ob that took an H, the H of
        # an Osih that became a water), from the atom of the intact
        # system in the same environment.
        self._contact_cache = {}
        first, second = self._environment(bonds_dict)
        by_first, by_second, before = templates
        ids = system.atoms.index.to_numpy()
        row = {i: n for n, i in enumerate(ids)}
        contacts = self._covalent_contacts(bonds_dict)
        partners: dict[int, list[int]] = {}
        for a, b in zip(ids[contacts["i"]], ids[contacts["j"]]):
            partners.setdefault(int(a), []).append(int(b))
            partners.setdefault(int(b), []).append(int(a))
        touched = {a for bond in new_bonds for a in bond}
        touched |= {p for a in list(touched) for p in partners.get(a, ())}
        retyped = untyped = 0
        relabel = set()
        for atom in sorted(touched):
            n = row[atom]
            is_cap = atom >= first_id
            if not is_cap and before.get(atom) == second[n]:
                continue
            match = by_second.get(second[n]) or by_first.get(first[n])
            old_type = system.atoms.at[atom, "type"]
            if match is None:
                untyped += 1
                atype = first[n][0]
                charge = 0.0 if is_cap else system.atoms.at[atom, "charge"]
                if atype not in system._masses:
                    system._masses[atype] = MASSES_DICT.get(atype, 1.0)
            else:
                atype, charge = match
            if not is_cap and atype != old_type:
                retyped += 1
                relabel.add(atom)
            system.atoms.at[atom, "type"] = atype
            system.atoms.at[atom, "charge"] = charge

        # The explicit connections of a retyped atom name its old type
        types = system.atoms["type"].astype(str)
        for name, kind, n_atoms in (("bonds", "bond", 2), ("angles", "angle", 3)):
            frame = getattr(system, name)
            if frame is None or frame.empty or not relabel:
                continue
            columns = [f"atom_{i + 1}" for i in range(n_atoms)]
            stale = frame[columns].isin(relabel).any(axis=1)
            known = explicit[0 if kind == "bond" else 1]
            frame.loc[stale, "type"] = [
                known.get(key, key)
                for key in (
                    canonical(tuple(types[a] for a in atoms), kind)
                    for atoms in frame.loc[stale, columns].itertuples(index=False)
                )
            ]

        bonds_added, angles_added = self._extend_topology(
            new_bonds, explicit, canonical
        )
        self._contact_cache = {}
        self._molecule_cache = None

        return {
            "broken": len(plan) // 2,
            "capped": len(caps),
            "crowded": crowded,
            "retyped": retyped,
            "untyped": untyped,
            "bonds_added": bonds_added,
            "angles_added": angles_added,
            "charge_change": float(system.atoms["charge"].sum()) - charge_before,
        }

    def _explicit_labels(self) -> tuple[dict, dict]:
        """The explicit bond and angle types of the system, by atom types.

        Keyed on the canonical form of the atom types they join, with the
        label the system itself uses as value: a builder may have written
        "O1-H1" where the canonical form reads "H1-O1".
        """
        from ..core._format import canonical_ff_type

        types = self.system.atoms["type"].astype(str)
        labels = []
        for name, kind, n_atoms in (("bonds", "bond", 2), ("angles", "angle", 3)):
            frame = getattr(self.system, name)
            table: dict[str, str] = {}
            if frame is not None and not frame.empty:
                columns = [f"atom_{i + 1}" for i in range(n_atoms)]
                for label, atoms in zip(
                    frame["type"], frame[columns].itertuples(index=False)
                ):
                    key = canonical_ff_type(tuple(types[a] for a in atoms), kind)
                    if isinstance(label, str):
                        table.setdefault(key, label)
            labels.append(table)
        return labels[0], labels[1]

    def _extend_topology(
        self, new_bonds: list[tuple[int, int]], explicit: tuple[set, set], canonical
    ) -> tuple[int, int]:
        """Give the new bonds the explicit bonds and angles the system has.

        A new pair gets an explicit bond when the intact system had one
        between the same two types (Osih-Hsi, Oh-Hh), and a new triple an
        angle likewise, with the system's label; pairs the force field
        leaves implicit (Si-O under ClayFF) stay so.
        """
        bond_types, angle_types = explicit
        system = self.system
        types = system.atoms["type"].astype(str)

        added_bonds = []
        for a, b in new_bonds:
            label = bond_types.get(canonical((types[a], types[b]), "bond"))
            if label is not None:
                added_bonds.append({"type": label, "atom_1": a, "atom_2": b})
        if added_bonds:
            frame = pd.DataFrame(added_bonds)
            system.bonds = (
                frame
                if system.bonds is None or system.bonds.empty
                else pd.concat([system.bonds, frame], ignore_index=True)
            )

        added_angles = []
        if added_bonds and angle_types:
            partners: dict[int, set] = {}
            for a, b in zip(system.bonds["atom_1"], system.bonds["atom_2"]):
                partners.setdefault(int(a), set()).add(int(b))
                partners.setdefault(int(b), set()).add(int(a))
            existing = set()
            if system.angles is not None and not system.angles.empty:
                existing = {
                    (int(a), int(c), int(e))
                    for a, c, e in zip(
                        system.angles["atom_1"],
                        system.angles["atom_2"],
                        system.angles["atom_3"],
                    )
                }
            for bond in added_bonds:
                a, b = bond["atom_1"], bond["atom_2"]
                for centre, end in ((a, b), (b, a)):
                    for other in partners.get(centre, ()):
                        if other == end:
                            continue
                        triple = (min(other, end), centre, max(other, end))
                        if triple in existing:
                            continue
                        label = angle_types.get(
                            canonical(
                                (types[triple[0]], types[centre], types[triple[2]]),
                                "angle",
                            )
                        )
                        if label is not None:
                            existing.add(triple)
                            added_angles.append(
                                {
                                    "type": label,
                                    "atom_1": triple[0],
                                    "atom_2": centre,
                                    "atom_3": triple[2],
                                }
                            )
        if added_angles:
            frame = pd.DataFrame(added_angles)
            system.angles = (
                frame
                if system.angles is None or system.angles.empty
                else pd.concat([system.angles, frame], ignore_index=True)
            )
        return len(added_bonds), len(added_angles)

    def split(
        self,
        repair: bool | None = None,
        bonds_dict: dict | None = None,
        oh_length: float = 1.0,
    ) -> AtomicSystem:
        """
        Execute the split and return the resulting AtomicSystem.

        Parameters
        ----------
        repair : bool, optional
            Rebuild the severed bonds (see the class description). Defaults
            to the splitter's own ``repair``; given here, it replaces it.
        bonds_dict : dict, optional
            Cutoff table used for that detection, see
            :data:`DEFAULT_BOND_CUTOFFS`.
        oh_length : float, default=1.0
            O-H distance in Å used when capping.

        Returns
        -------
        AtomicSystem
            The split system, optionally repaired and filled with solution:
            a new system, the one given to the Splitter is left unchanged.

        Raises
        ------
        ValueError
            With ``repair``, if the cut severs a bond without an O, which
            no hydrolysis rebuilds. Nothing has moved yet.

        Notes
        -----
        After a repair, ``repair_report`` gives ``broken`` (bonds
        rebuilt), ``capped`` (atoms added), ``crowded`` (caps that found
        no spot clear of other atoms), ``retyped`` (atoms of the system
        whose type changed, an Ob turned Osih), ``untyped`` (atoms with no
        model in the system: typed by element, a new one uncharged),
        ``bonds_added``, ``angles_added`` and ``charge_change``. The
        caps take their charges from the system, which leaves it charged
        in general: ``neutralize_charge`` afterwards.

        ``split_report`` compares, pair by pair, the bonds (explicit and
        ``bonds_dict``) and the ionic contacts (``ionic_dict``) of the
        intact geometry with the same pairs once split: any pair whose
        minimum-image distance changed is listed as broken.
        """
        if repair is not None:
            self.repair = repair
        self._resolve_coordinate(None)
        # Builders return a new system: the one given to the Splitter used
        # to be split in place and returned, so it came back opened (a
        # second cut of "the same" system cut the opened one) and shared its
        # charges with the result.
        self.system = self.system.copy()
        original_positions = self.system.atoms[["x", "y", "z"]].copy()
        old_box = np.array(self.system.box, dtype=np.float64)
        n_explicit = 0 if self.system.bonds is None else len(self.system.bonds)

        # Detection has to happen on the intact geometry, and so does the
        # reading of the types, charges and explicit topology to copy.
        if self.repair:
            plan = self._plan_repair(bonds_dict)
            templates = self._templates(bonds_dict)
            explicit = self._explicit_labels()
        bonds = self._covalent_contacts(bonds_dict)
        ionic = self._contacts(self._resolve_ionic())

        universe = self.system.to_mda()
        self._move_fragments(universe)
        self.system.atoms[["x", "y", "z"]] = universe.atoms.positions

        # The box must already reflect the gap before `_remove_crossing_bonds`
        # runs: it needs the *new* box to correctly apply the minimum-image
        # convention on the post-move geometry.
        self._update_box()

        self._remove_crossing_bonds(original_positions, old_box)
        self.split_report = self._compare_contacts(
            original_positions, old_box, bonds, ionic
        )
        self.split_report["explicit_bonds_removed"] = n_explicit - (
            0 if self.system.bonds is None else len(self.system.bonds)
        )
        if self.keep_molecules and self.split_report["bonds_broken"]:
            warnings.warn(
                f"{len(self.split_report['bonds_broken'])} bonds broken despite "
                "keep_molecules=True: see `split_report`."
            )

        # Positions (and, below, atoms) change: drop what was measured
        self._contact_cache = {}
        self._molecule_cache = None

        if self.repair:
            self.repair_report = self._repair(
                plan, templates, explicit, bonds_dict, oh_length
            )

        if self._has_solution:
            self._insert_solution()

        if hasattr(self.system, "metadata"):
            self.system.metadata["split_info"] = {
                "axis": self.axis,
                "axis_name": self._axis_name,
                "coordinate": self.coordinate,
                "gap_size": self.gap_size,
                "has_solution": self._has_solution,
                "solution_padding": self._padding if self._has_solution else None,
                "solution_thickness": self._solution_thickness,
                "solution_vacuum": self._vacuum if self._has_solution else None,
            }

        return self.system

    def _compare_contacts(
        self, original_positions: pd.DataFrame, old_box: np.ndarray, bonds: dict, ionic: dict
    ) -> dict:
        """Measure, on the split geometry, which of the given pairs broke."""
        from MDAnalysis.lib.distances import calc_bonds

        new_positions = self.system.atoms[["x", "y", "z"]].to_numpy(dtype=np.float64)
        old_positions = original_positions.to_numpy(dtype=np.float64)
        new_box = np.asarray(self.system.box, dtype=np.float64)
        ids = self.system.atoms.index.to_numpy()

        def broken_pairs(contacts: dict) -> list[tuple[int, int]]:
            i, j = contacts["i"], contacts["j"]
            if len(i) == 0:
                return []
            before = calc_bonds(old_positions[i], old_positions[j], box=old_box)
            after = calc_bonds(new_positions[i], new_positions[j], box=new_box)
            changed = np.abs(after - before) > 1e-3
            return [(int(ids[a]), int(ids[b])) for a, b in zip(i[changed], j[changed])]

        ionic_broken = broken_pairs(ionic)
        return {
            "coordinate": self.coordinate,
            "bonds": len(bonds["i"]),
            "bonds_broken": broken_pairs(bonds),
            "ionic_contacts": len(ionic["i"]),
            "ionic_broken": len(ionic_broken),
            "ionic_broken_pairs": ionic_broken,
        }

    def build(self, **kwargs) -> AtomicSystem:
        """Alias for :meth:`split`, to satisfy the ``BaseBuilder`` interface."""
        return self.split(**kwargs)
