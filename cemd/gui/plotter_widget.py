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

from typing import TYPE_CHECKING, Any

import numpy as np
import pyvista as pv
import vtk
from pyvistaqt import QtInteractor
from scipy.spatial import cKDTree

from .display_settings import BG_CYCLE, DisplaySettings

if TYPE_CHECKING:
    import pandas as pd

    from cemd.core.atomic_system import AtomicSystem

    from .main_window import AtomViewerGUI


class AtomicPlotter(QtInteractor):
    #: Sides of the extruded tube used to draw a bond. Eight is smooth
    #: enough at any zoom; each side costs cells on every rendered frame.
    BOND_TUBE_SIDES = 8

    #: Above this many bonds, tubes are replaced by hardware-drawn lines.
    #: Tubing 2400 bonds already produces ~48k cells, and the geometry is
    #: rebuilt whenever the bond settings change.
    MAX_TUBED_BONDS = 8000

    #: Atoms drawn while the view is rotated, panned or zoomed. Each point
    #: sprite costs on every frame -- 500 000 atoms take ~140 ms a frame on
    #: an integrated GPU, 50 000 ~30 ms -- so a larger system is drawn from
    #: a subset of its atoms, and without its bonds, while it moves, and
    #: whole again once it stops.
    LOD_ATOM_BUDGET = 100_000

    #: Types with no more atoms than this are drawn whole even while the
    #: view moves, so that ions and other minor species do not blink out.
    LOD_KEEP_TYPES_UP_TO = 1000

    bg_cycle = BG_CYCLE

    def __init__(
        self, parent: AtomViewerGUI = None, settings: DisplaySettings | None = None
    ):
        super().__init__(parent)
        #: Shared with every other tab's plotter, not copied: a change made
        #: through one is seen by all, with nothing to keep in step.
        self.settings = settings if settings is not None else DisplaySettings()

        self.update_background_style()
        # self.render_window.SetMultiSamples(4)
        # self.enable_anti_aliasing('msaa')
        self.enable_lightkit()
        self.enable_eye_dome_lighting()

        self.add_axes()

        #: The actors hidden for the moving view, to show again once it
        #: stops; see `_on_render_start`.
        self._hidden_while_moving: list = []
        self.render_window.AddObserver("StartEvent", self._on_render_start)

    # The drawing code reads these names; they all lead to the shared
    # settings.

    @property
    def color_map(self) -> dict[str, str]:
        return self.settings.color_map

    @property
    def radius_map(self) -> dict[str, float]:
        return self.settings.radius_map

    @property
    def bond_map(self) -> dict[str, float]:
        return self.settings.bond_map

    @property
    def global_scale(self) -> float:
        return self.settings.global_scale

    @property
    def global_bond_radius(self) -> float:
        return self.settings.global_bond_radius

    def hex_to_rgb(self, h: str) -> list[float]:
        """Converts colors from JSON to 3D format [0, 1]"""
        if h == "gray":
            return [0.5, 0.5, 0.5]
        try:
            h = h.lstrip("#")
            return [int(h[i : i + 2], 16) / 255.0 for i in (0, 2, 4)]
        except (ValueError, AttributeError):
            return [0.5, 0.5, 0.5]

    def add_axes(self) -> None:
        """Adds axes with white titles"""

        # Create axes actor
        axes_actor = vtk.vtkAxesActor()

        # Configure axes geometry
        axes_actor.SetTotalLength(3, 3, 3)
        axes_actor.SetShaftTypeToCylinder()
        axes_actor.SetCylinderRadius(0.1)
        axes_actor.SetConeRadius(1)

        # Configure labels
        for axis in [
            axes_actor.GetXAxisCaptionActor2D(),
            axes_actor.GetYAxisCaptionActor2D(),
            axes_actor.GetZAxisCaptionActor2D(),
        ]:
            prop = axis.GetCaptionTextProperty()
            prop.SetColor(1, 1, 1)

        # Create orientation widget
        self.axes_widget = vtk.vtkOrientationMarkerWidget()
        self.axes_widget.SetOrientationMarker(axes_actor)
        self.axes_widget.SetInteractor(self.interactor)
        self.axes_widget.SetViewport(0.0, 0.0, 0.3, 0.3)
        self.axes_widget.SetEnabled(1)
        self.axes_widget.InteractiveOff()

    def draw_atoms(self, df: pd.DataFrame) -> None:
        """Group atoms by type and render them as point-based spheres."""
        for actor in list(self.renderer.actors.keys()):
            if actor.startswith(("atoms_", "atom_single_", "lod_atoms_")):
                self.remove_actor(actor)

        if df is None or df.empty:
            return

        # Every stride-th atom of each type makes the moving view.
        stride = -(-len(df) // self.LOD_ATOM_BUDGET)

        for atype, group_df in df.groupby("type"):
            xyz = group_df[["x", "y", "z"]].values
            mesh = pv.PolyData(xyz)
            self._add_atom_mesh(mesh, atype, self.global_scale, name=f"atoms_{atype}")

            if stride > 1 and len(xyz) > self.LOD_KEEP_TYPES_UP_TO:
                lod_actor = self._add_atom_mesh(
                    pv.PolyData(xyz[::stride]),
                    atype,
                    self.global_scale,
                    name=f"lod_atoms_{atype}",
                    pickable=False,
                )
                lod_actor.visibility = False

    def _on_render_start(self, render_window: Any, event: str) -> None:
        """Swap in the light scene while the view moves, the full one after.

        VTK's interactor styles raise the render window's desired update
        rate for the length of a rotation, pan or zoom, and set it back to
        the still rate for the frame drawn when it ends. Read at the start
        of each render, that rate says which scene the frame needs, whatever
        the interactor style.
        """
        interactor = render_window.GetInteractor()
        if interactor is None:
            return
        moving = render_window.GetDesiredUpdateRate() > interactor.GetStillUpdateRate()

        if not moving:
            # An actor removed since (a redraw) is gone from the scene, and
            # showing it again does nothing.
            for actor in self._hidden_while_moving:
                actor.visibility = True
            self._hidden_while_moving = []
            for name, actor in self.renderer.actors.items():
                if name.startswith("lod_atoms_"):
                    actor.visibility = False
            return

        if self._hidden_while_moving:
            return
        actors = self.renderer.actors
        for name, actor in actors.items():
            if not name.startswith("lod_atoms_"):
                continue
            full = actors.get(name[len("lod_") :])
            if full is not None and full.visibility:
                full.visibility = False
                actor.visibility = True
                self._hidden_while_moving.append(full)
        bonds = actors.get("bonds")
        if self._hidden_while_moving and bonds is not None and bonds.visibility:
            bonds.visibility = False
            self._hidden_while_moving.append(bonds)

    def _add_atom_mesh(
        self,
        mesh: pv.PolyData,
        atype: str | int,
        scale: float,
        name: str,
        pickable: bool = True,
    ) -> pv.Actor:
        """Helper to style and add an atom mesh to the renderer."""
        atype_str = str(atype).strip()

        # get the color from the object map (updated by logic.py)
        color_hex = self.color_map.get(atype_str, "gray")

        # make sure that the radius is there too
        base_radius = self.radius_map.get(atype_str, 1.5)
        display_size = base_radius * 10 * scale

        # sphere_source = pv.Sphere(radius=display_size, theta_resolution=16, phi_resolution=16)
        # styled_mesh = mesh.glyph(geom=sphere_source)

        if name in self.renderer.actors:
            self.remove_actor(name)

        # use color=color_hex directly.
        # PyVista accepts names 'white', 'red' or HEX codes '#FFFFFF'
        return self.add_mesh(
            mesh,
            color=color_hex,
            render_points_as_spheres=True,
            point_size=display_size,
            name=name,
            pickable=pickable,
            reset_camera=False,
            lighting=True,
            style="points",
            ambient=0.5,
        )

    def add_bonds_smart(
        self,
        xyz: np.ndarray,
        types: np.ndarray | list,
        exceptions: dict[str, Any] | None = None,
    ) -> None:
        """Calculate and render bonds between atoms based on distance rules."""
        if not exceptions or len(xyz) < 2:
            return

        if "bonds" in self.renderer.actors:
            self.remove_actor("bonds")

        bond_thickness = exceptions.get("global_bond_radius", 0.1)

        str_types = np.asarray([str(t).strip() for t in types])
        max_search = max(exceptions.values())
        tree = cKDTree(xyz)

        # output_type="ndarray" avoids materializing a set of tuples, and
        # lets everything below stay in numpy.
        pairs = tree.query_pairs(max_search, output_type="ndarray")
        if len(pairs) == 0:
            return

        first, second = pairs[:, 0], pairs[:, 1]

        # Cutoffs as a type-by-type matrix, so the per-pair test becomes a
        # single fancy-indexed comparison instead of a Python loop building
        # a "A-B" key and calling np.linalg.norm once per pair. The keys are
        # composed exactly as before, which keeps type names containing a
        # dash unambiguous.
        uniq_types, codes = np.unique(str_types, return_inverse=True)
        n_types = len(uniq_types)
        cutoffs = np.zeros((n_types, n_types))
        for a_index, a_type in enumerate(uniq_types):
            for b_index, b_type in enumerate(uniq_types):
                key = "-".join(sorted([a_type, b_type]))
                if key in exceptions:
                    cutoffs[a_index, b_index] = exceptions[key]

        limits = cutoffs[codes[first], codes[second]]
        distances = np.linalg.norm(xyz[first] - xyz[second], axis=1)

        # A zero cutoff means "these two types are not bonded", which is
        # distinct from a pair that is simply too far apart.
        keep = (limits > 0) & (distances <= limits)
        first, second = first[keep], second[keep]
        if len(first) == 0:
            return

        # Two half-bonds per pair, so each half can carry its own atom's
        # color: point 0 -> midpoint, midpoint -> point 1.
        start, end = xyz[first], xyz[second]
        middle = (start + end) / 2.0

        points = np.empty((4 * len(first), 3), dtype=float)
        points[0::4], points[1::4] = start, middle
        points[2::4], points[3::4] = middle, end

        offsets = np.arange(len(first)) * 4
        pair_header = np.full(len(first), 2)
        lines = np.column_stack(
            [
                pair_header,
                offsets,
                offsets + 1,
                pair_header,
                offsets + 2,
                offsets + 3,
            ]
        ).ravel()

        rgb_table = np.array(
            [self.hex_to_rgb(self.color_map.get(t, "gray")) for t in uniq_types]
        )
        colors = np.empty((4 * len(first), 3))
        colors[0::4] = colors[1::4] = rgb_table[codes[first]]
        colors[2::4] = colors[3::4] = rgb_table[codes[second]]

        b_mesh = pv.PolyData(points)
        b_mesh.lines = lines
        b_mesh.point_data["colors"] = (colors * 255).astype(np.uint8)

        # Tubes are real geometry: at n_sides=8 a few thousand bonds turn
        # into tens of thousands of cells that VTK re-renders every frame.
        # Past a threshold the tube is dropped for plain lines, which look
        # near-identical at the zoom level a large system is viewed at and
        # cost roughly a third of the cells.
        if len(first) <= self.MAX_TUBED_BONDS:
            mesh = b_mesh.tube(radius=bond_thickness, n_sides=self.BOND_TUBE_SIDES)
            self.add_mesh(
                mesh,
                scalars="colors",
                rgb=True,
                name="bonds",
                reset_camera=False,
                pickable=False,
            )
        else:
            self.add_mesh(
                b_mesh,
                scalars="colors",
                rgb=True,
                name="bonds",
                reset_camera=False,
                pickable=False,
                line_width=max(1, int(bond_thickness * 30)),
                render_lines_as_tubes=True,
            )

    @staticmethod
    def _box_corners(system: AtomicSystem) -> np.ndarray:
        """The eight corners of the cell, in VTK hexahedron order."""
        v = np.array(
            [
                [0, 0, 0],
                [1, 0, 0],
                [1, 1, 0],
                [0, 1, 0],
                [0, 0, 1],
                [1, 0, 1],
                [1, 1, 1],
                [0, 1, 1],
            ]
        )
        return (v @ np.array(system._box_vectors)) + np.array(
            [system._box_lmp[0][0], system._box_lmp[1][0], system._box_lmp[2][0]]
        )

    def draw_box(self, system: AtomicSystem) -> None:
        """Draws the simulation box with thickness and adaptive color"""
        if system is None or not hasattr(system, "_box_vectors"):
            return
        try:
            pts = self._box_corners(system)

            faces = np.hstack(
                [
                    [4, 0, 1, 2, 3],
                    [4, 4, 5, 6, 7],
                    [4, 0, 1, 5, 4],
                    [4, 1, 2, 6, 5],
                    [4, 2, 3, 7, 6],
                    [4, 3, 0, 4, 7],
                ]
            )

            # Dark edges on a light background, light ones on a dark one.
            box_color = "black" if self.settings.is_light_background else "white"

            mesh = pv.PolyData(pts, faces)
            # box_tubes = mesh.extract_all_edges().tube(radius=0.03)

            self.add_mesh(
                mesh,
                color=box_color,
                style="wireframe",
                line_width=2,  # Thicker for visibility
                opacity=0.5,  # A little more opaque to clearly see the edges
                pickable=False,
                name="box",  # CRUCIAL NAME for update
                reset_camera=False,
            )
        except Exception as e:
            print(f"Could not draw the box: {e}")

    def draw_cut_preview(
        self,
        system: AtomicSystem,
        normal: np.ndarray,
        coordinate: float,
        broken: np.ndarray | None = None,
        shifted: np.ndarray | None = None,
    ) -> None:
        """Show a split plane and the bonds it would sever.

        The plane is the set of points whose projection on the unit vector
        `normal` equals `coordinate` (the convention of `Splitter`),
        clipped to the cell. `broken` holds one (start, end) segment per
        severed bond, shape (n, 2, 3). `shifted` holds the positions of
        the atoms a cut that follows the molecules puts on the other side
        of the plane, haloed in orange.
        """
        self.clear_cut_preview()

        corners = self._box_corners(system)
        cell = pv.UnstructuredGrid(
            {pv.CellType.HEXAHEDRON: np.arange(8)[None]}, corners
        )
        normal = np.asarray(normal, dtype=float)
        # A plane lying on a face of the cell slices nothing: the best cut
        # often sits right on the periodic boundary, so nudge it inside.
        projections = corners @ normal
        low, high = projections.min(), projections.max()
        tolerance = 1e-6 * (high - low)
        plane = None
        if low - tolerance <= coordinate <= high + tolerance:
            coordinate = np.clip(coordinate, low + tolerance, high - tolerance)
            plane = cell.slice(normal=normal, origin=coordinate * normal)
        if plane is not None and plane.n_points:
            self.add_mesh(
                plane,
                color="#3b82f6",
                opacity=0.35,
                show_edges=True,
                edge_color="#3b82f6",
                name="cut_plane",
                pickable=False,
                reset_camera=False,
            )

        if shifted is not None and len(shifted):
            self.add_mesh(
                pv.PolyData(np.asarray(shifted, dtype=float)),
                color="orange",
                style="points",
                render_points_as_spheres=True,
                point_size=1.5 * 20 * self.global_scale,
                opacity=0.6,
                name="cut_shifted",
                pickable=False,
                reset_camera=False,
            )

        if broken is not None and len(broken):
            n_bonds = len(broken)
            mesh = pv.PolyData(np.asarray(broken, dtype=float).reshape(-1, 3))
            offsets = np.arange(n_bonds) * 2
            mesh.lines = np.column_stack(
                [np.full(n_bonds, 2), offsets, offsets + 1]
            ).ravel()
            # Well thicker than the regular bonds, and the atoms at both
            # ends haloed like a selection: the atom spheres hide a tube
            # of bond thickness.
            radius = max(0.25, 2.5 * self.global_bond_radius)
            self.add_mesh(
                mesh.copy(),
                color="red",
                style="points",
                render_points_as_spheres=True,
                point_size=1.5 * 20 * self.global_scale,
                opacity=0.6,
                name="cut_atoms",
                pickable=False,
                reset_camera=False,
            )
            if n_bonds <= self.MAX_TUBED_BONDS:
                self.add_mesh(
                    mesh.tube(radius=radius, n_sides=self.BOND_TUBE_SIDES),
                    color="red",
                    name="cut_bonds",
                    pickable=False,
                    reset_camera=False,
                )
            else:
                self.add_mesh(
                    mesh,
                    color="red",
                    name="cut_bonds",
                    pickable=False,
                    reset_camera=False,
                    line_width=max(2, int(radius * 30)),
                    render_lines_as_tubes=True,
                )

        self.render()

    def clear_cut_preview(self) -> None:
        """Remove what `draw_cut_preview` drew."""
        for name in ("cut_plane", "cut_bonds", "cut_atoms", "cut_shifted"):
            if name in self.renderer.actors:
                self.remove_actor(name)

    def update_bonds_only(
        self, df: pd.DataFrame, bond_settings: dict[str, Any]
    ) -> None:
        """Updates bindings only"""
        if df is None or df.empty:
            return

        xyz = df[["x", "y", "z"]].values
        atom_types = df["type"].values

        self.add_bonds_smart(xyz, atom_types, bond_settings)

    def update_background_style(self) -> None:
        """Changes the background and forces the box to change color"""
        self.set_background(self.settings.bg_color)

        # The same rule as the box in `draw_box`. This one used to check
        # for "white" alone, so on the light blue background the axis
        # labels turned white while the box was drawn black.
        contrast_color = (0, 0, 0) if self.settings.is_light_background else (1, 1, 1)

        # Axes Update
        if hasattr(self, "axes_widget"):
            axes_actor = self.axes_widget.GetOrientationMarker()
            for ax in [
                axes_actor.GetXAxisCaptionActor2D(),
                axes_actor.GetYAxisCaptionActor2D(),
                axes_actor.GetZAxisCaptionActor2D(),
            ]:
                ax.GetCaptionTextProperty().SetColor(contrast_color)

        # REAL box color update
        if "box" in self.renderer.actors:
            box_actor = self.renderer.actors["box"]

            # force the color via the mapper object (deeper in VTK)
            box_actor.prop.SetColor(contrast_color)

            # deactivate the scalar mode which could block the fixed color
            box_actor.mapper.scalar_visibility = False

            # reset the thickness just in case
            box_actor.prop.line_width = 2

        self.render()

    def update_atom_sizes(self) -> None:
        """Only updates the point_size of existing atom actors."""
        for atype in self.radius_map:
            for actor_name in (f"atoms_{atype}", f"lod_atoms_{atype}"):
                if actor_name in self.renderer.actors:
                    actor = self.renderer.actors[actor_name]
                    base_radius = self.radius_map.get(str(atype), 1.5)
                    actor.prop.point_size = base_radius * 10 * self.global_scale

    def update_atom_visibility(self, visible_types: set[str]) -> None:
        """Show/hide whole per-type atom actors without rebuilding geometry."""
        for actor_name, actor in self.renderer.actors.items():
            if actor_name.startswith("atoms_"):
                atype = actor_name[len("atoms_") :]
                actor.visibility = atype in visible_types

    def update_atom_colors(self) -> None:
        """Only updates the colors of existing atom actors."""
        for atype in self.color_map:
            for actor_name in (f"atoms_{atype}", f"lod_atoms_{atype}"):
                if actor_name in self.renderer.actors:
                    actor = self.renderer.actors[actor_name]
                    color_hex = self.color_map.get(str(atype), "gray")
                    actor.prop.color = color_hex
