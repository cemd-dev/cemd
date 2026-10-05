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

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
import pyvista as pv
from PySide6 import QtCore, QtWidgets

from .plotter_widget import AtomicPlotter
from .ui.atom_table import AtomTable

if TYPE_CHECKING:
    from PySide6.QtGui import QKeyEvent

    from ..core.atomic_system import AtomicSystem
    from .main_window import AtomViewerGUI


class StructureTabWidget(QtWidgets.QWidget):
    def __init__(
        self, system: AtomicSystem, parent_gui: AtomViewerGUI, file_path: str = None
    ) -> None:
        super().__init__()
        self.system = system
        self.file_path = file_path
        self.parent_gui = parent_gui
        self.df_visible = pd.DataFrame()
        self.selected_real_indices = []

        self.setup_ui()

    def setup_ui(self) -> None:
        """Initializes the user interface components, splitter, plotter, and atom table."""
        layout = QtWidgets.QVBoxLayout(self)
        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)

        # move the creation of the plotter and the table here
        self.plotter = AtomicPlotter(self, config=self.parent_gui.global_config)
        self.plotter.enable_point_picking(
            callback=None, show_point=False, left_clicking=True
        )
        self.plotter.render_window.GetInteractor().AddObserver(
            "LeftButtonPressEvent", self._force_pick_callback
        )

        # Add the keyboard listener (for the Del key)
        self.plotter.render_window.GetInteractor().AddObserver(
            "KeyPressEvent", self._handle_plotter_key
        )

        self.table = AtomTable()
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)

        self.splitter.addWidget(self.plotter)
        self.splitter.addWidget(self.table)
        self.splitter.setSizes([700, 300])
        layout.addWidget(self.splitter)

    def _handle_plotter_key(self, obj: Any, event: str) -> None:
        """Intercepts keyboard keys coming from the 3D Plotter"""
        key = self.plotter.interactor.GetKeySym()
        if key in ["Delete", "BackSpace"]:
            self.delete_selected()

    def _force_pick_callback(self, obj: Any, event: str) -> None:
        """Calculates the 3D coordinates of a mouse click and performs a proximity search to select the nearest atom."""
        click_pos = obj.GetEventPosition()

        # Use the Picker (Point or Cell)
        picker = obj.GetPicker()
        # slightly increase the native tolerance of VTK
        picker.SetTolerance(0.005)
        picker.Pick(click_pos[0], click_pos[1], 0, self.plotter.renderer)

        picked_point = picker.GetPickPosition()

        # If the point is (0,0,0) and VTK found nothing, try a manual Ray Cast
        # If any(...) is False, the click did not even touch the rendering plane
        if not any(v != 0 for v in picked_point):
            return

        # IMPROVED PROXIMITY SEARCH
        if not self.df_visible.empty:
            xyz = self.df_visible[["x", "y", "z"]].values

            # Calculation of distances between the clicked point in 3D space and all atoms
            distances = np.linalg.norm(xyz - picked_point, axis=1)
            min_dist = np.min(distances)
            closest_idx = np.argmin(distances)

            # Dynamic tolerance: accept a click if are close to an atom.
            # 2.5 Å is a good value for "click next"
            if min_dist < 2.5:
                closest_xyz = xyz[closest_idx]
                # call on_pick with the coordinates of the real atom
                self.on_pick(closest_xyz)
            else:
                # If are too far away, empty the selection (click in the empty space voluntarily)
                mods = QtWidgets.QApplication.keyboardModifiers()
                if not mods == QtCore.Qt.ControlModifier:
                    self.selected_real_indices = []
                    self.update_selection_plot()
                    self.sync_table_selection()

    def on_pick(self, point: np.ndarray | Sequence[float]) -> None:
        """Manages the selection and updates the display"""
        if self.df_visible.empty:
            return

        xyz = self.df_visible[["x", "y", "z"]].values
        distances = np.linalg.norm(xyz - point, axis=1)

        idx = np.argmin(distances)
        real_idx = int(self.df_visible.index[idx])

        # Selection logic (CTRL or simple selection)
        mods = QtWidgets.QApplication.keyboardModifiers()
        if mods == QtCore.Qt.ControlModifier:
            if real_idx in self.selected_real_indices:
                self.selected_real_indices.remove(real_idx)
            else:
                self.selected_real_indices.append(real_idx)
        else:
            self.selected_real_indices = [real_idx]

        # Update of the TABLE (3D Sense -> Table)
        # block the signals to prevent the table from sending a refresh order to 3D
        self.table.blockSignals(True)
        self.sync_table_selection()
        self.table.blockSignals(False)
        self.table.setFocus()

        # Updated 3D HIGHLIGHT (Variable Direction -> Plotter)
        self.update_selection_plot()

        if hasattr(self.parent_gui, "update_protonate_state"):
            self.parent_gui.update_protonate_state()

        # Display in status bar
        atom = self.system.atoms.loc[real_idx]
        self.parent_gui.statusBar().showMessage(
            f"Atom {real_idx} | Type: {atom['type']} | Pos: ({atom['x']:.2f}, {atom['y']:.2f}, {atom['z']:.2f})",
            3000,
        )

    def on_table_selection_changed(
        self, selected: QtCore.QItemSelection, deselected: QtCore.QItemSelection
    ) -> None:
        """Synchronizes the 3D plotter selection highlight when rows are selected or deselected in the data table."""
        model = self.table.model()
        if model is None:
            return

        rows = self.table.selectionModel().selectedRows()
        # Tab Truth List Update
        self.selected_real_indices = [
            int(model.system_obj.atoms.index[r.row()]) for r in rows
        ]

        # Visual update (the pink halo)
        self.update_selection_plot()

        # ---CRUCIAL CALL (NOW VIA PARENT) ---
        # tell the MainWindow: "Hey, my selection has changed, check if you need to activate the button"
        if hasattr(self.parent_gui, "update_protonate_state"):
            self.parent_gui.update_protonate_state()

        # Visual refresh of the table
        self.table.viewport().update()

    def update_selection_plot(self) -> None:
        """Updates the 3D visualization by adding or removing highlight halos around selected atoms."""
        # retrieve all the actor names present in the renderer
        for actor_name in list(self.plotter.renderer.actors.keys()):
            if actor_name.startswith("sel_") or actor_name == "selection_highlight":
                self.plotter.remove_actor(actor_name)

        if hasattr(self, "topo_win") and self.topo_win.isVisible():
            self.topo_win.highlight_atoms_in_tables(self.selected_real_indices)

        if not self.selected_real_indices:
            self.plotter.render()
            return

        global_scale = self.parent_gui.filter_panel.get_scale_value()

        # Collect all positions at once
        positions = []
        point_sizes = []

        for idx in self.selected_real_indices:
            atom = self.system.atoms.loc[idx]
            atype = str(atom["type"])
            positions.append([atom["x"], atom["y"], atom["z"]])
            base_radius = self.plotter.radius_map.get(atype, 1.5)
            point_sizes.append(base_radius * 20 * global_scale)

        # A single mesh for all selected atoms
        mesh = pv.PolyData(np.array(positions))

        # PyVista does not support different point_sizes per point
        # We take the maximum size as a compromise
        selection_size = max(point_sizes)

        self.plotter.add_mesh(
            mesh,
            color="pink",
            point_size=selection_size,
            render_points_as_spheres=True,
            name="selection_highlight",  # ← only one name
            reset_camera=False,
            opacity=0.9,
        )

        self.plotter.render()

    def sync_table_selection(self) -> None:
        """Synchronizes the table's visual selection and scroll position based on the current list of selected atom indices."""
        model = self.table.model()
        sel_model = self.table.selectionModel()
        if model is None or sel_model is None:
            return

        sel_model.blockSignals(True)
        sel_model.clearSelection()

        selected_set = set(self.selected_real_indices)
        if not selected_set:
            sel_model.blockSignals(False)
            return

        selection = QtCore.QItemSelection()
        first_index = None
        index_array = model.system_obj.atoms.index

        for row, real_id in enumerate(index_array):
            if real_id in selected_set:  # ← O(1) instead of O(n)
                left = model.index(row, 0)
                right = model.index(row, model.columnCount() - 1)
                selection.select(left, right)
                if first_index is None:
                    first_index = left

        sel_model.select(selection, QtCore.QItemSelectionModel.SelectionFlag.Select)

        # ---THE CRUCIAL ADDITION: Scroll to the atom ---
        if first_index:
            self.table.scrollTo(
                first_index, QtWidgets.QAbstractItemView.PositionAtCenter
            )

        sel_model.blockSignals(False)
        self.table.viewport().update()

    def delete_selected(self) -> None:
        """Removes the currently selected atoms from the atomic system and triggers a global UI synchronization."""
        # Retrieving selected rows
        selection = self.table.selectionModel().selectedRows()
        if not selection:
            return

        # Retrieving real IDs from the DataFrame
        model = self.table.model()
        indices_to_remove = [model.system_obj.atoms.index[s.row()] for s in selection]

        # Delete action in Data object
        self.system.remove_atoms(indices_to_remove)

        # THE FIX: call the function that really exists!
        # put full_rebuild=True because the number of atoms has changed,
        # you must therefore update the FilterPanel counters if necessary.
        self.parent_gui.sync_ui(full_rebuild=True)

    def keyPressEvent(self, event: QKeyEvent) -> None:
        """Handles keyboard shortcuts, specifically intercepting the Delete and Backspace keys for atom removal."""
        # Security: If are editing a cell, let the table do its work
        if self.table.state() == QtWidgets.QAbstractItemView.EditingState:
            super().keyPressEvent(event)
            return

        # intercept Del (and Backspace for Mac/Laptop keyboards)
        if event.key() in [QtCore.Qt.Key_Delete, QtCore.Qt.Key_Backspace]:
            if self.selected_real_indices:
                self.delete_selected()
            else:
                self.parent_gui.statusBar().showMessage(
                    "Aucun atome sélectionné.", 2000
                )
        else:
            # Very important for other keys (Escape, etc.) to work
            super().keyPressEvent(event)

    def on_data_changed(self, top_left, bottom_right=None):
        """Processes manual edits in the table cells to update atom coordinates or types within the underlying atomic system."""
        # Basic security
        if getattr(self.parent_gui, "_is_syncing", False) or self.system is None:
            return

        # CASE A: Forced refresh (exit BEFORE reading .row())
        if top_left is None:
            # request global synchronization (which will do the update_info)
            self.parent_gui.sync_ui(full_rebuild=True)
            return

        # CASE B: Modification of a specific cell
        try:
            # Now can read .row() and .column() safely
            row = top_left.row()
            col = top_left.column()
            model = self.table.model()

            # Retrieving the real ID via column 0 (hidden or not)
            real_idx = model.index(row, 0).data(QtCore.Qt.ItemDataRole.DisplayRole)

            cols_map = {1: "type", 2: "x", 3: "y", 4: "z", 5: "charge"}
            col_name = cols_map.get(col)
            if not col_name:
                return

            # ---UPDATE LOGIC ---
            if col_name in ["x", "y", "z"]:
                # `set_atom_position` takes a single atom index and the
                # full (x, y, z) position, not a per-axis keyword -- start
                # from the atom's current position and only change the
                # edited axis.
                new_val = model.data(top_left, QtCore.Qt.ItemDataRole.EditRole)
                position = self.system.atoms.loc[real_idx, ["x", "y", "z"]].to_numpy(
                    dtype=float, copy=True
                )
                position["xyz".index(col_name)] = float(new_val)
                self.system.set_atom_position(real_idx, position)
                QtCore.QTimer.singleShot(
                    10, lambda: self.parent_gui.sync_ui(full_rebuild=False)
                )

            elif col_name == "type":
                new_val = model.data(top_left, QtCore.Qt.ItemDataRole.EditRole)
                new_type = str(new_val)
                self.system.set_type2atoms([real_idx], new_type)
                # Gray by default [2026-03-08]
                if new_type not in self.plotter.color_map:
                    self.plotter.color_map[new_type] = "gray"
                QtCore.QTimer.singleShot(
                    10, lambda: self.parent_gui.sync_ui(full_rebuild=True)
                )

            # NOTE: For 'charge', do nothing!
            # The model (setData) has already written to the DataFrame.

            # ---DELAYED UI SYNC ---
            # call sync_ui which contains the single call to update_info

        except Exception as e:
            print(f"Erreur d'édition dans l'onglet : {e}")

    def refresh_tab_view(
        self,
        full_rebuild: bool = False,
        reset_camera: bool = False,
        refresh_bonds: bool = True,
    ) -> None:
        """Manages 3D drawing and synchronization between the UI and the Plotter."""
        if self.system is None or self.system.atoms.empty:
            return

        # # ---UI TABLE UPDATE ---
        if full_rebuild:
            self.table.update_data(self.system)
            try:
                self.table.model().dataChanged.disconnect()
            except (RuntimeError, TypeError):  # nothing connected yet
                pass
            self.table.model().dataChanged.connect(self.on_data_changed)
            try:
                self.table.selectionModel().selectionChanged.disconnect()
            except (RuntimeError, TypeError):  # nothing connected yet
                pass
            self.table.selectionModel().selectionChanged.connect(
                self.on_table_selection_changed
            )

        # --- PLOTTER SETTINGS UPDATE ---
        f_panel = self.parent_gui.filter_panel
        f_settings = f_panel.get_settings()

        self.plotter.global_scale = f_settings.get("global_scale", 100) / 100.0
        if "radii" in f_settings:
            self.plotter.radius_map.update(f_settings["radii"])
        if "colors" in f_settings:
            self.plotter.color_map.update(f_settings["colors"])

        b_panel = self.parent_gui.bond_manager
        bond_settings = b_panel.get_settings()
        if bond_settings:
            self.plotter.global_bond_radius = b_panel.get_bond_radius()
            self.plotter.bond_map.update(
                {k: v for k, v in bond_settings.items() if k != "global_bond_radius"}
            )

        # Sync global config
        self.parent_gui.global_config = self.plotter.get_current_config_dict()
        self.plotter.color_map.update(
            self.parent_gui._config_cache.get("color_map", {})
        )
        self.plotter.radius_map.update(
            self.parent_gui._config_cache.get("radius_map", {})
        )

        # Visible atoms filter
        sel_types = [str(t) for t, cb in f_panel.checkboxes.items() if cb.isChecked()]
        if not sel_types:
            sel_types = self.system.atoms["type"].astype(str).unique().tolist()
        self.df_visible = self.system.atoms[
            self.system.atoms["type"].astype(str).isin(sel_types)
        ].copy()

        # --- 3D RENDERING ---
        if full_rebuild:
            # Draw every type, hidden ones included: an actor that a
            # rebuild never creates cannot be shown again later by simply
            # toggling its checkbox back on.
            self.plotter.clear()
            self.plotter.draw_atoms(self.system.atoms)
            self.plotter.draw_box(self.system)
        else:
            self.plotter.update_atom_sizes()
            self.plotter.update_atom_colors()

        self.plotter.update_atom_visibility(set(sel_types))

        if refresh_bonds and not self.df_visible.empty:
            self.plotter.update_bonds_only(self.df_visible, bond_settings)

        if reset_camera:
            self.plotter.reset_camera()

        self.plotter.render()
