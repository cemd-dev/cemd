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

import json
import os
import sys

# setdefault, not assignment: xcb is the right default for VTK, which
# renders poorly under Wayland, but forcing it made the module unusable
# anywhere without an X server -- a headless CI runner, or a user who has
# deliberately chosen another platform plugin.
os.environ.setdefault("QT_QPA_PLATFORM", "xcb")
os.environ.setdefault("QT_X11_NO_MITSHM", "1")
os.environ.setdefault("QT_LINUX_ACCESSIBILITY_ALWAYS_ON", "0")

from typing import Any

from PySide6 import QtCore, QtGui, QtWidgets
from PySide6Qlementine import QlementineStyle

from cemd.core.atomic_system import AtomicSystem
from cemd.gui import _userdata
from cemd.gui.logic.build import (
    on_protonate,
    open_add_droplet,
    open_add_liquid,
    open_add_structure,
    open_make_cash,
    open_make_glass,
    open_make_solution,
    open_make_surface,
    open_replicate,
    open_smiles_builder,
    open_split,
    open_translate_atoms,
)
from cemd.gui.logic.file_handler import open_file, save_file, save_file_as
from cemd.gui.tabs import StructureTabWidget
from cemd.gui.ui.analysis_view import RDFDialog, SilicateDialog
from cemd.gui.ui.cod import CODBrowserDialog
from cemd.gui.ui.gui_utils import get_icon
from cemd.gui.ui.managers import ConnectivityDialog, TypeManagerDialog
from cemd.gui.ui.panels import BondManagerPanel, FilterPanel, SystemSummaryPanel
from cemd.gui.ui.pubchem import PubChemBrowserDialog


def deep_update(base_dict: dict, update_with: dict) -> None:
    """Recursively merge dictionaries (for color_map, etc.)"""
    for key, value in update_with.items():
        if isinstance(value, dict) and key in base_dict:
            deep_update(base_dict[key], value)
        else:
            base_dict[key] = value


def get_config_diff(default: dict, current: dict) -> dict:
    """Compares two dicts and returns only modified or new values."""
    diff = {}
    for key, val in current.items():
        if key not in default:
            diff[key] = val
        elif isinstance(val, dict):
            nested_diff = get_config_diff(default.get(key, {}), val)
            if nested_diff:
                diff[key] = nested_diff
        elif val != default[key]:
            diff[key] = val
    return diff


class AtomViewerGUI(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()

        self._is_syncing = False

        self._config_cache = {}
        self.load_initial_config_from_disk()

        self.setWindowTitle("C.E.M.D | Computational Elementary Matter Design")
        self.resize(1400, 900)

        self.setup_ui()
        self.setup_menu()

    @property
    def system(self) -> AtomicSystem | None:
        """Dynamically retrieves system from the active tab."""
        if not hasattr(self, "tabs"):
            return None
        tab = self.tabs.currentWidget()
        return tab.system if tab else None

    @system.setter
    def system(self, value: AtomicSystem) -> None:
        tab = self.tabs.currentWidget()
        if tab:
            tab.system = value

    @property
    def current_file_path(self) -> str | None:
        """Dynamically retrieve the file path of the active tab."""
        tab = self.tabs.currentWidget()
        return getattr(tab, "file_path", None) if tab else None

    @current_file_path.setter
    def current_file_path(self, value: str) -> None:
        """Updates the path of the active tab."""
        tab = self.tabs.currentWidget()
        if tab:
            tab.file_path = value

    @property
    def global_config(self) -> dict[str, Any]:
        """Reads the current merged state (Cache + Plotters + UI)."""

        return self._build_config()

    @global_config.setter
    def global_config(self, value: dict[str, Any]) -> None:
        """Injects a config dictionary DIRECTLY into the UI widgets."""
        if not value or not isinstance(value, dict) or self._is_syncing:
            return

        self._config_cache.update(value)

        # PHYSICAL update of the FilterPanel (Sliders/Labels)
        fp = self.filter_panel
        fp.blockSignals(True)  # Avoid looping refreshes during the update

        scale = value.get("global_scale", 1.0)
        fp.global_scale_slider.setValue(int(scale * 100))
        fp.lbl_scale_val.setText(f"x{scale:.1f}")

        fp.blockSignals(False)

        # PHYSICAL update of the BondManager
        bm = self.bond_manager
        bm.blockSignals(True)

        b_radius = value.get("global_bond_radius", 0.1)
        bm.bond_radius_slider.setValue(int(b_radius * 100))
        # bm.refresh_from_map(value.get("bond_map", {}))

        bm.blockSignals(False)

        # In the global_config setter (MainWindow)
        new_bg = value.get("bg_color")
        if new_bg:
            for i in range(self.tabs.count()):
                tab = self.tabs.widget(i)
                p = tab.plotter
                # apply the color physically
                p.set_background(new_bg)
                if new_bg in p.bg_cycle:
                    p.bg_idx = p.bg_cycle.index(new_bg)
                p.update_background_style()

    def _build_config(self) -> dict[str, Any]:
        """Construit et retourne la config mergée sans modifier l'état."""
        cfg = self._config_cache.copy()

        for i in range(self.tabs.count()):
            tab = self.tabs.widget(i)
            p = tab.plotter
            p_cfg = p.get_current_config_dict()
            cfg.setdefault("color_map", {}).update(p_cfg["color_map"])
            cfg.setdefault("radius_map", {}).update(p_cfg["radius_map"])
            cfg.setdefault("bond_map", {}).update(p_cfg["bond_map"])
            cfg["global_scale"] = p_cfg["global_scale"]
            cfg["global_bond_radius"] = p_cfg["global_bond_radius"]
            cfg["bg_color"] = p_cfg["bg_color"]

        f_settings = self.filter_panel.get_settings()
        cfg["global_scale"] = f_settings.get("global_scale", 100) / 100.0
        cfg.setdefault("radius_map", {}).update(f_settings.get("radii", {}))

        b_settings = self.bond_manager.get_settings()
        cfg.setdefault("bond_map", {}).update(
            {k: v for k, v in b_settings.items() if k != "global_bond_radius"}
        )
        cfg["global_bond_radius"] = self.bond_manager.get_bond_radius()

        return cfg

    def load_initial_config_from_disk(self) -> None:
        """Loads the default THEN the user modifications."""
        base_dir = os.path.dirname(os.path.realpath(__file__))

        default_path = os.path.join(base_dir, "default_config.json")
        # Preferences live under the user's own config directory. Fall back
        # to the copy earlier versions left inside the package, so that an
        # upgrade does not discard settings already saved.
        user_path = _userdata.config_file()
        if not user_path.exists():
            legacy = _userdata.legacy_config_file()
            if legacy.exists():
                user_path = legacy

        merged_cfg = {}

        # Load the default base (Unchangeable)
        if os.path.exists(default_path):
            try:
                with open(default_path, encoding="utf-8") as f:
                    merged_cfg = json.load(f) or {}
            except Exception as e:
                print(f"Erreur lecture default_config: {e}")

        # Overwrite with user preferences (Priority)
        if user_path.exists():
            try:
                with open(user_path, encoding="utf-8") as f:
                    user_cfg = json.load(f) or {}
                    # use update to merge the first level dictionaries
                    # Note: if you have nested dicts (ex: color_map),
                    # it may need a recursive merge.
                    deep_update(merged_cfg, user_cfg)
            except Exception as e:
                print(f"Erreur lecture config.json: {e}")

        self._config_cache = merged_cfg

    def setup_ui(self) -> None:
        """Initialize the main UI layout, tab manager, and side control panels."""
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        main_layout = QtWidgets.QHBoxLayout(central)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.setTabsClosable(True)
        self.tabs.setMovable(True)

        # Tab Signal Connections
        self.tabs.currentChanged.connect(self.on_tab_changed)
        self.tabs.tabCloseRequested.connect(self.request_close_tab)

        main_layout.addWidget(self.tabs, stretch=3)

        right_panel = QtWidgets.QWidget()
        right_layout = QtWidgets.QVBoxLayout(right_panel)
        right_layout.setContentsMargins(5, 5, 5, 5)
        right_layout.setSpacing(10)

        self.system_summary = SystemSummaryPanel()
        right_layout.addWidget(self.system_summary)

        self.filter_panel = FilterPanel()
        # synchronize the UI when check/uncheck a type
        self.filter_panel.type_changed.connect(
            lambda settings: self.sync_ui(full_rebuild=False)
        )
        right_layout.addWidget(self.filter_panel, stretch=2)

        # Bonds manager
        self.bond_manager = BondManagerPanel()
        self.bond_manager.bond_settings_changed.connect(self.refresh_bonds_view)
        right_layout.addWidget(self.bond_manager, stretch=2)

        right_layout.addStretch()

        # Setting up scrolling for the right panel
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(right_panel)
        scroll.setMinimumWidth(350)
        scroll.setFrameShape(QtWidgets.QFrame.NoFrame)

        main_layout.addWidget(scroll)

        # setup tab renaming
        self.tabs.tabBar().tabBarDoubleClicked.connect(self._edit_tab_text)

    def setup_menu(self) -> None:
        """Create toolbar actions, assign icons, and connect signals to handlers."""
        self.tools_toolbar = self.addToolBar("Manipulation Tools")
        self.tools_toolbar.setIconSize(QtCore.QSize(24, 24))
        self._tab_sensitive_actions = []

        def add_action(
            attr_name, name, icon, slot, sensitive=True, tip=None, shortcut=None
        ):
            """Create the action, bind it to self.attr_name and add it to the toolbar.

            `tip` is what the toolbar shows on hover and the status bar
            shows on focus: the icons carry no label, so without it the
            only clue to what a button does is its one-word name.
            """
            action = QtGui.QAction(get_icon(icon), name, self)
            action.triggered.connect(slot)

            if shortcut is not None:
                action.setShortcut(QtGui.QKeySequence(shortcut))
                # Qt appends the shortcut to the tooltip only when one is
                # set explicitly, so build it here to get "Open (Ctrl+O)".
                action.setShortcutContext(QtCore.Qt.ApplicationShortcut)

            if tip is not None:
                hint = action.shortcut().toString()
                action.setToolTip(f"{tip} ({hint})" if hint else tip)
                action.setStatusTip(tip)

            self.tools_toolbar.addAction(action)

            setattr(self, attr_name, action)

            if sensitive:
                self._tab_sensitive_actions.append(action)
            return action

        add_action(
            "action_open",
            "Open",
            "open",
            self.open_file_clicked,
            sensitive=False,
            tip="Open a structure file (LAMMPS data, PDB, CIF, moltemplate, SDF)",
            shortcut="Ctrl+O",
        )
        add_action(
            "action_save",
            "Save",
            "save",
            self.save_file_clicked,
            tip="Save the active structure to its current file",
            shortcut="Ctrl+S",
        )
        add_action(
            "action_save_as",
            "Save as",
            "save-all",
            self.save_file_as_clicked,
            tip="Save the active structure under a new name",
            shortcut="Ctrl+Shift+S",
        )
        add_action(
            "action_bg_color",
            "Background color",
            "bg-color",
            self.on_cycle_bg_clicked,
            tip="Cycle the 3D view background color",
        )
        add_action(
            "action_reset_camera",
            "Reset camera",
            "camera",
            lambda: self.sync_ui(reset_camera=True),
            tip="Reset the camera to fit the whole structure",
            shortcut="Ctrl+R",
        )

        self.tools_toolbar.addSeparator()

        add_action(
            "action_cod",
            "COD",
            "cod",
            self.open_cod_browser,
            sensitive=False,
            tip="Search the Crystallography Open Database and import a structure",
        )
        add_action(
            "action_pubchem",
            "PubChem",
            "pubchem",
            self.open_pubchem_browser,
            sensitive=False,
            tip="Search PubChem and import a molecule by name or formula",
        )
        add_action(
            "action_smiles",
            "SMILES",
            "mol",
            lambda: open_smiles_builder(self),
            sensitive=False,
            tip="Build a molecule from a SMILES string",
        )

        self.tools_toolbar.addSeparator()

        add_action(
            "action_solution",
            "Solution",
            "solution",
            lambda: open_make_solution(self),
            sensitive=False,
            tip="Pack a solution of given density and composition",
        )
        add_action(
            "action_glass",
            "Glass",
            "glass",
            lambda: open_make_glass(self),
            sensitive=False,
            tip="Build an amorphous oxide glass",
        )
        add_action(
            "action_cash",
            "Cement hydrates",
            "cash",
            lambda: open_make_cash(self),
            sensitive=False,
            tip="Build a C-S-H or C-A-S-H structure from tobermorite",
        )

        self.tools_toolbar.addSeparator()

        add_action(
            "action_replicate",
            "Replicate",
            "replicate",
            lambda: open_replicate(self),
            tip="Replicate the cell along a, b and c",
        )
        add_action(
            "action_orthogonalize",
            "Orthogonalize",
            "orthogonalize",
            self.on_orthogonalize_clicked,
            tip="Convert the cell to an orthogonal box",
        )
        add_action(
            "action_translate",
            "Translate",
            "move",
            lambda: open_translate_atoms(self),
            tip="Translate selected atoms by a vector",
        )
        add_action(
            "action_center",
            "Center/Wrap",
            "wrap",
            self.on_center_clicked,
            tip="Centre the structure and wrap atoms into the box",
        )
        add_action(
            "action_surface",
            "Surface",
            "surface",
            lambda: open_make_surface(self),
            tip="Cut a surface slab from Miller indices",
        )
        add_action(
            "action_add_structure",
            "Add structure",
            "add_mol",
            lambda: open_add_structure(self),
            tip="Insert another structure into the active system",
        )
        add_action(
            "action_add_liquid",
            "Add liquid",
            "interface",
            lambda: open_add_liquid(self),
            tip="Add a liquid layer on top of the structure",
        )
        add_action(
            "action_add_droplet",
            "Add droplet",
            "droplet",
            lambda: open_add_droplet(self),
            tip="Add a hemispherical droplet on the structure",
        )
        add_action(
            "action_split",
            "Split",
            "channel",
            lambda: open_split(self),
            tip="Open a gap along an axis and optionally fill it with a solution",
        )
        add_action(
            "action_protonate",
            "Protonate",
            "protonate",
            lambda: on_protonate(self),
            tip="Add hydrogens to under-coordinated oxygens",
        )

        self.tools_toolbar.addSeparator()

        add_action(
            "action_undo",
            "Undo",
            "reference",
            self.undo_clicked,
            tip="Undo the last structural change",
            shortcut="Ctrl+Z",
        )
        add_action(
            "action_type_manager",
            "Types",
            "atom",
            self.on_type_manager_clicked,
            tip="Edit atom types, masses and charges",
        )
        add_action(
            "action_connectivity_manager",
            "Connectivity",
            "connectivity",
            self.on_connectivity_manager_clicked,
            tip="Edit bonds, angles and topology rules",
        )

        self.tools_toolbar.addSeparator()

        add_action(
            "action_rdf_analysis",
            "RDF",
            "rdf",
            self.open_rdf_analysis,
            tip="Compute a radial distribution function",
        )
        add_action(
            "action_silicate_analysis",
            "Silicate",
            "silicate",
            self.open_silicate_analysis,
            tip="Analyse the silicate network (Ca/Si, Qn, MCL)",
        )

        self.setup_menubar()
        self.set_tools_enabled(False)

    def setup_menubar(self) -> None:
        """Group the toolbar actions into a menu bar.

        The very same QAction objects are reused, so a menu entry and its
        toolbar button stay in lockstep -- same icon, same tooltip, same
        enabled state, and a single shortcut that works from either. The
        toolbar shows 25 unlabelled icons and hides the overflow behind a
        chevron on a narrow window; the menus are where the features can
        be found by name.
        """
        menus = {
            "&File": [
                self.action_open,
                self.action_save,
                self.action_save_as,
                None,
                self.action_cod,
                self.action_pubchem,
                self.action_smiles,
            ],
            "&Edit": [
                self.action_undo,
                None,
                self.action_type_manager,
                self.action_connectivity_manager,
            ],
            "&Build": [
                self.action_solution,
                self.action_glass,
                self.action_cash,
                None,
                self.action_surface,
                self.action_add_structure,
                self.action_add_liquid,
                self.action_add_droplet,
                self.action_split,
                self.action_protonate,
            ],
            "&Transform": [
                self.action_replicate,
                self.action_orthogonalize,
                self.action_translate,
                self.action_center,
            ],
            "&Analyze": [
                self.action_rdf_analysis,
                self.action_silicate_analysis,
            ],
            "&View": [
                self.action_reset_camera,
                self.action_bg_color,
            ],
        }

        bar = self.menuBar()
        for title, entries in menus.items():
            menu = bar.addMenu(title)
            for entry in entries:
                if entry is None:
                    menu.addSeparator()
                else:
                    menu.addAction(entry)

    # ------------------------------------------------------------------
    # Undo
    # ------------------------------------------------------------------

    #: How many structures back the undo history reaches. Snapshots are
    #: whole systems, so the bound is a memory bound, not a preference.
    UNDO_DEPTH = 5

    def push_undo(self) -> None:
        """Snapshot the active structure before a destructive edit.

        Call this *before* mutating a system, from every handler that
        changes one in place: replicating, protonating, orthogonalizing,
        adding a liquid, splitting. Operations that open a new tab leave
        the original untouched and need no snapshot.
        """
        tab = self.tabs.currentWidget()
        if tab is None or getattr(tab, "system", None) is None:
            return

        stack = getattr(tab, "undo_stack", None)
        if stack is None:
            stack = tab.undo_stack = []

        stack.append(tab.system.copy())
        # Keep only the last UNDO_DEPTH snapshots.
        del stack[: -self.UNDO_DEPTH]
        self.refresh_undo_state()

    @QtCore.Slot()
    def undo_clicked(self) -> None:
        """Restore the structure as it was before the last edit."""
        tab = self.tabs.currentWidget()
        stack = getattr(tab, "undo_stack", None) if tab else None

        if not stack:
            self.statusBar().showMessage("Nothing to undo.", 3000)
            return

        tab.system = stack.pop()
        self.sync_ui(full_rebuild=True)
        self.statusBar().showMessage("Undone.", 3000)
        self.refresh_undo_state()

    def refresh_undo_state(self) -> None:
        """Grey out Undo when the active tab has no history."""
        tab = self.tabs.currentWidget()
        stack = getattr(tab, "undo_stack", None) if tab else None
        self.action_undo.setEnabled(bool(stack))

    def set_tools_enabled(self, state: bool = False) -> None:
        """Enable or disable all manipulation tools based on the system state."""
        for action in self._tab_sensitive_actions:
            action.setEnabled(state)

        # Undo follows its own history rather than the "a tab exists" rule.
        self.refresh_undo_state()

        self.update_protonate_state()

    @QtCore.Slot()
    def open_file_clicked(self) -> None:
        """The user clicked Open."""
        system, path = open_file(self)
        if system:
            self.add_structure_tab(system, path)
            self.statusBar().showMessage(f"Chargé : {path}", 3000)

    @QtCore.Slot()
    def save_file_clicked(self) -> None:
        """The user clicked Save."""
        if not self.system:
            return

        if self.current_file_path:
            save_file(self, self.system, self.current_file_path)
        else:
            self.save_file_as_clicked()

    @QtCore.Slot()
    def save_file_as_clicked(self) -> None:
        """The user clicked Save As."""
        if not self.system:
            return

        new_path = save_file_as(self, self.system, self.current_file_path)
        if new_path:
            self.current_file_path = new_path
            self.tabs.setTabText(self.tabs.currentIndex(), os.path.basename(new_path))

    def on_tab_changed(self, index: int) -> None:
        """Synchronize global configuration and UI settings when switching between structure tabs."""
        if index == -1 or self._is_syncing:
            return
        current_tab = self.tabs.widget(index)

        if current_tab and hasattr(current_tab, "plotter"):
            # Complete synchro
            # recover EVERYTHING stored in the MainWindow
            cfg = self.global_config

            # inject the global settings into the plotter of the incoming tab
            current_tab.plotter.color_map.update(cfg.get("color_map", {}))
            current_tab.plotter.radius_map.update(cfg.get("radius_map", {}))
            current_tab.plotter.bond_map.update(cfg.get("bond_map", {}))

            # also synchronize the scalar values
            current_tab.plotter.global_scale = cfg.get(
                "global_scale", current_tab.plotter.global_scale
            )
            current_tab.plotter.global_bond_radius = cfg.get(
                "global_bond_radius", current_tab.plotter.global_bond_radius
            )

            # UI UPDATE (Sliders <---Plotter)
            # The setter will now place the cursors in the correct positions
            self.global_config = current_tab.plotter.get_current_config_dict()

            # DEFINITION OF ACTIVE SYSTEM AND DRAWING
            self.system = current_tab.system
            self.sync_ui(full_rebuild=True, reset_camera=True)

    def add_structure_tab(
        self, system: AtomicSystem, path: str = None, title: str = None
    ) -> None:
        """Create and append a new StructureTabWidget to the tab manager."""

        if title:
            name = title
        elif path:
            filename = os.path.basename(path)
            name = os.path.splitext(filename)[0]
        else:
            name = f"# {self.tabs.count() + 1}"

        new_tab = StructureTabWidget(system, self, file_path=path)

        self._is_syncing = True

        try:
            self.tabs.blockSignals(True)
            index = self.tabs.addTab(new_tab, name)

            self.tabs.setCurrentIndex(index)
            self.set_tools_enabled(True)

            self.tabs.blockSignals(False)
        except Exception as e:
            print(f"CRASH dans sync_ui: {e}")
        finally:
            self._is_syncing = False

        self.sync_ui(full_rebuild=True, reset_camera=True)

    def request_close_tab(self, index: int) -> None:
        """Request confirmation before closing a tab."""
        name = self.tabs.tabText(index)

        msg = QtWidgets.QMessageBox.question(
            self,
            "Close structure",
            f"Do you want to close '{name}'?\nMake sure you have saved your changes.",
            QtWidgets.QMessageBox.StandardButton.Yes
            | QtWidgets.QMessageBox.StandardButton.No,
        )

        if msg == QtWidgets.QMessageBox.StandardButton.Yes:
            self.tabs.removeTab(index)
            # If we closed everything, we gray out the tools
            if self.tabs.count() == 0:
                self.set_tools_enabled(False)

    def sync_ui(
        self, full_rebuild=True, reset_camera=False, refresh_bonds=True
    ) -> None:
        """Update side panels and trigger a refresh of the active tab's 3D visualization."""
        active_tab = self.tabs.currentWidget()

        # If there is no tab (at the first start), we shut down immediately
        if active_tab is None or self._is_syncing:
            return

        # Now we can work with active_tab.system
        if self.system is None:
            return

        try:
            self._is_syncing = True

            # SYNCHRO OF SHARED PANELS (MAIN WINDOW)
            if full_rebuild:
                self.filter_panel.refresh(
                    self.system,
                    active_tab.plotter.color_map,
                    active_tab.plotter.radius_map,
                    active_tab.plotter.global_scale,
                )
                self.system_summary.update_info(self.system)

                # Updating the Global Topology Manager
                if refresh_bonds:
                    self.bond_manager.refresh(self.system, active_tab.plotter.bond_map)

            # DELEGATION TO ACTIVE TAB
            active_tab.refresh_tab_view(
                full_rebuild=full_rebuild,
                reset_camera=reset_camera,
                refresh_bonds=refresh_bonds,
            )

        finally:
            self._is_syncing = False

    def refresh_bonds_view(self) -> None:
        """Ask the active tab to update the joint drawing."""
        active_tab = self.tabs.currentWidget()
        # retrieve the settings from the right panel
        settings = self.bond_manager.get_settings()

        # ask the tab plotter to redraw only the jumps
        # using its own df_visible
        active_tab.plotter.update_bonds_only(active_tab.df_visible, settings)
        active_tab.plotter.render()

    def update_protonate_state(self) -> None:
        """Activates the protonate icon if at least one atom is selected in the active tab."""
        active_tab = self.tabs.currentWidget()

        # By default, it is disabled (if there is no tab)
        has_selection = False

        # If we have a tab open, we look at its local selection list
        if active_tab and hasattr(active_tab, "selected_real_indices"):
            has_selection = len(active_tab.selected_real_indices) > 0

        self.action_protonate.setEnabled(has_selection)

        # force the visual update of the toolbar
        self.tools_toolbar.update()

    @QtCore.Slot()
    def on_orthogonalize_clicked(self) -> None:
        self.push_undo()
        try:
            self.system.orthogonalize()
        except ValueError:
            self.system.unskew()
            self.statusBar().showMessage(
                "Could not orthogonalize — unskew applied instead.", 5000
            )
        self.sync_ui(refresh_bonds=True)

    @QtCore.Slot()
    def on_center_clicked(self) -> None:
        self.push_undo()
        self.system.center_on_com()
        self.sync_ui(refresh_bonds=True)

    @QtCore.Slot()
    def on_type_manager_clicked(self) -> None:
        TypeManagerDialog(self).exec_()

    @QtCore.Slot()
    def on_connectivity_manager_clicked(self) -> None:
        ConnectivityDialog(self).exec_()

    @QtCore.Slot()
    def on_cycle_bg_clicked(self) -> None:
        active_tab = self.tabs.currentWidget()
        if not active_tab:
            return

        # change the index
        new_idx = (active_tab.plotter.bg_idx + 1) % len(active_tab.plotter.bg_cycle)

        # synchronize all the tabs
        for i in range(self.tabs.count()):
            tab = self.tabs.widget(i)
            tab.plotter.bg_idx = new_idx
            tab.plotter.update_background_style()

    @QtCore.Slot()
    def open_rdf_analysis(self) -> None:
        """Opens the RDF calculation dialog on the active tab."""
        # Retrieve system from the active tab via property
        current_system = self.system

        if current_system is None:
            QtWidgets.QMessageBox.warning(
                self,
                "Missing system",
                "Please load a structure before running RDF parsing.",
            )
            return

        try:
            rdf_diag = RDFDialog(self, current_system)
            rdf_diag.exec()

        except Exception as e:
            QtWidgets.QMessageBox.critical(
                self,
                "Conversion error",
                f"Unable to convert the atomic system for MDAnalysis:\n{e}",
            )

    @QtCore.Slot()
    def open_silicate_analysis(self):
        current_system = self.system

        if current_system:
            dialog = SilicateDialog(self, current_system)
            dialog.exec()

    @QtCore.Slot()
    def open_cod_browser(self):
        dialog = CODBrowserDialog(self)
        if dialog.exec() == QtWidgets.QDialog.Accepted:
            system = dialog.selected_system  # ← already an AtomicSystem
            cod_id = dialog.last_cod_id
            title = f"COD_{cod_id}"
            self.add_structure_tab(system, title=title)
            self.statusBar().showMessage(f"Structure COD {cod_id} imported.", 5000)

    @QtCore.Slot()
    def open_pubchem_browser(self) -> None:
        """Launches the PubChem browser and creates a tab if a structure is chosen."""
        dialog = PubChemBrowserDialog(self)

        if dialog.exec() == QtWidgets.QDialog.Accepted:
            # retrieve the topology_dict generated by the SDF reader
            system = dialog.selected_system

            selected_row = dialog.table.currentRow()
            cid = dialog.table.item(selected_row, 0).text()
            name = dialog.table.item(selected_row, 1).text()
            formula = dialog.table.item(selected_row, 2).text()

            self.add_structure_tab(system, title=f"PubChem_{cid}")

            self.statusBar().showMessage(
                f"Compound {name} (CID {cid}) imported from PubChem.", 5000
            )

    def save_global_to_json(self) -> None:
        """Only saves differences from the default."""
        base_dir = os.path.dirname(os.path.realpath(__file__))
        default_path = os.path.join(base_dir, "default_config.json")
        user_path = _userdata.config_file()

        # Load fault to compare
        default_data = {}
        if os.path.exists(default_path):
            with open(default_path) as f:
                default_data = json.load(f)

        # Calculate only what is different from the original
        current_data = self.global_config
        diff_data = get_config_diff(default_data, current_data)

        # Write delta to config.json
        try:
            with open(user_path, "w") as f:
                json.dump(diff_data, f, indent=4)
            print(f"Saved user configuration ({user_path})")
        except Exception as e:
            print(f"JSON writing error: {e}")

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        """Centralized global backup on shutdown."""
        try:
            self.save_global_to_json()

            for i in range(self.tabs.count()):
                tab = self.tabs.widget(i)
                if tab and hasattr(tab, "plotter"):
                    tab.plotter.close()

            for name in ["cod_cache.json", "pubchem_cache.json"]:
                _userdata.cache_file(name).unlink(missing_ok=True)

        except Exception as e:
            print(f"Error closing globally: {e}")
            import traceback

            traceback.print_exc()

        event.accept()

    def _edit_tab_text(self, index):
        if index < 0:
            return

        # get the current name
        bar = self.tabs.tabBar()
        current_text = self.tabs.tabText(index)

        # Creating the edit field
        line_edit = QtWidgets.QLineEdit(current_text)
        line_edit.setFrame(False)
        line_edit.setMaximumWidth(150)
        line_edit.setStyleSheet(
            "background: palette(window); border: 1px solid palette(highlight);"
        )

        # inject the widget into the tab
        bar.setTabButton(index, QtWidgets.QTabBar.LeftSide, line_edit)
        self.tabs.setTabText(index, "")  # Hide static text

        line_edit.setFocus()
        line_edit.selectAll()

        # Validation logic
        self._renaming_active = True

        def save_name():
            if not hasattr(self, "_renaming_active") or not self._renaming_active:
                return
            self._renaming_active = False

            new_name = line_edit.text().strip() or current_text

            # remove the widget and put the text back
            bar.setTabButton(index, QtWidgets.QTabBar.LeftSide, None)
            self.tabs.setTabText(index, new_name)

            # ---System update ---
            tab_widget = self.tabs.widget(index)
            tab_widget.system.name = new_name

        line_edit.returnPressed.connect(save_name)
        line_edit.editingFinished.connect(save_name)


def main() -> int:
    """Launch the graphical interface.

    Exposed as the ``cemd-gui`` command so the interface can be started
    without knowing where the package was installed.
    """
    app = QtWidgets.QApplication(sys.argv)

    style = QlementineStyle(app)
    app.setStyle(style)

    # Force dot as decimal separator
    qt_locale = QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates)
    QtCore.QLocale.setDefault(qt_locale)

    window = AtomViewerGUI()
    window.show()

    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
