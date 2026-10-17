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

from cemd import __version__
from cemd.core.atomic_system import AtomicSystem
from cemd.gui import _userdata, display_settings
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


class AtomViewerGUI(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()

        self._is_syncing = False

        #: How every tab draws its structure; see `display_settings`.
        self.display = display_settings.load()

        self.setWindowTitle(
            f"cemd {__version__} \u2014 atomistic model builder for LAMMPS"
        )
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
        self.filter_panel.set_scale_value(self.display.global_scale)
        self.filter_panel.type_changed.connect(self.on_view_settings_changed)
        right_layout.addWidget(self.filter_panel, stretch=2)

        # Bonds manager
        self.bond_manager = BondManagerPanel()
        self.bond_manager.set_bond_radius(self.display.global_bond_radius)
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
        adding a liquid, splitting, and the type and connectivity
        managers. Operations that open a new tab leave the original
        untouched and need no snapshot.

        The tab's force-field choices (`session_ff_dict`) are kept with
        the system: they live on the tab, not in the system, and the
        type manager rewrites both.
        """
        tab = self.tabs.currentWidget()
        if tab is None or getattr(tab, "system", None) is None:
            return

        stack = getattr(tab, "undo_stack", None)
        if stack is None:
            stack = tab.undo_stack = []

        stack.append((tab.system.copy(), dict(getattr(tab, "session_ff_dict", {}))))
        # Keep only the last UNDO_DEPTH snapshots.
        del stack[: -self.UNDO_DEPTH]
        self.refresh_undo_state()
        # Every destructive edit comes through here first, which makes it
        # the one place that has to know the structure is about to change.
        self.set_tab_modified(tab, True)

    def discard_undo(self) -> None:
        """Drop the snapshot `push_undo` just took, for an edit that failed."""
        tab = self.tabs.currentWidget()
        stack = getattr(tab, "undo_stack", None) if tab else None
        if stack:
            stack.pop()
        self.refresh_undo_state()

    @QtCore.Slot()
    def undo_clicked(self) -> None:
        """Restore the structure as it was before the last edit."""
        tab = self.tabs.currentWidget()
        stack = getattr(tab, "undo_stack", None) if tab else None

        if not stack:
            self.statusBar().showMessage("Nothing to undo.", 3000)
            return

        tab.system, session_ff = stack.pop()
        tab.session_ff_dict = session_ff
        self.sync_ui(full_rebuild=True)
        self.statusBar().showMessage("Undone.", 3000)
        self.refresh_undo_state()
        # Undo may or may not land back on the saved state; nothing records
        # which snapshot was saved, so stay on the safe side.
        self.set_tab_modified(tab, True)

    #: Shown on the tab of a structure with unsaved changes. The title
    #: itself stays clean: other dialogs list structures by tab title.
    _MODIFIED_TIP = "Unsaved changes"

    def set_tab_modified(self, tab: QtWidgets.QWidget, modified: bool) -> None:
        """Record whether `tab` holds changes that are not on disk."""
        if tab is None:
            return
        tab.modified = modified
        index = self.tabs.indexOf(tab)
        if index < 0:
            return
        self.tabs.setTabIcon(
            index, self._modified_icon() if modified else QtGui.QIcon()
        )
        self.tabs.setTabToolTip(index, self._MODIFIED_TIP if modified else "")

    def _modified_icon(self) -> QtGui.QIcon:
        """A small dot in the theme's accent color."""
        if getattr(self, "_modified_icon_cache", None) is None:
            size = 10
            pixmap = QtGui.QPixmap(size, size)
            pixmap.fill(QtCore.Qt.transparent)
            painter = QtGui.QPainter(pixmap)
            painter.setRenderHint(QtGui.QPainter.Antialiasing)
            painter.setPen(QtCore.Qt.NoPen)
            painter.setBrush(self.palette().color(QtGui.QPalette.Highlight))
            painter.drawEllipse(1, 1, size - 2, size - 2)
            painter.end()
            self._modified_icon_cache = QtGui.QIcon(pixmap)
        return self._modified_icon_cache

    def modified_tabs(self) -> list[QtWidgets.QWidget]:
        """The open tabs whose structure has unsaved changes."""
        tabs = (self.tabs.widget(i) for i in range(self.tabs.count()))
        return [tab for tab in tabs if getattr(tab, "modified", False)]

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
            self.statusBar().showMessage(f"Loaded: {path}", 3000)

    @QtCore.Slot()
    def save_file_clicked(self) -> None:
        """The user clicked Save."""
        if not self.system:
            return

        if self.current_file_path:
            if save_file(self, self.system, self.current_file_path):
                self.set_tab_modified(self.tabs.currentWidget(), False)
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
            self.set_tab_modified(self.tabs.currentWidget(), False)

    def on_tab_changed(self, index: int) -> None:
        """Show the incoming tab's structure in the side panels and redraw it.

        The display settings are shared, so there is nothing to carry
        over from one tab to the next: the incoming plotter only has to
        catch up with changes made while it was hidden.
        """
        if index == -1 or self._is_syncing:
            return
        current_tab = self.tabs.widget(index)

        if current_tab and hasattr(current_tab, "plotter"):
            current_tab.plotter.update_background_style()
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
        # A structure read from a file matches it; one built or downloaded
        # here exists nowhere else until it is saved.
        new_tab.modified = path is None

        self._is_syncing = True

        try:
            self.tabs.blockSignals(True)
            index = self.tabs.addTab(new_tab, name)
            self.set_tab_modified(new_tab, new_tab.modified)

            self.tabs.setCurrentIndex(index)
            self.set_tools_enabled(True)
        except Exception as e:
            print(f"Could not add the tab: {e}")
        finally:
            # In `finally`: a failure above used to leave the tab bar mute,
            # and switching tabs then stopped refreshing anything.
            self.tabs.blockSignals(False)
            self._is_syncing = False

        self.sync_ui(full_rebuild=True, reset_camera=True)

    def request_close_tab(self, index: int) -> None:
        """Close a tab, offering to save it first if it has unsaved changes.

        A saved structure closes at once: asking every time trained users
        to click Yes without reading, which is what lost changes.
        """
        tab = self.tabs.widget(index)
        if getattr(tab, "modified", False) and not self._confirm_discard(tab):
            return

        # By widget, not by the index passed in: saving may have run in
        # between.
        self.tabs.removeTab(self.tabs.indexOf(tab))
        # removeTab only hides the page. Without this, each closed tab
        # kept its render window, its GPU buffers and its system (and
        # undo snapshots) alive until the application quit.
        if tab is not None:
            if hasattr(tab, "plotter"):
                tab.plotter.close()
            tab.deleteLater()
        # If we closed everything, we gray out the tools
        if self.tabs.count() == 0:
            self.set_tools_enabled(False)

    def _confirm_discard(self, tab: QtWidgets.QWidget) -> bool:
        """Ask what to do with `tab`'s unsaved changes; True to go ahead."""
        buttons = QtWidgets.QMessageBox.StandardButton
        name = self.tabs.tabText(self.tabs.indexOf(tab))
        answer = QtWidgets.QMessageBox.warning(
            self,
            "Unsaved changes",
            f"'{name}' has changes that are not saved.\nSave them before closing?",
            buttons.Save | buttons.Discard | buttons.Cancel,
            buttons.Save,
        )
        if answer == buttons.Discard:
            return True
        if answer == buttons.Save:
            self.tabs.setCurrentWidget(tab)
            self.save_file_clicked()
            # Cancelling the file dialog leaves the tab modified: keep it.
            return not getattr(tab, "modified", False)
        return False

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
                    self.display.color_map,
                    self.display.radius_map,
                    self.display.global_scale,
                )
                # A type the settings do not know yet gets the panel's
                # default color and radius; record them so the drawing
                # uses the values the panel shows.
                self._store_view_settings(self.filter_panel.get_settings())
                self.system_summary.update_info(self.system)

                # Updating the Global Topology Manager
                if refresh_bonds:
                    self.bond_manager.refresh(self.system, self.display.bond_map)

            # DELEGATION TO ACTIVE TAB
            active_tab.refresh_tab_view(
                full_rebuild=full_rebuild,
                reset_camera=reset_camera,
                refresh_bonds=refresh_bonds,
            )

        finally:
            self._is_syncing = False

    def _store_view_settings(self, panel_settings: dict[str, Any]) -> None:
        """Copy the type panel's colors, radii and scale into the settings."""
        self.display.color_map.update(panel_settings.get("colors", {}))
        self.display.radius_map.update(panel_settings.get("radii", {}))
        self.display.global_scale = panel_settings.get("global_scale", 100) / 100.0

    def on_view_settings_changed(self, panel_settings: dict[str, Any]) -> None:
        """The user changed a color, a radius, the scale or a visible type."""
        self._store_view_settings(panel_settings)
        self.sync_ui(full_rebuild=False)

    def refresh_bonds_view(self, panel_settings: dict[str, Any] | None = None) -> None:
        """Store the bond panel's rules and redraw the active tab's bonds."""
        active_tab = self.tabs.currentWidget()
        if active_tab is None or self.system is None:
            return

        panel_settings = dict(panel_settings or self.bond_manager.get_settings())
        self.display.global_bond_radius = panel_settings.pop(
            "global_bond_radius", self.display.global_bond_radius
        )
        self.display.set_bond_rules(self.system.atom_types, panel_settings)

        active_tab.refresh_tab_view(full_rebuild=False)

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

        self.display.next_background()
        # The other tabs catch up when they are shown (`on_tab_changed`).
        active_tab.plotter.update_background_style()

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
        """Save the display settings that differ from the shipped defaults."""
        try:
            path = display_settings.save(self.display)
            print(f"Saved user configuration ({path})")
        except OSError as e:
            print(f"Could not save the user configuration: {e}")

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        """Centralized global backup on shutdown."""
        unsaved = self.modified_tabs()
        if unsaved:
            buttons = QtWidgets.QMessageBox.StandardButton
            names = "\n".join(
                f"  \u2022 {self.tabs.tabText(self.tabs.indexOf(t))}" for t in unsaved
            )
            answer = QtWidgets.QMessageBox.warning(
                self,
                "Unsaved changes",
                f"{len(unsaved)} structure(s) have changes that are not saved:\n"
                f"{names}\n\nQuit anyway and lose them?",
                buttons.Discard | buttons.Cancel,
                buttons.Cancel,
            )
            if answer != buttons.Discard:
                event.ignore()
                return

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


class _TooltipColorFixer(QtCore.QObject):
    """Give each tooltip the theme's colors as it appears.

    Qlementine paints the tooltip background with its dark theme color,
    but Qt keeps resetting the tooltip palette to the system default --
    black text -- after startup, so ``QToolTip.setPalette`` does not
    stick. Setting the palette on the label itself when it is shown is
    the only point that wins.
    """

    def __init__(self, style: QlementineStyle) -> None:
        super().__init__(style)
        self._style = style

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:
        if (
            event.type() == QtCore.QEvent.Type.Show
            and obj.metaObject().className() == "QTipLabel"
        ):
            theme = self._style.theme()
            palette = obj.palette()
            palette.setColor(QtGui.QPalette.ColorRole.ToolTipBase, theme.secondaryColor)
            palette.setColor(
                QtGui.QPalette.ColorRole.ToolTipText, theme.secondaryColorForeground
            )
            obj.setPalette(palette)
        return False


def main() -> int:
    """Launch the graphical interface.

    Exposed as the ``cemd-gui`` command so the interface can be started
    without knowing where the package was installed.
    """
    app = QtWidgets.QApplication(sys.argv)

    style = QlementineStyle(app)
    app.setStyle(style)
    # app.installEventFilter(_TooltipColorFixer(style))

    # Force dot as decimal separator
    qt_locale = QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates)
    QtCore.QLocale.setDefault(qt_locale)

    window = AtomViewerGUI()
    window.show()

    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
