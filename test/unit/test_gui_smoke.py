"""Smoke test for the graphical interface.

Deliberately shallow: it opens the main window and checks that the
actions are wired, nothing more. Simulating clicks would cost far more
than it is worth, but at the start of this suite's life the GUI package
was not even importable -- `tabs.py` imported `plotter_widget` as a
top-level module -- and several Qt overrides were spelled in snake_case
so Qt never called them. Three lines of import would have caught all of
it.

The 3D view is not exercised: `StructureTabWidget` builds a PyVista
`QtInteractor`, which needs a real GL context and dies on a headless
runner. Opening a structure is therefore out of scope here.
"""

from __future__ import annotations

import os

import pytest

# Must be set before the Qt platform plugin is chosen, i.e. before the
# first QApplication is built.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6", reason="the GUI extra is not installed")
pytest.importorskip("pyvistaqt", reason="the GUI extra is not installed")


@pytest.fixture(scope="module")
def qt_app():
    from PySide6 import QtWidgets

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield app


@pytest.fixture(scope="module")
def window(qt_app):
    from cemd.gui.main_window import AtomViewerGUI

    return AtomViewerGUI()


def test_gui_package_is_importable():
    # The whole point: an import error here means the interface cannot
    # start at all, however well the library underneath behaves.
    import cemd.gui.logic.build  # noqa: F401
    import cemd.gui.main_window  # noqa: F401
    import cemd.gui.plotter_widget  # noqa: F401
    import cemd.gui.tabs  # noqa: F401
    import cemd.gui.ui.build  # noqa: F401


def test_entry_point_target_exists():
    from cemd.gui.main_window import main

    assert callable(main)


def test_main_window_opens(window):
    assert window.windowTitle()
    assert window.tabs.count() == 0


def test_every_toolbar_action_has_a_tooltip(window):
    actions = [a for a in window.tools_toolbar.actions() if not a.isSeparator()]
    assert len(actions) > 20

    without = [a.text() for a in actions if not a.toolTip()]
    assert without == [], f"actions with no tooltip: {without}"


def test_file_shortcuts_are_bound(window):
    bound = {
        a.text(): a.shortcut().toString()
        for a in window.tools_toolbar.actions()
        if not a.isSeparator() and not a.shortcut().isEmpty()
    }
    assert bound.get("Open") == "Ctrl+O"
    assert bound.get("Save") == "Ctrl+S"
    assert bound.get("Undo") == "Ctrl+Z"


def test_menu_bar_groups_the_actions(window):
    # The QAction list is kept alive for the whole loop: dropping it and
    # holding only the QMenu it returns lets PySide6 collect the Python
    # wrapper and invalidate the menu underneath.
    menu_actions = list(window.menuBar().actions())
    titles = [m.text() for m in menu_actions]
    assert set(titles) == {
        "&File", "&Edit", "&Build", "&Transform", "&Analyze", "&View"
    }

    # Menu entries are the toolbar's own QAction objects, not copies, so
    # the two can never drift apart.
    toolbar_actions = set(window.tools_toolbar.actions())
    for menu_action in menu_actions:
        entries = [a for a in menu_action.menu().actions() if not a.isSeparator()]
        assert entries, f"{menu_action.text()} is empty"
        assert all(a in toolbar_actions for a in entries), menu_action.text()


def test_tools_are_disabled_until_a_structure_is_open(window):
    window.set_tools_enabled(False)
    assert not window.action_replicate.isEnabled()
    assert not window.action_undo.isEnabled()


def test_undo_on_an_empty_window_is_a_noop(window):
    # No tab, no history: this must not raise.
    window.undo_clicked()


def test_user_files_are_written_outside_the_package(tmp_path, monkeypatch):
    """The interface must not write into its own installation.

    It used to: preferences and the COD/PubChem caches were saved next to
    `main_window.py`. In a source tree that works and hides the problem.
    Installed, that path is `site-packages/cemd/gui` -- shared between the
    users of a machine, replaced on every upgrade, and read-only as often
    as not.
    """
    from pathlib import Path

    from cemd.gui import _userdata

    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))

    package = Path(_userdata.__file__).resolve().parent
    for path in (_userdata.config_file(), _userdata.cache_file("cod_cache.json")):
        assert package not in path.resolve().parents
        assert path.parent.is_dir(), "the directory must be created for the write"
        assert path.parent.name == _userdata.APP_DIR


def test_every_writer_goes_through_userdata():
    """`__file__` may locate shipped data, never a file the GUI writes.

    A tripwire rather than a proof: it catches a new `json.dump` or
    `open(..., "w")` added to a module that has no idea where the user's
    directories are.
    """
    from pathlib import Path

    import cemd.gui

    package = Path(cemd.gui.__file__).parent
    for source in package.rglob("*.py"):
        if source.name == "_userdata.py":
            continue  # it resolves the legacy path on purpose, read-only
        text = source.read_text()
        writes = 'open(' in text and '"w"' in text or "json.dump(" in text
        if writes and "_userdata" not in text:
            raise AssertionError(
                f"{source.relative_to(package)} writes a file but does not use "
                "_userdata; a path built from __file__ lands in site-packages"
            )


# ----------------------------------------------------------------------
# Undo through the type and connectivity managers
#
# The 3D view cannot be built headless (see the module docstring), so these
# use a bare QWidget as the tab: the managers only need `.system`,
# `.session_ff_dict` and the main window's `push_undo` / `sync_ui`.
# ----------------------------------------------------------------------


@pytest.fixture
def water_tab(window, monkeypatch):
    from PySide6 import QtWidgets

    from conftest import make_water_system

    tab = QtWidgets.QWidget()
    tab.system = make_water_system()
    tab.session_ff_dict = {}
    monkeypatch.setattr(window, "sync_ui", lambda *a, **k: None)
    window.tabs.addTab(tab, "water")
    window.tabs.setCurrentWidget(tab)
    yield tab
    window.tabs.removeTab(window.tabs.indexOf(tab))


def test_undo_restores_the_force_field_choices_with_the_system(window, water_tab):
    water_tab.session_ff_dict = {"Ow": {"type": "OW", "model": "spce"}}
    window.push_undo()

    water_tab.session_ff_dict.clear()
    water_tab.system.set_types({"Ow": "O"})
    window.undo_clicked()

    assert "Ow" in water_tab.system.atom_types
    assert water_tab.session_ff_dict == {"Ow": {"type": "OW", "model": "spce"}}


def test_type_manager_edits_can_be_undone(window, water_tab):
    from cemd.gui.ui.managers import TypeManagerDialog

    dialog = TypeManagerDialog(window)
    original_charges = dict(water_tab.system.charges)

    # Rename a type and change a charge through the table, then Apply.
    table = dialog.table
    for row in range(table.rowCount()):
        if table.item(row, 0).text() == "Ow":
            table.item(row, 0).setText("Ox")
            table.item(row, 3).setText("-0.9000")
    dialog.apply_and_refresh()

    assert "Ox" in water_tab.system.atom_types
    assert len(water_tab.undo_stack) == 1

    # Nothing changed the second time: no snapshot is spent on it.
    dialog.apply_and_refresh()
    assert len(water_tab.undo_stack) == 1

    window.undo_clicked()
    assert "Ow" in water_tab.system.atom_types
    assert "Ox" not in water_tab.system.atom_types
    assert dict(water_tab.system.charges) == original_charges


def test_neutralize_and_guess_types_can_be_undone(window, water_tab):
    from cemd.gui.ui.managers import TypeManagerDialog

    dialog = TypeManagerDialog(window)
    system = water_tab.system
    system.set_charges({"Ow": -0.5, "Hw": 0.4})
    charges = dict(system.charges)

    dialog.neutralize_logic({"Ow": 1.0})
    assert abs(system.atoms["charge"].sum()) < 1e-6
    window.undo_clicked()
    assert dict(water_tab.system.charges) == charges

    # An already-neutral system takes no snapshot.
    system = water_tab.system
    system.set_charges({"Ow": -0.8, "Hw": 0.4})
    depth = len(water_tab.undo_stack)
    dialog.neutralize_logic({"Ow": 1.0})
    assert len(water_tab.undo_stack) == depth

    dialog.reset_via_masses()
    assert len(water_tab.undo_stack) == depth + 1


def test_connectivity_edits_can_be_undone(window, water_tab, monkeypatch):
    from PySide6 import QtWidgets

    from cemd.gui.ui.managers import ConnectivityDialog

    dialog = ConnectivityDialog(window)
    n_bonds = len(water_tab.system.bonds)
    bond_type = water_tab.system.bonds["type"].iloc[0]

    monkeypatch.setattr(
        QtWidgets.QMessageBox,
        "question",
        staticmethod(lambda *a, **k: QtWidgets.QMessageBox.Yes),
    )
    dialog.delete_by_type("Bonds", bond_type)
    assert water_tab.system.bonds is None or len(water_tab.system.bonds) < n_bonds
    window.undo_clicked()
    assert len(water_tab.system.bonds) == n_bonds


def test_a_declined_deletion_takes_no_snapshot(window, water_tab, monkeypatch):
    from PySide6 import QtWidgets

    from cemd.gui.ui.managers import ConnectivityDialog

    dialog = ConnectivityDialog(window)
    monkeypatch.setattr(
        QtWidgets.QMessageBox,
        "question",
        staticmethod(lambda *a, **k: QtWidgets.QMessageBox.No),
    )
    dialog.delete_by_type("Bonds", water_tab.system.bonds["type"].iloc[0])
    assert not getattr(water_tab, "undo_stack", [])


# ----------------------------------------------------------------------
# Atom table and side panels
#
# These run without a 3D view: the table, its selection and the panels are
# plain Qt widgets.
# ----------------------------------------------------------------------


def test_atom_table_shows_the_atoms(qt_app):
    from conftest import make_water_system

    from cemd.gui.ui.atom_table import AtomTable

    table = AtomTable()
    table.update_data(make_water_system())
    model = table.model()

    assert model.rowCount() == 6
    assert model.index(0, 0).data() == "1"
    assert model.index(0, 1).data() == "Ow"
    assert model.index(1, 2).data() == "0.9600"


def test_table_selection_follows_the_picked_atoms(qt_app):
    from types import SimpleNamespace

    from conftest import make_water_system

    from cemd.gui.tabs import StructureTabWidget
    from cemd.gui.ui.atom_table import AtomTable

    table = AtomTable()
    table.update_data(make_water_system())
    # Atom ids 1-3 are a contiguous run, 5 stands alone.
    tab = SimpleNamespace(table=table, selected_real_indices=[1, 2, 3, 5])
    StructureTabWidget.sync_table_selection(tab)

    rows = sorted(i.row() for i in table.selectionModel().selectedRows())
    assert rows == [0, 1, 2, 4]


def test_filter_panel_counts_atoms_per_type(qt_app):
    from conftest import make_water_system

    from cemd.gui.ui.panels import FilterPanel

    panel = FilterPanel()
    panel.refresh(make_water_system(), {}, {}, 1.0)
    assert "4 atoms of type Hw" in panel.checkboxes["Hw"].toolTip()


def test_a_color_change_refreshes_the_view_once(qt_app):
    """The color button opens its own dialog; the panel must not add one."""
    from conftest import make_water_system
    from PySide6 import QtGui

    from cemd.gui.ui.panels import FilterPanel

    panel = FilterPanel()
    panel.refresh(make_water_system(), {}, {}, 1.0)
    emitted = []
    panel.type_changed.connect(emitted.append)

    panel.color_buttons["Ow"].setColor(QtGui.QColor("#00ff00"))
    assert len(emitted) == 1
    assert emitted[0]["colors"]["Ow"] == "#00ff00"


def test_slider_moves_refresh_the_view_after_a_burst(qt_app):
    """Any move applies, not just a drag release, and a burst applies once."""
    from conftest import make_water_system

    from cemd.gui.ui.panels import FilterPanel

    panel = FilterPanel()
    panel.refresh(make_water_system(), {}, {}, 1.0)
    emitted = []
    panel.type_changed.connect(emitted.append)

    for value in (110, 120, 130):
        panel.global_scale_slider.setValue(value)
    assert emitted == []

    panel._emit_timer.timeout.emit()
    assert len(emitted) == 1
    assert emitted[0]["global_scale"] == 130


# ----------------------------------------------------------------------
# Long builds off the GUI thread
# ----------------------------------------------------------------------


def test_background_work_returns_its_result_and_keeps_the_gui_alive(qt_app):
    import time

    from PySide6 import QtCore, QtWidgets

    from cemd.gui.logic.build import run_in_background

    parent = QtWidgets.QWidget()
    ticks = []
    timer = QtCore.QTimer()
    timer.timeout.connect(lambda: ticks.append(1))
    timer.start(10)

    def work():
        time.sleep(0.3)
        return 42

    assert run_in_background(parent, "Working...", work) == 42
    timer.stop()
    # The event loop ran while the work slept: the window would repaint.
    assert len(ticks) >= 5


def test_background_work_raises_its_error_in_the_gui_thread(qt_app):
    from PySide6 import QtWidgets

    from cemd.gui.logic.build import run_in_background

    def work():
        raise RuntimeError("Packmol failed")

    with pytest.raises(RuntimeError, match="Packmol failed"):
        run_in_background(QtWidgets.QWidget(), "Working...", work)


def test_background_work_that_finishes_at_once_does_not_hang(qt_app):
    from PySide6 import QtWidgets

    from cemd.gui.logic.build import run_in_background

    parent = QtWidgets.QWidget()
    for _ in range(20):
        assert run_in_background(parent, "Working...", lambda: "ok") == "ok"


# ----------------------------------------------------------------------
# Unsaved changes
# ----------------------------------------------------------------------


def _forbid_message_boxes(monkeypatch):
    from PySide6 import QtWidgets

    def fail(*args, **kwargs):
        raise AssertionError("no dialog expected")

    monkeypatch.setattr(QtWidgets.QMessageBox, "warning", staticmethod(fail))


def _answer_message_boxes(monkeypatch, answer):
    from PySide6 import QtWidgets

    monkeypatch.setattr(
        QtWidgets.QMessageBox,
        "warning",
        staticmethod(lambda *a, **k: getattr(QtWidgets.QMessageBox, answer)),
    )


def test_an_edit_marks_the_tab_modified(window, water_tab):
    window.set_tab_modified(water_tab, False)
    window.push_undo()
    assert water_tab.modified
    assert window.tabs.tabToolTip(window.tabs.indexOf(water_tab)) == "Unsaved changes"
    # The title stays clean: other dialogs list structures by tab title.
    assert window.tabs.tabText(window.tabs.indexOf(water_tab)) == "water"


def test_a_saved_tab_closes_without_asking(window, water_tab, monkeypatch):
    _forbid_message_boxes(monkeypatch)
    window.set_tab_modified(water_tab, False)
    window.request_close_tab(window.tabs.indexOf(water_tab))
    assert window.tabs.indexOf(water_tab) == -1


def test_a_modified_tab_stays_open_on_cancel(window, water_tab, monkeypatch):
    _answer_message_boxes(monkeypatch, "Cancel")
    window.set_tab_modified(water_tab, True)
    window.request_close_tab(window.tabs.indexOf(water_tab))
    assert window.tabs.indexOf(water_tab) != -1


def test_a_modified_tab_closes_on_discard(window, water_tab, monkeypatch):
    _answer_message_boxes(monkeypatch, "Discard")
    window.set_tab_modified(water_tab, True)
    window.request_close_tab(window.tabs.indexOf(water_tab))
    assert window.tabs.indexOf(water_tab) == -1


def test_quitting_with_unsaved_changes_can_be_cancelled(window, water_tab, monkeypatch):
    from PySide6 import QtGui

    _answer_message_boxes(monkeypatch, "Cancel")
    monkeypatch.setattr(
        window, "save_global_to_json", lambda: pytest.fail("must not get this far")
    )
    window.set_tab_modified(water_tab, True)

    event = QtGui.QCloseEvent()
    window.closeEvent(event)
    assert not event.isAccepted()


def test_a_coordinate_typed_in_the_table_reaches_the_system(qt_app):
    from types import SimpleNamespace

    from conftest import make_water_system
    from PySide6 import QtCore

    from cemd.gui.tabs import StructureTabWidget
    from cemd.gui.ui.atom_table import AtomTable

    system = make_water_system()
    table = AtomTable()
    table.update_data(system)
    marked = []
    tab = SimpleNamespace(
        table=table,
        system=system,
        plotter=SimpleNamespace(color_map={}),
        _bonds_key="stale",
        parent_gui=SimpleNamespace(
            _is_syncing=False,
            set_tab_modified=lambda t, m: marked.append(m),
            sync_ui=lambda **k: None,
        ),
    )

    model = table.model()
    cell = model.index(1, 2)  # atom id 2, column x
    model.setData(cell, "4.5", QtCore.Qt.ItemDataRole.EditRole)
    StructureTabWidget.on_data_changed(tab, cell)

    assert system.atoms.loc[2, "x"] == 4.5
    assert tab._bonds_key is None, "the bonds must be searched again"
    assert marked == [True]


# ----------------------------------------------------------------------
# Display settings: one object shared by every tab
# ----------------------------------------------------------------------


@pytest.fixture
def user_config_dir(tmp_path, monkeypatch):
    from cemd.gui import _userdata

    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    # A development tree may still hold the config.json earlier versions
    # wrote into the package, which `load` falls back to.
    monkeypatch.setattr(
        _userdata, "legacy_config_file", lambda: tmp_path / "no-legacy.json"
    )
    assert tmp_path in _userdata.config_file().parents
    return tmp_path / "config"


def test_saved_preferences_survive_a_session_with_no_structure(
    qt_app, user_config_dir
):
    """Opening and closing the window used to reset the scale sliders' values."""
    import json

    from cemd.gui import _userdata, display_settings
    from cemd.gui.main_window import AtomViewerGUI

    saved = {
        "global_scale": 1.5,
        "global_bond_radius": 0.25,
        "bg_color": "white",
        "bond_map": {"Hw-Ow": 1.1},
    }
    _userdata.config_file().write_text(json.dumps(saved))

    window = AtomViewerGUI()
    assert window.filter_panel.get_scale_value() == 1.5
    assert window.bond_manager.get_bond_radius() == 0.25
    window.save_global_to_json()

    reloaded = display_settings.load()
    assert reloaded.global_scale == 1.5
    assert reloaded.global_bond_radius == 0.25
    assert reloaded.bg_color == "white"
    assert reloaded.bond_map == {"Hw-Ow": 1.1}


def test_only_changes_from_the_defaults_are_saved(qt_app, user_config_dir):
    import json

    from cemd.gui import display_settings

    settings = display_settings.load()
    some_type, default_color = next(iter(settings.color_map.items()))
    settings.color_map[some_type] = "#123456"

    path = display_settings.save(settings)
    written = json.loads(path.read_text())
    assert written["color_map"] == {some_type: "#123456"}
    assert default_color != "#123456"


def test_a_deleted_bond_rule_stays_deleted():
    """Rules were merged into the map, so the next panel rebuild revived them."""
    from cemd.gui.display_settings import DisplaySettings

    settings = DisplaySettings(bond_map={"Hw-Ow": 1.1, "Ca-O": 2.6, "Si-O": 1.8})

    # The water structure's panel now holds a single rule, Ow-Ow.
    settings.set_bond_rules(["Ow", "Hw"], {"Ow-Ow": 3.0})

    assert settings.bond_map == {"Ow-Ow": 3.0, "Ca-O": 2.6, "Si-O": 1.8}
    assert settings.bond_rules_for(["Ow", "Hw"]) == {"Ow-Ow": 3.0}
    assert settings.bond_rules_for(["Ca", "O", "Si"]) == {"Ca-O": 2.6, "Si-O": 1.8}


def test_the_background_cycles_and_knows_when_it_is_light():
    from cemd.gui.display_settings import BG_CYCLE, DisplaySettings

    settings = DisplaySettings()
    seen = [settings.bg_color] + [settings.next_background() for _ in BG_CYCLE]
    assert seen[: len(BG_CYCLE)] == list(BG_CYCLE)
    assert seen[-1] == BG_CYCLE[0]

    settings.bg_color = "#B0C4DE"
    assert settings.is_light_background


def test_bond_panel_changes_reach_the_shared_settings(window, water_tab):
    water_tab.refresh_tab_view = lambda **kwargs: None
    window.display.bond_map = {"Hw-Ow": 1.1, "Ca-O": 2.6}

    window.refresh_bonds_view({"Ow-Ow": 3.0, "global_bond_radius": 0.2})

    assert window.display.bond_map == {"Ow-Ow": 3.0, "Ca-O": 2.6}
    assert window.display.global_bond_radius == 0.2


def test_type_panel_changes_reach_the_shared_settings(window, water_tab):
    window._store_view_settings(
        {"colors": {"Ow": "#ff0000"}, "radii": {"Ow": 2.0}, "global_scale": 130}
    )
    assert window.display.color_map["Ow"] == "#ff0000"
    assert window.display.radius_map["Ow"] == 2.0
    assert window.display.global_scale == 1.3
