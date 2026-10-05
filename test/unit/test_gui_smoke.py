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
