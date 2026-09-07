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
