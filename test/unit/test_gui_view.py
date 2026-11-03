"""Tests for the 3D view: opening a structure, and the moving-view LOD.

The PyVista `QtInteractor` needs a GL context, which the offscreen Qt
platform used by test_gui_smoke.py cannot give, so these tests need a
display. CI runs this module on its own, under Xvfb:

    xvfb-run -a env QT_QPA_PLATFORM=xcb pytest test/unit/test_gui_view.py

In a run that also collects test_gui_smoke.py, Qt is set to offscreen
before this module is imported, and the tests skip themselves.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PySide6", reason="the GUI extra is not installed")
pytest.importorskip("pyvistaqt", reason="the GUI extra is not installed")

if os.environ.get("QT_QPA_PLATFORM") == "offscreen" or not (
    os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")
):
    pytest.skip(
        "the 3D view needs a display; run this module on its own, "
        "under xvfb-run if headless",
        allow_module_level=True,
    )


@pytest.fixture(scope="module")
def qt_app():
    from PySide6 import QtWidgets

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    if app.platformName() == "offscreen":
        pytest.skip("Qt already runs on the offscreen platform")
    yield app


@pytest.fixture(scope="module")
def window(qt_app, tmp_path_factory):
    from cemd.gui import _userdata
    from cemd.gui.main_window import AtomViewerGUI

    # The window loads its display settings from the user's config, and
    # the shipped defaults hold no bond rules: a test that relied on the
    # developer's own config.json passed locally and failed on CI. Point
    # it at an empty directory so it sees the defaults only, wherever it
    # runs, and never writes over a real config.
    home = tmp_path_factory.mktemp("home")
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("HOME", str(home))
        mp.setenv("XDG_CONFIG_HOME", str(home / "config"))
        mp.setattr(_userdata, "legacy_config_file", lambda: home / "no-legacy.json")
        window = AtomViewerGUI()
        window.show()
        yield window
        for index in reversed(range(window.tabs.count())):
            window.set_tab_modified(window.tabs.widget(index), False)
        window.close()


def test_opening_a_structure_draws_it(window):
    from conftest import make_water_system

    # The defaults define no bond rules, so give the one this system needs.
    window.display.bond_map["Hw-Ow"] = 1.1
    window.add_structure_tab(make_water_system(), title="water")
    tab = window.tabs.currentWidget()
    actors = tab.plotter.renderer.actors

    assert {"atoms_Ow", "atoms_Hw"} <= set(actors)
    assert "bonds" in actors
    # Six atoms: far below the budget, so no moving-view copy.
    assert not any(name.startswith("lod_") for name in actors)


@pytest.fixture
def plotter(qt_app, monkeypatch):
    from cemd.gui.plotter_widget import AtomicPlotter

    monkeypatch.setattr(AtomicPlotter, "LOD_ATOM_BUDGET", 100)
    monkeypatch.setattr(AtomicPlotter, "LOD_KEEP_TYPES_UP_TO", 50)
    plotter = AtomicPlotter()
    plotter.show()
    rng = np.random.default_rng(0)
    xyz = rng.random((1030, 3)) * 30.0
    df = pd.DataFrame(xyz, columns=["x", "y", "z"])
    # 1000 "O" (subsampled while moving) and 30 "Na" (kept whole).
    df["type"] = ["O"] * 1000 + ["Na"] * 30
    plotter.draw_atoms(df)
    plotter.add_bonds_smart(xyz, df["type"].values, {"O-O": 1.5})
    plotter.render()
    yield plotter
    plotter.close()


def _visible(plotter):
    return {name for name, a in plotter.renderer.actors.items() if a.visibility}


def test_a_large_system_gets_a_hidden_moving_view_copy(plotter):
    actors = plotter.renderer.actors
    # 1030 atoms for a budget of 100: every 11th atom.
    assert actors["lod_atoms_O"].GetMapper().GetInput().n_points == 91
    assert "lod_atoms_Na" not in actors
    assert {"atoms_O", "atoms_Na", "bonds"} <= _visible(plotter)
    assert "lod_atoms_O" not in _visible(plotter)


def test_a_rotation_draws_the_copy_and_restores_the_full_scene(plotter):
    style = plotter.iren.style
    style.SetCurrentRenderer(plotter.renderer)

    # What a mouse drag does: StartRotate raises the desired update rate,
    # every frame of the drag renders with it, EndRotate drops it back and
    # renders once more.
    style.StartRotate()
    plotter.render()
    moving = _visible(plotter)
    style.EndRotate()

    assert "lod_atoms_O" in moving
    assert "atoms_Na" in moving
    assert "atoms_O" not in moving
    assert "bonds" not in moving

    still = _visible(plotter)
    assert {"atoms_O", "atoms_Na", "bonds"} <= still
    assert "lod_atoms_O" not in still


def test_a_hidden_type_stays_hidden_through_a_rotation(plotter):
    plotter.update_atom_visibility({"Na"})
    style = plotter.iren.style
    style.SetCurrentRenderer(plotter.renderer)

    style.StartRotate()
    plotter.render()
    moving = _visible(plotter)
    style.EndRotate()

    assert "lod_atoms_O" not in moving
    assert "atoms_O" not in _visible(plotter)


def test_a_redraw_during_a_rotation_is_not_hidden_after_it(plotter):
    style = plotter.iren.style
    style.SetCurrentRenderer(plotter.renderer)

    style.StartRotate()
    plotter.render()
    # A refresh replaces the actors while the old ones are hidden.
    df = pd.DataFrame(
        np.random.default_rng(1).random((1030, 3)) * 30.0, columns=["x", "y", "z"]
    )
    df["type"] = "O"
    plotter.draw_atoms(df)
    style.EndRotate()

    assert "atoms_O" in _visible(plotter)
    assert "lod_atoms_O" not in _visible(plotter)


def test_size_and_color_changes_reach_the_moving_view_copy(plotter):
    plotter.settings.color_map["O"] = "#00ff00"
    plotter.settings.radius_map["O"] = 2.0
    plotter.update_atom_colors()
    plotter.update_atom_sizes()

    full = plotter.renderer.actors["atoms_O"].prop
    lod = plotter.renderer.actors["lod_atoms_O"].prop
    assert lod.color == full.color
    assert lod.point_size == full.point_size


def test_the_split_dialog_previews_the_cut(window):
    from cemd.gui.ui.build import SplitterDialog

    from conftest import make_water_system

    window.add_structure_tab(make_water_system(), title="water to split")
    dialog = SplitterDialog(window)
    actors = dialog.view.renderer.actors

    # It opens on the best cut, which breaks no O-H bond, and draws it.
    assert dialog.scan is not None
    assert dialog.coord.value() == dialog.best
    assert "cut_plane" in actors and "cut_bonds" not in actors

    # Along a, between an O (x=0) and its H (x=0.96): the bond is counted
    # and drawn in red.
    dialog.axis.setCurrentIndex(0)
    dialog.coord.setValue(0.5)
    dialog.update_preview()
    expected = dialog.splitter.count_broken_bonds(0.5)
    assert expected > 0
    assert dialog.severed.text().startswith(f"<b>{expected}</b> bond")
    assert "cut_bonds" in dialog.view.renderer.actors
    assert dialog.slider.value() == round(dialog.coord.value() * 100)

    # Kept whole, the water molecule crosses the plane instead: no bond is
    # severed, and its atoms on the wrong side of the plane are haloed.
    dialog.cut_mode.setCurrentIndex(2)
    dialog.coord.setValue(0.5)
    dialog.update_preview()
    actors = dialog.view.renderer.actors
    assert dialog.severed.text().startswith("<b>0</b> bonds")
    assert "atoms off the plane" in dialog.severed.text()
    assert "cut_bonds" not in actors and "cut_shifted" in actors
    assert dialog.get_values()["keep_molecules"] is True

    # Repaired: the O-H stays whole (the H goes with its O), nothing to
    # rebuild here.
    dialog.cut_mode.setCurrentIndex(1)
    dialog.coord.setValue(0.5)
    dialog.update_preview()
    assert dialog.severed.text().startswith("<b>0</b> bonds to rebuild")
    assert dialog.get_values()["repair"] is True

    dialog.reject()
    assert dialog.view is None
