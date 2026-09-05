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

from typing import TYPE_CHECKING

from PySide6 import QtCore, QtWidgets

from cemd.build import (
    GlassBuilder,
    SolutionBuilder,
    Splitter,
)
from cemd.gui.ui.build import (
    AddDropletDialog,
    AddLiquidLayerDialog,
    AddStructureDialog,
    CASHBuilderDialog,
    GlassBuilderDialog,
    ReplicateDialog,
    SmilesDialog,
    SolutionDialog,
    SplitterDialog,
    SurfaceDialog,
    TranslateAtomsDialog,
)

if TYPE_CHECKING:
    from cemd.core.atomic_system import AtomicSystem

    from ..main_window import AtomViewerGUI


def handle_error(parent: AtomViewerGUI, title: str, error: Exception | str) -> None:
    QtWidgets.QApplication.restoreOverrideCursor()
    parent.setUpdatesEnabled(True)
    QtWidgets.QMessageBox.critical(parent, title, f"An error occurred:\n{error}")


def open_make_solution(parent: AtomViewerGUI) -> None:
    dialog = SolutionDialog(parent)
    if dialog.exec_():
        try:
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            box, density, solutes, structures = dialog.get_values()
            builder = SolutionBuilder(
                density=density, counts=solutes, structures=structures
            )
            new_data = builder.build(box)
            parent.add_structure_tab(new_data, title="Aqueous Solution")
        except Exception as e:
            handle_error(parent, "Solution Error", e)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()


def open_make_glass(parent: AtomViewerGUI) -> None:
    dialog = GlassBuilderDialog(parent)
    if dialog.exec() == QtWidgets.QDialog.Accepted:
        box, density, stoich = dialog.get_values()

        try:
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            parent.statusBar().showMessage("Generating glass structure with Packmol...")

            builder = GlassBuilder(density=density, composition=stoich)
            new_system = builder.build(box)

            if new_system:
                parent.add_structure_tab(new_system, f"Glass_{density}gcm3")

        except Exception as e:
            QtWidgets.QMessageBox.critical(parent, "Error", f"Generation failed: {e}")
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
            parent.statusBar().showMessage("Ready")


def open_make_cash(parent: AtomViewerGUI) -> None:
    dialog = CASHBuilderDialog(parent)
    if dialog.exec_():
        system = dialog.get_system()
        if system:
            try:
                QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

                nsi = system.get_count("Si")
                nca = system.get_count("Ca")

                system.set_topology("cshff")

                ratio_name = f"{nca / nsi:.2f}" if nsi > 0 else "Custom"

                parent.add_structure_tab(system, title=f"C-S-H {ratio_name}")
                parent.statusBar().showMessage(
                    f"C-(A)-S-H (Ca/Si: {ratio_name}) model imported!", 5000
                )

            except Exception as e:
                handle_error(parent, "Import Error", e)
            finally:
                QtWidgets.QApplication.restoreOverrideCursor()


def open_make_surface(parent: AtomViewerGUI) -> None:
    active_tab = parent.tabs.currentWidget()
    if not active_tab or not hasattr(active_tab, "system"):
        return QtWidgets.QMessageBox.warning(
            parent, "Selection Error", "Please select a solid system first."
        )

    dialog = SurfaceDialog(active_tab.system, parent=parent)
    if dialog.exec_() == QtWidgets.QDialog.Accepted:
        systems = dialog.selected_systems
        if not systems:
            return
        try:
            parent.setUpdatesEnabled(False)
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

            miller_indices = dialog.get_values()[0]
            m_str = "".join(map(str, miller_indices))

            for i, s in enumerate(systems):
                parent.add_structure_tab(s, title=f"slab_{m_str}_{i + 1}")

            parent.sync_ui(full_rebuild=True, reset_camera=True)
            parent.statusBar().showMessage(f"Generated {len(systems)} surfaces.", 5000)
        except Exception as e:
            handle_error(parent, "Surface Error", e)
        finally:
            parent.setUpdatesEnabled(True)
            QtWidgets.QApplication.restoreOverrideCursor()


def open_add_liquid(parent: AtomViewerGUI) -> None:
    active_tab = parent.tabs.currentWidget()
    if not active_tab or not hasattr(active_tab, "system"):
        return

    dialog = AddLiquidLayerDialog(parent)
    if dialog.exec() == QtWidgets.QDialog.Accepted:
        p = dialog.get_values()
        print(p)
        try:
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            parent.statusBar().showMessage("Creating interface...")

            builder = SolutionBuilder(
                density=p["density"],
                counts=p["solutes_dict"],
                structures=p["structures_dict"],
            )
            parent.push_undo()
            active_tab.system.add_liquid_layer(
                blueprint=builder,
                thickness=p["thickness"],
                distance=2.0,
                vacuum=p["vacuum"],
                axis=p["axis"],
            )

            parent.sync_ui(full_rebuild=True, reset_camera=True)
            parent.statusBar().showMessage("Interface created!", 5000)
        except Exception as e:
            handle_error(parent, "Interface Error", e)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()


def open_add_droplet(parent: AtomViewerGUI) -> None:
    active_tab = parent.tabs.currentWidget()
    if not active_tab or not hasattr(active_tab, "system"):
        return

    dialog = AddDropletDialog(parent)
    if dialog.exec() == QtWidgets.QDialog.Accepted:
        p = dialog.get_values()
        try:
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            parent.statusBar().showMessage("Adding droplet...")
            builder = SolutionBuilder(
                density=p["density"],
                counts=p["solutes_dict"],
                structures=p["structures_dict"],
            )
            parent.push_undo()
            active_tab.system.add_droplet(
                blueprint=builder,
                radius=p["radius"],
                distance=2.0,
                vacuum=p["vacuum"],
                axis="z",
            )

            parent.sync_ui(full_rebuild=True, reset_camera=True)
            parent.statusBar().showMessage("Droplet added successfully!", 5000)
        except Exception as e:
            QtWidgets.QMessageBox.critical(parent, "Droplet Error", str(e))
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()


def open_add_structure(parent: AtomViewerGUI) -> None:
    active_tab = parent.tabs.currentWidget()
    if not active_tab or not hasattr(active_tab, "system"):
        return

    dialog = AddStructureDialog(parent)
    if dialog.exec() == QtWidgets.QDialog.Accepted:
        p = dialog.get_values()
        if p is None:
            return
        try:
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            parent.statusBar().showMessage("Adding structure to surface...")

            parent.push_undo()
            active_tab.system.add_structure(
                structure_to_add=p["structure_to_add"],
                distance=p["distance"],
                axis=p["axis"],
                vacuum=p["vacuum"],
            )
            parent.sync_ui(full_rebuild=True, reset_camera=True)
            parent.statusBar().showMessage("Structure added successfully!", 5000)
        except Exception as e:
            QtWidgets.QMessageBox.critical(parent, "Add Structure Error", str(e))
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()


def open_split(parent: AtomViewerGUI) -> None:
    active_tab = parent.tabs.currentWidget()
    if not active_tab or not hasattr(active_tab, "system"):
        return

    dialog = SplitterDialog(parent)
    if dialog.exec() == QtWidgets.QDialog.Accepted:
        p = dialog.get_values()
        try:
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            parent.statusBar().showMessage("Creating channel...")
            parent.push_undo()

            splitter = Splitter(
                active_tab.system,
                coordinate=p["coordinate"],
                axis=p["axis"],
                gap_size=p["gap_size"],
            )

            if p["add_solution"]:
                blueprint = SolutionBuilder(
                    density=p["density"],
                    counts=p["solutes_dict"],
                    structures=p["structures_dict"],
                )
                splitter.add_solution(blueprint, padding=p["tolerance"])

            active_tab.system = splitter.split()

            parent.sync_ui(full_rebuild=True, reset_camera=True)
            parent.statusBar().showMessage("Channel created!", 5000)
        except Exception as e:
            handle_error(parent, "Channel Error", e)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()


def open_replicate(parent: AtomViewerGUI) -> None:
    active_tab = parent.tabs.currentWidget()
    if not active_tab or not hasattr(active_tab, "system"):
        return

    dialog = ReplicateDialog(parent)
    if dialog.exec_():
        try:
            factors = dialog.get_values()
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            parent.push_undo()
            active_tab.system.replicate(factors)
            parent.sync_ui(full_rebuild=True, reset_camera=True)
            parent.statusBar().showMessage(f"Replicated {factors}", 3000)
        except Exception as e:
            handle_error(parent, "Replication Error", e)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()


def on_protonate(parent: AtomViewerGUI) -> None:
    active_tab = parent.tabs.currentWidget()
    if not active_tab or not active_tab.system:
        return

    indices = active_tab.selected_real_indices.copy()
    if not indices:
        return

    try:
        parent.push_undo()
        system: AtomicSystem = active_tab.system

        system.protonate_atoms(indices)

        active_tab.system = system

        parent._is_syncing = False

        parent.sync_ui(full_rebuild=True, reset_camera=False)

    except Exception as e:
        print(f"Erreur Protonate: {e}")
    finally:
        parent._is_syncing = False
        active_tab.selected_real_indices = []
        parent.update_protonate_state()


def open_smiles_builder(parent: AtomViewerGUI) -> None:
    dialog = SmilesDialog(parent)
    if dialog.exec() == QtWidgets.QDialog.Accepted:
        new_system = dialog.system
        parent.add_structure_tab(new_system)


def open_translate_atoms(parent: AtomViewerGUI) -> None:
    active_tab = parent.tabs.currentWidget()
    if not active_tab or not hasattr(active_tab, "system"):
        return

    dialog = TranslateAtomsDialog(parent)
    if dialog.exec() == QtWidgets.QDialog.Accepted:
        dx, dy, dz = dialog.get_values()

        if all(v == 0.0 for v in [dx, dy, dz]):
            return

        try:
            parent.statusBar().showMessage("Translating and wrapping atoms...")

            active_tab.system.atoms["x"] += dx
            active_tab.system.atoms["y"] += dy
            active_tab.system.atoms["z"] += dz

            # Wrap is always applied after a translation
            active_tab.system.wrap()

            # Don't reset the camera zoom
            parent.sync_ui(full_rebuild=True, reset_camera=False)
            parent.statusBar().showMessage("Translation successful!", 3000)

        except Exception as e:
            QtWidgets.QMessageBox.critical(parent, "Translation Error", str(e))
