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

import webbrowser
from typing import Any

import pandas as pd
from PySide6 import QtCore, QtGui, QtWidgets
from PySide6Qlementine import ActionButton

from ...forcefield import ForceFieldDatabase
from ...topology import NeighborCriterion, TopologyRule
from .base_dialog import BaseBuilderDialog
from .gui_utils import get_icon


class FFComboDelegate(QtWidgets.QStyledItemDelegate):
    def __init__(self, parent, ff_db) -> None:
        super().__init__(parent)
        self.ff_db: Any = ff_db

    def createEditor(
        self,
        parent: QtWidgets.QWidget,
        option: QtWidgets.QStyleOptionViewItem,
        index: QtCore.QModelIndex,
    ) -> QtWidgets.QWidget:

        combo = QtWidgets.QComboBox(parent)
        combo.addItem("None", None)

        row_data = index.data(QtCore.Qt.ItemDataRole.UserRole)
        if not row_data:
            return combo

        elem: str = row_data.get("element", "X")
        mass: float = row_data.get("mass", 1.0)

        detected_elem: str = str(elem).strip().upper()
        if mass > 10.0 and detected_elem == "H":
            detected_elem = "O"

        if "element" in self.ff_db.columns:
            mask = self.ff_db["element"].astype(str).str.upper() == detected_elem
            models_sorted = self.ff_db[mask].sort_values(by="model")

            last_model = None
            for _, row in models_sorted.iterrows():
                current_model = str(row.get("model", ""))

                if last_model is not None and current_model != last_model:
                    combo.insertSeparator(combo.count())

                label = f"{row['type']} [{current_model}]"
                combo.addItem(label, row.to_dict())

                if "environment" in row:
                    combo.setItemData(
                        combo.count() - 1,
                        str(row["environment"]),
                        QtCore.Qt.ItemDataRole.ToolTipRole,
                    )

                last_model: str = current_model

        return combo

    def setEditorData(
        self, editor: QtWidgets.QWidget, index: QtCore.QModelIndex
    ) -> None:
        value = index.data(QtCore.Qt.ItemDataRole.DisplayRole)
        idx = editor.findText(value)
        if idx >= 0:
            editor.setCurrentIndex(idx)

    def setModelData(
        self,
        editor: QtWidgets.QWidget,
        model: QtCore.QAbstractItemModel,
        index: QtCore.QModelIndex,
    ) -> None:
        model.setData(index, editor.currentText(), QtCore.Qt.ItemDataRole.DisplayRole)
        model.setData(index, editor.currentData(), QtCore.Qt.ItemDataRole.UserRole)


class TypeManagerDialog(BaseBuilderDialog):
    def __init__(self, parent=None):
        super().__init__(parent, title="Types Manager")
        self.resize(800, 550)

        self.load_ff_database()
        self.setup_ui()

    @property
    def system_obj(self) -> None:
        """Dynamically retrieves the active tab's system via parent."""
        if self.parent() and hasattr(self.parent(), "system"):
            return self.parent().system
        return None

    @property
    def session_ff_dict(self) -> None:
        active_tab = self.parent().tabs.currentWidget()
        if not hasattr(active_tab, "session_ff_dict"):
            active_tab.session_ff_dict = {}
        return active_tab.session_ff_dict

    def load_ff_database(self) -> None:
        """Loads the Excel file and its different tabs"""
        try:
            self.all_sheets = ForceFieldDatabase().to_dataframes()
            self.ff_db = self.all_sheets["atoms"]
            self.ff_db.columns = self.ff_db.columns.str.strip()
        except Exception as e:
            print(f"Error while reading the forcefield database : {e}")
            self.ff_db = pd.DataFrame()

    def setup_ui(self) -> None:
        tools_layout = QtWidgets.QHBoxLayout()

        self.btn_reset = self.create_icon_button("Guess types from masses", "question")
        self.btn_reset.clicked.connect(self.reset_via_masses)

        self.btn_neutralize = self.create_icon_button("Neutralize system", "zap")
        self.btn_neutralize.clicked.connect(self.open_neutralize_dialog)

        self.btn_apply = self.create_icon_button("Apply changes", "save", primary=True)
        self.btn_apply.clicked.connect(self.apply_and_refresh)

        self.btn_copy_coeffs = self.create_icon_button(
            "Copy pair coeff. to clipboard", "paste"
        )
        self.btn_copy_coeffs.clicked.connect(self.copy_pair_coeffs_to_clipboard)

        tools_layout.addWidget(self.btn_reset)
        tools_layout.addWidget(self.btn_neutralize)
        tools_layout.addWidget(self.btn_apply)
        tools_layout.addWidget(self.btn_copy_coeffs)
        self.main_layout.addLayout(tools_layout)

        self.lbl_alert = QtWidgets.QLabel("")
        self.lbl_alert.setWordWrap(True)
        self.lbl_alert.setVisible(False)
        self.main_layout.addWidget(self.lbl_alert)

        self.table = self.create_table(
            ["Type", "FF Parameters", "Mass", "Charge", "Ref."],
            editable=True,
            selectable=True,
        )
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QtWidgets.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QtWidgets.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(4, QtWidgets.QHeaderView.Fixed)
        self.table.setColumnWidth(4, 50)

        self.table.setItemDelegateForColumn(1, FFComboDelegate(self.table, self.ff_db))
        self.table.itemChanged.connect(self.on_item_changed)

        self.table.cellClicked.connect(self.on_cell_clicked)

        self.main_layout.addWidget(self.table, stretch=1)

        self.fill_table()

        self.dialog_buttons = self.create_dialog_buttons(ok_text="Accept All")
        self.dialog_buttons.accepted.connect(self.accept_changes)
        self.dialog_buttons.rejected.connect(self.reject)
        self.main_layout.addWidget(self.dialog_buttons)

    def fill_table(self) -> None:
        """Populate the table and restore saved session parameters."""
        types = self.system_obj.atom_types
        elements_list = self.system_obj.elements

        self.table.setRowCount(len(types))

        for i, t_name in enumerate(types):
            m_val = self.system_obj.masses[t_name]
            c_val = self.system_obj.charges[t_name]
            elem = elements_list[t_name]

            item_id = QtWidgets.QTableWidgetItem(str(t_name))
            item_id.setData(QtCore.Qt.ItemDataRole.UserRole, t_name)
            self.table.setItem(i, 0, item_id)

            ff_item = QtWidgets.QTableWidgetItem("None")
            ff_item.setData(
                QtCore.Qt.ItemDataRole.UserRole, {"element": elem, "mass": m_val}
            )
            self.table.setItem(i, 1, ff_item)

            self.table.setItem(i, 2, QtWidgets.QTableWidgetItem(f"{float(m_val):.5f}"))
            self.table.setItem(i, 3, QtWidgets.QTableWidgetItem(f"{float(c_val):.4f}"))

            ref_item = QtWidgets.QTableWidgetItem()
            ref_item.setFlags(
                QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable
            )
            self.table.setItem(i, 4, ref_item)

            # Restore saved FF
            if t_name in self.session_ff_dict:
                saved_row = self.session_ff_dict[t_name]
                label = f"{saved_row['type']} [{saved_row['model']}]"
                ff_item.setText(label)
                ff_item.setData(QtCore.Qt.ItemDataRole.UserRole, saved_row)
                self.on_ff_changed(i)

        for row in range(self.table.rowCount()):
            self.table.openPersistentEditor(self.table.item(row, 1))

    def on_item_changed(self, item: QtWidgets.QTableWidgetItem) -> None:
        if item.column() == 1:
            self.on_ff_changed(item.row())

    def reset_via_masses(self) -> None:
        self.parent().push_undo()
        self.system_obj.set_types_from_elements()
        self.fill_table()

    def apply_and_refresh(self) -> None:
        self.apply_logic()
        self.fill_table()

    def apply_logic(self) -> None:
        self.table.blockSignals(True)
        temp_names = {}
        temp_charges = {}
        new_session_map = {}

        for i in range(self.table.rowCount()):
            name_item = self.table.item(i, 0)
            if not name_item:
                continue

            original_id = name_item.data(QtCore.Qt.ItemDataRole.UserRole)
            new_name = name_item.text().strip()
            temp_names[original_id] = new_name

            combo = self.table.cellWidget(i, 1)
            ff_data = combo.currentData() if combo else None

            if ff_data:
                new_session_map[new_name] = ff_data
                q_val = float(ff_data.get("charge", 0.0))
            else:
                charge_item = self.table.item(i, 3)
                q_val = float(charge_item.text()) if charge_item else 0.0

            temp_charges[original_id] = q_val

            active_tab = self.parent().tabs.currentWidget()
            color_map = {}
            radius_map = {}
            if hasattr(active_tab, "plotter"):
                color_map = active_tab.plotter.color_map
                radius_map = active_tab.plotter.radius_map

            for original_id, new_name in temp_names.items():
                if original_id != new_name:
                    if original_id in color_map:
                        color_map[new_name] = color_map[original_id]
                    if original_id in radius_map:
                        radius_map[new_name] = radius_map[original_id]

        # Apply and Accept All both land here; only snapshot when the table
        # differs from the system, so a second click does not eat a slot of
        # the undo history.
        system = self.system_obj
        if (
            any(temp_names[t] != t for t in temp_names)
            or any(temp_charges[t] != system.charges[t] for t in temp_charges)
            or new_session_map != self.session_ff_dict
        ):
            self.parent().push_undo()

        # `set_types`/`set_charges` take a {old_value: new_value} mapping,
        # not a bare list of new values -- `temp_names`/`temp_charges` are
        # already built as such dicts above.
        self.system_obj.set_types(temp_names)
        self.system_obj.set_charges(temp_charges)

        self.session_ff_dict.clear()
        self.session_ff_dict.update(new_session_map)
        self.sync_forcefield_parameters()

        self.table.blockSignals(False)

    def sync_forcefield_parameters(self) -> None:
        assignments = {}
        for atom_label, ff_info in self.session_ff_dict.items():
            if ff_info and "type" in ff_info:
                assignments[atom_label] = ff_info["type"]

        if not assignments:
            return

        try:
            self.system_obj.set_ff_from_database(atom_assignments=assignments)
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                self,
                "Error ForceField",
                f"Impossible to assign parameters automatically :\n{str(e)}",
            )

    def on_cell_clicked(self, row: int, col: int) -> None:
        if col == 4:
            item = self.table.item(row, col)
            url = item.data(QtCore.Qt.ItemDataRole.UserRole)
            if url:
                QtGui.QDesktopServices.openUrl(QtCore.QUrl(url))

    def accept_changes(self) -> None:
        self.apply_logic()
        self.accept()

    def on_ff_changed(self, row_idx: int) -> None:
        if self.table.cellWidget(row_idx, 1) is None:
            return
        combo = self.table.cellWidget(row_idx, 1)
        data = combo.currentData()
        charge_item = self.table.item(row_idx, 3)

        # Delete the old widget if it exists to start from scratch
        self.table.removeCellWidget(row_idx, 4)
        # Recreate an empty item for the background/style just in case
        ref_item = QtWidgets.QTableWidgetItem("")
        self.table.setItem(row_idx, 4, ref_item)

        if data:
            # Load management
            charge_col = [c for c in data.keys() if "charge" in c]
            if charge_col:
                charge_item.setText(f"{float(data[charge_col[0]]):.4f}")

            env = data.get("environment", "No description available")
            combo.setToolTip(f"<b>Chemical environment :</b><br>{env}")

            # Reference button
            url = data.get("ref")
            if isinstance(url, str) and url.startswith("http"):
                # Creating the container to center the button
                container = QtWidgets.QWidget()
                btn_layout = QtWidgets.QHBoxLayout(container)
                btn_layout.setContentsMargins(0, 0, 0, 0)
                btn_layout.setAlignment(QtCore.Qt.AlignCenter)

                btn_ref = self.create_icon_button("", "reference")
                btn_ref.setIconSize(QtCore.QSize(18, 18))
                btn_ref.setFixedSize(24, 24)
                btn_ref.setFlat(True)
                btn_ref.setCursor(QtCore.Qt.PointingHandCursor)
                btn_ref.setToolTip(f"Open reference:\n{url}")

                btn_ref.clicked.connect(lambda checked=False, u=url: webbrowser.open(u))

                btn_layout.addWidget(btn_ref)
                self.table.setCellWidget(row_idx, 4, container)
            else:
                ref_item.setText("-")
                ref_item.setTextAlignment(QtCore.Qt.AlignCenter)

            self.check_cross_interactions()
        else:
            combo.setToolTip("No model selected")
            ref_item.setText("-")
            ref_item.setTextAlignment(QtCore.Qt.AlignCenter)

    def check_cross_interactions(self) -> None:
        if not hasattr(self, "lbl_alert"):
            return
        df_lj = self.all_sheets.get("lj_12-6")
        if df_lj is None:
            self.lbl_alert.setVisible(False)
            return

        selected_ff_types = []
        for i in range(self.table.rowCount()):
            combo = self.table.cellWidget(i, 1)
            if combo:
                data = combo.currentData()
                if data and "type" in data:
                    selected_ff_types.append(data["type"])

        unique_selected = list(set(selected_ff_types))
        if len(unique_selected) < 2:
            return

        found_cross_terms = []
        import itertools

        for t1, t2 in itertools.combinations(unique_selected, 2):
            mask = ((df_lj["type 1"] == t1) & (df_lj["type 2"] == t2)) | (
                (df_lj["type 1"] == t2) & (df_lj["type 2"] == t1)
            )
            if not df_lj[mask].empty:
                found_cross_terms.append(f"{t1} ↔ {t2}")

        if found_cross_terms:
            self.lbl_alert.setText(f"⚠ Cross-terms: {', '.join(found_cross_terms)}")
            self.lbl_alert.setVisible(True)
        else:
            self.lbl_alert.setVisible(False)

    def copy_pair_coeffs_to_clipboard(self) -> None:
        df_lj = self.all_sheets.get("lj_12-6")
        if df_lj is None:
            return

        label_to_ff = {}
        for i in range(self.table.rowCount()):
            name_item = self.table.item(i, 0)
            combo = self.table.cellWidget(i, 1)
            if name_item and combo:
                atom_label = name_item.text().strip()
                data = combo.currentData()
                if data and "type" in data:
                    label_to_ff[atom_label] = data["type"]

        if not label_to_ff:
            QtWidgets.QMessageBox.warning(
                self, "Export", "No Force Field types assigned to atoms."
            )
            return

        lines = ["# --- LAMMPS Pair Coefficients (using atom labels) ---"]
        labels = list(label_to_ff.keys())
        count = 0
        import itertools

        for label1, label2 in itertools.combinations_with_replacement(labels, 2):
            ff_t1 = label_to_ff[label1]
            ff_t2 = label_to_ff[label2]
            mask = ((df_lj["type 1"] == ff_t1) & (df_lj["type 2"] == ff_t2)) | (
                (df_lj["type 1"] == ff_t2) & (df_lj["type 2"] == ff_t1)
            )
            row = df_lj[mask]
            if not row.empty:
                eps = float(row.iloc[0].iloc[2])
                sig = float(row.iloc[0].iloc[3])
                line = f"pair_coeff {label1} {label2} {eps:.6f} {sig:.6f}  # ({ff_t1} - {ff_t2})"
                lines.append(line)
                count += 1

        full_text = "\n".join(lines)
        clipboard = QtWidgets.QApplication.clipboard()
        clipboard.setText(full_text)
        QtWidgets.QMessageBox.information(
            self, "Export Successful", f"Copied {count} interaction lines."
        )

    def open_neutralize_dialog(self) -> None:
        """Opens the selection popup and applies the override."""
        if not self.system_obj:
            return

        types = self.system_obj.atom_types
        diag = NeutralizeDialog(types, self)

        if diag.exec() == QtWidgets.QDialog.Accepted:
            weights = diag.get_selected_weights()
            if not weights:
                QtWidgets.QMessageBox.warning(
                    self, "Neutralize", "Please select at least one atom type."
                )
                return

            try:
                self.neutralize_logic(weights)
                self.fill_table()
                QtWidgets.QMessageBox.information(
                    self, "Success", "System neutralized successfully."
                )
            except Exception as e:
                QtWidgets.QMessageBox.critical(
                    self, "Error", f"Neutralization failed: {e}"
                )

    def neutralize_logic(self, weights: dict[str, float]) -> None:
        """Built-in version of neutralization function using set_charges and get_count."""
        system = self.system_obj
        df_atoms = system.atoms
        total_charge = df_atoms["charge"].sum()

        if abs(total_charge) < 1e-6:
            return

        total_weight_atoms = 0
        for atype, w in weights.items():
            num_atoms = system.get_count(atype)
            total_weight_atoms += w * num_atoms

        if total_weight_atoms == 0:
            raise ValueError("Selected atom types are not present in the system.")

        charge_delta = -total_charge / total_weight_atoms

        new_charges = {}
        for atype in system.atom_types:
            current_q = system._charges.get(atype, 0.0)
            if atype in weights:
                new_charges[atype] = current_q + (weights[atype] * charge_delta)
            else:
                new_charges[atype] = current_q

        self.parent().push_undo()
        system.set_charges(new_charges)


class ConnectivityDialog(BaseBuilderDialog):
    def __init__(self, parent=None):
        # Base class initialization
        super().__init__(parent, title="LAMMPS Topology Manager")

        self.resize(700, 600)

        self.load_ff_database()
        self.setup_ui()

    @property
    def system_obj(self):
        """Dynamically retrieves the active tab's system via parent."""
        if self.parent() and hasattr(self.parent(), "system"):
            return self.parent().system
        return None

    @property
    def session_ff_dict(self):
        active_tab = self.parent().tabs.currentWidget()
        if not hasattr(active_tab, "session_ff_dict"):
            active_tab.session_ff_dict = {}
        return active_tab.session_ff_dict

    @property
    def plotter(self):
        """Retrieves the plotter (for colors) of the active tab."""
        active_tab = self.parent().tabs.currentWidget()
        return getattr(active_tab, "plotter", None)

    def load_ff_database(self) -> None:
        """Loads the Excel file and its different tabs"""
        try:
            self.all_sheets = ForceFieldDatabase().to_dataframes()
            self.ff_db = self.all_sheets["atoms"]
            self.ff_db.columns = self.ff_db.columns.str.strip()
        except Exception as e:
            print(f"Error while reading the forcefield database : {e}")
            self.ff_db = pd.DataFrame()

    def setup_ui(self) -> None:
        """Dynamically rebuilds tabs and buttons"""

        current_idx = 0
        if hasattr(self, "tabs"):
            current_idx = self.tabs.currentIndex()

        protected_widgets = [
            getattr(self, "btn_build_topology", None),
            getattr(self, "btn_copy_topo", None),
            getattr(self, "dialog_buttons", None),
        ]

        while self.main_layout.count() > 0:
            item = self.main_layout.takeAt(0)
            widget = item.widget()
            if widget:
                if widget in protected_widgets:
                    widget.setParent(None)
                else:
                    widget.deleteLater()
            elif item.layout():
                self.clear_layout(item.layout())

        if not hasattr(self, "btn_build_topology"):
            action = QtGui.QAction(
                get_icon("plus-square-white"), "Add a topology rule", self
            )
            action.triggered.connect(self.open_rule_builder)
            self.btn_build_topology = ActionButton(self)
            self.btn_build_topology.setAction(action)

        if not hasattr(self, "btn_copy_topo"):
            self.btn_copy_topo = self.create_icon_button(
                "Copy Bond/Angle Coeffs", "copy"
            )
            self.btn_copy_topo.clicked.connect(self.copy_topology_coeffs_to_clipboard)

        # Added button layout AT TOP
        top_buttons_layout = QtWidgets.QHBoxLayout()
        top_buttons_layout.addWidget(self.btn_build_topology)
        top_buttons_layout.addWidget(self.btn_copy_topo)
        self.main_layout.addLayout(top_buttons_layout)

        # CREATION OF CONTENT (TABS)
        self.tabs = QtWidgets.QTabWidget()
        has_data = False
        topo_types = [
            ("Bonds", self.system_obj.bonds),
            ("Angles", self.system_obj.angles),
            ("Dihedrals", self.system_obj.dihedrals),
        ]

        for name, df in topo_types:
            if df is not None and not df.empty:
                self.tabs.addTab(self.create_table_tab(name, df), name)
                has_data = True

        if not has_data:
            placeholder = QtWidgets.QLabel("Aucune donnée topologique trouvée.")
            placeholder.setAlignment(QtCore.Qt.AlignCenter)
            self.main_layout.addWidget(placeholder)
        else:
            self.main_layout.addWidget(self.tabs)
            if current_idx < self.tabs.count():
                self.tabs.setCurrentIndex(current_idx)

        if not hasattr(self, "dialog_buttons"):
            self.dialog_buttons = self.create_dialog_buttons(ok_text="Close")
            self.dialog_buttons.button(QtWidgets.QDialogButtonBox.Cancel).hide()

        self.main_layout.addWidget(self.dialog_buttons)

    def clear_layout(self, layout) -> None:
        if layout is not None:
            while layout.count():
                item = layout.takeAt(0)
                if item.widget():
                    item.widget().setParent(None)
                elif item.layout():
                    self.clear_layout(item.layout())

    def create_table_tab(self, name: str, df: pd.DataFrame):
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)

        layout.addWidget(QtWidgets.QLabel(f"<b>{name} types summary</b>"))

        summary = df.copy()
        summary["type"] = summary["type"].astype(str)
        counts = summary.groupby("type").size().reset_index(name="count")

        type_table = QtWidgets.QTableWidget(len(counts), 4)
        type_table.setHorizontalHeaderLabels(["Type", "Count", "FF Parameters", ""])
        type_table.horizontalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.Stretch
        )

        view_only_flags = QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable

        sheet_map = {"Bonds": "bond", "Angles": "angle", "Dihedrals": "dihedral"}
        db_sheet = self.all_sheets.get(sheet_map.get(name), pd.DataFrame())

        for i, row in counts.iterrows():
            t_id = str(row["type"])

            item_type = QtWidgets.QTableWidgetItem(t_id)
            item_type.setFlags(view_only_flags)
            type_table.setItem(i, 0, item_type)

            item_count = QtWidgets.QTableWidgetItem(str(row["count"]))
            item_count.setFlags(view_only_flags)
            type_table.setItem(i, 1, item_count)

            # New button via basic method for styling
            param_label = QtWidgets.QLabel()
            param_label.setAlignment(QtCore.Qt.AlignCenter)

            found_params = self.search_ff_parameters(name, t_id, db_sheet)

            if found_params:
                param_label.setText(f"✓ {found_params['label']}")
                param_label.setToolTip(found_params["details"])
            else:
                param_label.setText("✕ Not defined")

            type_table.setCellWidget(i, 2, param_label)

            # Delete type button
            btn_del_type = QtWidgets.QPushButton("Delete type")
            btn_del_type.clicked.connect(
                lambda ch, n=name, t=t_id: self.delete_by_type(n, t)
            )
            type_table.setCellWidget(i, 3, btn_del_type)

            for col in range(4):
                it = type_table.item(i, col)
                if it:
                    it.setBackground(QtGui.QColor("lightgray"))

        layout.addWidget(type_table)
        return widget

    def search_ff_parameters(
        self, category: str, type_str: str, db_df: pd.DataFrame
    ) -> None:
        if db_df.empty:
            return None
        parts = type_str.split("-")
        translated_parts = []
        for p in parts:
            if p in self.session_ff_dict:
                translated_parts.append(str(self.session_ff_dict[p].get("type", p)))
            else:
                translated_parts.append(p)

        if category == "Bonds" and len(translated_parts) == 2:
            p1, p2 = translated_parts
            mask = ((db_df.iloc[:, 0] == p1) & (db_df.iloc[:, 1] == p2)) | (
                (db_df.iloc[:, 0] == p2) & (db_df.iloc[:, 1] == p1)
            )
        elif category == "Angles" and len(translated_parts) == 3:
            p1, p2, p3 = translated_parts
            mask = (
                (db_df.iloc[:, 0] == p1)
                & (db_df.iloc[:, 1] == p2)
                & (db_df.iloc[:, 2] == p3)
            ) | (
                (db_df.iloc[:, 0] == p3)
                & (db_df.iloc[:, 1] == p2)
                & (db_df.iloc[:, 2] == p1)
            )
        else:
            return None

        match = db_df[mask]
        if not match.empty:
            row = match.iloc[0]
            if category == "Bonds":
                return {
                    "label": "Harmonic",
                    "details": f"K: {row.iloc[2]} | r0: {row.iloc[3]}",
                }
            else:
                return {
                    "label": "Harmonic",
                    "details": f"K: {row.iloc[3]} | theta0: {row.iloc[4]}",
                }
        return None

    def delete_by_type(self, category: str, type_id: str | int) -> None:
        confirm = QtWidgets.QMessageBox.question(
            self, "Confirmation", f"Delete all {category} of type '{type_id}' ?"
        )
        if confirm == QtWidgets.QMessageBox.Yes:
            self.parent().push_undo()
            if category == "Bonds":
                self.system_obj.remove_connection_types(bond_types=[type_id])
            elif category == "Angles":
                self.system_obj.remove_connection_types(angle_types=[type_id])
        main_win = self.parent()
        main_win.sync_ui()

    def open_rule_builder(self):
        if not self.system_obj:
            return
        dialog = RuleBuilderDialog(self.system_obj, parent=self)
        if dialog.exec_():
            new_rule = dialog.get_rule()
            try:
                QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
                main_win = self.parent()
                main_win.push_undo()
                try:
                    self.system_obj.set_topology(new_rule)
                except Exception:
                    main_win.discard_undo()
                    raise
                print(new_rule)
                self.setup_ui()
                main_win = self.parent()
                if hasattr(main_win, "sync_ui"):
                    main_win.sync_ui()
                QtWidgets.QMessageBox.information(self, "Success", "Topology updated!")
            except Exception as e:
                QtWidgets.QMessageBox.critical(
                    self, "Error", f"Failed to apply topology: {e}"
                )
            finally:
                QtWidgets.QApplication.restoreOverrideCursor()

    def highlight_atoms_in_tables(self, selected_indices: list[int]) -> None:
        if not hasattr(self, "tabs") or self.tabs.count() == 0:
            return
        current_widget = self.tabs.currentWidget()
        tables = current_widget.findChildren(QtWidgets.QTableWidget)
        if not tables:
            return
        list_table = tables[-1]

        for row in range(list_table.rowCount()):
            for col in range(list_table.columnCount()):
                item = list_table.item(row, col)
                if item:
                    item.setBackground(QtGui.QColor("white"))

        if not selected_indices:
            return

        highlight_color = QtGui.QColor("#fff9c4")
        first_match = -1
        tab_name = self.tabs.tabText(self.tabs.currentIndex())
        atom_cols = [2, 3]
        if tab_name != "Bonds":
            atom_cols.append(4)

        for row in range(list_table.rowCount()):
            found = False
            for col in atom_cols:
                item = list_table.item(row, col)
                if item:
                    try:
                        if int(item.text()) in selected_indices:
                            found = True
                            break
                    except ValueError:
                        continue

            if found:
                if first_match == -1:
                    first_match = row
                for col in range(list_table.columnCount() - 1):
                    item = list_table.item(row, col)
                    if item:
                        item.setBackground(highlight_color)

        if first_match != -1:
            item_to_scroll = list_table.item(first_match, 0)
            if item_to_scroll:
                list_table.scrollToItem(
                    item_to_scroll, QtWidgets.QAbstractItemView.PositionAtTop
                )

    def copy_topology_coeffs_to_clipboard(self) -> None:
        lines = ["# --- LAMMPS Topology Coefficients ---"]
        count = 0
        if self.system_obj.bonds is not None and not self.system_obj.bonds.empty:
            lines.append("\n# Bond Coefficients")
            db_bond = self.all_sheets.get("bond", pd.DataFrame())
            bond_types = sorted(self.system_obj.bonds["type"].unique())
            for b_type in bond_types:
                params = self.search_ff_parameters("Bonds", str(b_type), db_bond)
                if params:
                    try:
                        k = params["details"].split("|")[0].split(":")[1].strip()
                        r0 = params["details"].split("|")[1].split(":")[1].strip()
                        lines.append(f"bond_coeff {b_type} {k} {r0}  # Harmonic")
                        count += 1
                    except (IndexError, KeyError, TypeError):
                        lines.append(f"# Error parsing bond type {b_type}")

        if self.system_obj.angles is not None and not self.system_obj.angles.empty:
            lines.append("\n# Angle Coefficients")
            db_angle = self.all_sheets.get("angle", pd.DataFrame())
            angle_types = sorted(self.system_obj.angles["type"].unique())
            for a_type in angle_types:
                params = self.search_ff_parameters("Angles", str(a_type), db_angle)
                if params:
                    try:
                        k = params["details"].split("|")[0].split(":")[1].strip()
                        t0 = params["details"].split("|")[1].split(":")[1].strip()
                        lines.append(f"angle_coeff {a_type} {k} {t0}  # Harmonic")
                        count += 1
                    except (IndexError, KeyError, TypeError):
                        lines.append(f"# Error parsing angle type {a_type}")

        full_text = "\n".join(lines)
        clipboard = QtWidgets.QApplication.clipboard()
        clipboard.setText(full_text)
        if count > 0:
            QtWidgets.QMessageBox.information(
                self,
                "Export successful",
                f"Copied {count} topology lines to clipboard.",
            )
        else:
            QtWidgets.QMessageBox.warning(
                self, "Export", "No matching parameters found in database to copy."
            )


class RuleBuilderDialog(QtWidgets.QDialog):
    def __init__(self, atomic_system, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Build Connectivity Rule")
        self.resize(650, 450)
        self.system = atomic_system
        self.neighbor_rules = []
        self.setup_ui()

    def setup_ui(self) -> None:
        self.layout = QtWidgets.QVBoxLayout(self)

        center_group = QtWidgets.QGroupBox("Central atom type")
        c_lay = QtWidgets.QFormLayout(center_group)
        self.combo_center = QtWidgets.QComboBox()
        self.combo_center.addItems([str(t) for t in self.system.atom_types])
        self.new_type_center = QtWidgets.QLineEdit()
        self.new_type_center.setPlaceholderText("New type (optional)")
        c_lay.addRow("Current type :", self.combo_center)
        c_lay.addRow("New type (optional) :", self.new_type_center)
        self.layout.addWidget(center_group)

        self.neighbors_group = QtWidgets.QGroupBox("Neighbors required")
        self.neighbors_layout = QtWidgets.QVBoxLayout(self.neighbors_group)
        self.layout.addWidget(self.neighbors_group)

        btn_add_neighbor = QtWidgets.QPushButton("Add a neighbor rule")
        btn_add_neighbor.clicked.connect(self.add_neighbor_row)
        self.layout.addWidget(btn_add_neighbor)

        opt_group = QtWidgets.QGroupBox("Actions")
        opt_lay = QtWidgets.QHBoxLayout(opt_group)
        self.cb_bonds = QtWidgets.QCheckBox("Create Bonds")
        self.cb_angles = QtWidgets.QCheckBox("Create Angles")
        self.cb_bonds.setChecked(True)
        opt_lay.addWidget(self.cb_bonds)
        opt_lay.addWidget(self.cb_angles)
        self.layout.addWidget(opt_group)

        self.btns = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        )
        self.btns.accepted.connect(self.accept)
        self.btns.rejected.connect(self.reject)
        self.layout.addWidget(self.btns)

    def add_neighbor_row(self) -> None:
        row_widget = QtWidgets.QWidget()
        row_lay = QtWidgets.QHBoxLayout(row_widget)
        row_lay.setContentsMargins(0, 0, 0, 0)

        # Inputs
        combo = QtWidgets.QComboBox()
        combo.addItems([str(t) for t in self.system.atom_types])

        n_spin = QtWidgets.QSpinBox()
        n_spin.setRange(1, 12)
        n_spin.setValue(1)

        dist_spin = QtWidgets.QDoubleSpinBox()
        dist_spin.setRange(0.5, 5.0)
        dist_spin.setValue(1.2)
        dist_spin.setSuffix(" Å")
        dist_spin.setSingleStep(0.1)

        new_t = QtWidgets.QLineEdit()
        new_t.setPlaceholderText("New type")

        # Delete button
        btn_del = QtWidgets.QPushButton("🗑️")
        btn_del.setFixedWidth(30)

        # Storage
        rule_data = {
            "widget": row_widget,
            "type_cb": combo,
            "n_spin": n_spin,
            "dist_spin": dist_spin,
            "new_type": new_t,
        }
        self.neighbor_rules.append(rule_data)
        btn_del.clicked.connect(lambda: self.remove_neighbor_row(rule_data))

        # Assembly
        row_lay.addWidget(QtWidgets.QLabel("type:"))
        row_lay.addWidget(combo)
        row_lay.addWidget(QtWidgets.QLabel("n min:"))
        row_lay.addWidget(n_spin)
        row_lay.addWidget(QtWidgets.QLabel("dist:"))
        row_lay.addWidget(dist_spin)
        row_lay.addWidget(new_t)
        row_lay.addWidget(btn_del)

        self.neighbors_layout.addWidget(row_widget)

    def remove_neighbor_row(self, rule_data: dict) -> None:
        rule_data["widget"].deleteLater()
        self.neighbor_rules.remove(rule_data)

    def get_rule(self) -> TopologyRule:
        """Builds the TopologyRule expected by AtomicSystem.set_topology()."""
        neighbors = [
            NeighborCriterion(
                selection=f"type {r['type_cb'].currentText()}",
                cutoff=r["dist_spin"].value(),
                count=r["n_spin"].value(),
                new_type=r["new_type"].text() if r["new_type"].text() else None,
            )
            for r in self.neighbor_rules
        ]

        return TopologyRule(
            center=f"type {self.combo_center.currentText()}",
            neighbors=neighbors,
            new_type=self.new_type_center.text()
            if self.new_type_center.text()
            else None,
            bonds=self.cb_bonds.isChecked(),
            angles=self.cb_angles.isChecked(),
        )


class NeutralizeDialog(QtWidgets.QDialog):
    def __init__(self, atom_types, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Neutralize System")
        self.setMinimumWidth(300)
        layout = QtWidgets.QVBoxLayout(self)

        layout.addWidget(
            QtWidgets.QLabel("<b>Select types and weights for distribution:</b>")
        )

        self.rows = {}
        grid = QtWidgets.QGridLayout()

        for i, atype in enumerate(atom_types):
            cb = QtWidgets.QCheckBox(str(atype))
            cb.setChecked(True)

            sb = QtWidgets.QDoubleSpinBox()
            sb.setRange(0.01, 10.0)
            sb.setValue(1.0)
            sb.setSingleStep(0.1)
            sb.setEnabled(True)

            cb.toggled.connect(sb.setEnabled)

            grid.addWidget(cb, i, 0)
            grid.addWidget(sb, i, 1)
            self.rows[atype] = (cb, sb)

        layout.addLayout(grid)

        btns = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        )
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

    def get_selected_weights(self) -> None:
        return {
            atype: sb.value() for atype, (cb, sb) in self.rows.items() if cb.isChecked()
        }
