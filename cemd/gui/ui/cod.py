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

import json
import os
import webbrowser
from typing import TYPE_CHECKING, Any

from PySide6 import QtCore, QtWidgets

if TYPE_CHECKING:
    from ..main_window import AtomViewerGUI

from ...core._io.sources.cod import (
    cod_search_by_elements,
    cod_search_by_id,
    cod_search_by_name,
    formula_to_elements,
    get_structure_by_cod_id,
)
from .. import _userdata
from .base_dialog import BaseBuilderDialog


class CODBrowserDialog(BaseBuilderDialog):
    def __init__(self, parent: AtomViewerGUI = None):
        # Use the BaseBuilderDialog constructor (parent, title, width)
        super().__init__(parent, "COD Database Explorer", 850)
        self.setMinimumHeight(600)

        self.cache_file = _userdata.cache_file("cod_cache.json")

        self.selected_system = None
        self.last_cod_id = None
        self.results_data = []

        self.setup_ui()

        QtCore.QTimer.singleShot(100, self.load_last_search)

    def setup_ui(self) -> None:
        # Use self.main_layout provided by BaseBuilderDialog
        layout = self.main_layout

        # Search section
        search_group = QtWidgets.QGroupBox("Search in COD")
        search_layout = QtWidgets.QVBoxLayout(search_group)

        top_row = QtWidgets.QHBoxLayout()
        self.combo_type = QtWidgets.QComboBox()
        self.combo_type.addItems(["Mineral Name", "Chemical Formula", "COD ID"])
        self.combo_type.setFixedWidth(150)

        self.search_input = QtWidgets.QLineEdit()
        self.search_input.setPlaceholderText("Enter name, formula or ID...")
        self.search_input.returnPressed.connect(self.run_search)

        self.btn_search = self.create_icon_button("Search", "search", primary=True)
        self.btn_search.clicked.connect(self.run_search)

        top_row.addWidget(self.combo_type)
        top_row.addWidget(self.search_input)
        top_row.addWidget(self.btn_search)
        search_layout.addLayout(top_row)
        layout.addWidget(search_group)

        # Result table
        self.table = self.create_table(
            ["COD ID", "Common Name", "Formula", "Space Group", "Cell", "Ref"],
            selectable=True,
        )

        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.Fixed)
        self.table.setColumnWidth(0, 100)
        header.setSectionResizeMode(3, QtWidgets.QHeaderView.Fixed)
        self.table.setColumnWidth(3, 100)
        header.setSectionResizeMode(5, QtWidgets.QHeaderView.Fixed)
        self.table.setColumnWidth(5, 50)

        self.table.verticalHeader().setDefaultSectionSize(32)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.doubleClicked.connect(self.load_selected)

        self.table.itemSelectionChanged.connect(
            lambda: self.btn_import.setEnabled(True)
        )
        self.table.setSortingEnabled(True)

        layout.addWidget(self.table)

        # Action section
        bottom_layout = QtWidgets.QHBoxLayout()

        self.btn_import = self.create_action_button(
            "Open selected structure", primary=True
        )
        self.btn_import.setMinimumHeight(35)
        self.btn_import.setEnabled(False)
        self.btn_import.clicked.connect(self.load_selected)

        bottom_layout.addStretch()
        bottom_layout.addWidget(self.btn_import)

        layout.addLayout(bottom_layout)

        self.status_bar = QtWidgets.QStatusBar()
        self.status_bar.setSizeGripEnabled(False)
        self.status_bar.showMessage("Ready")
        layout.addWidget(self.status_bar)

    def run_search(self) -> None:
        query = self.search_input.text().strip()
        if not query:
            return

        search_mode = self.combo_type.currentText()
        self.table.setRowCount(0)
        self.status_bar.showMessage(f"⏱ Searching in COD for '{query}'...")
        QtCore.QCoreApplication.processEvents()  # to show the message

        try:
            results = []
            if search_mode == "Mineral Name":
                results = cod_search_by_name(query)
            elif search_mode == "Chemical Formula":
                elements = formula_to_elements(query)
                results = cod_search_by_elements(elements)
            elif search_mode == "COD ID":
                results = cod_search_by_id(query)

            self.display_results(results)
            self.save_last_search(query, search_mode, results)
            self.status_bar.showMessage(f"{len(results)} structures found")

        except Exception as e:
            self.status_bar.showMessage("⚠ Error during search")
            QtWidgets.QMessageBox.critical(self, "Search Error", str(e))

    def save_last_search(
        self, query: str, mode: str, results: list[dict[str, Any]]
    ) -> None:
        """Saves the search to a local JSON file."""
        data = {"query": query, "mode": mode, "results": results}
        try:
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=4)
        except Exception as e:
            print(f"Error saving JSON cache: {e}")

    def load_last_search(self) -> None:
        """Reads the JSON file and populates the interface."""
        if not os.path.exists(self.cache_file):
            return

        try:
            with open(self.cache_file, encoding="utf-8") as f:
                data = json.load(f)

            query = data.get("query", "")
            mode = data.get("mode", "Mineral Name")
            results = data.get("results", [])

            if query:
                self.search_input.setText(query)
                self.combo_type.setCurrentText(mode)
                self.results_data = results
                self.display_results(results)
                self.status_bar.showMessage(
                    f"Last search restored ({len(results)} results)"
                )
        except Exception as e:
            print(f"Error loading JSON cache: {e}")

    def display_results(self, results: dict[str, Any]) -> None:
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)

        for r in results:
            row = self.table.rowCount()
            self.table.insertRow(row)

            id_item = QtWidgets.QTableWidgetItem(str(r.get("id", "")))
            id_item.setTextAlignment(QtCore.Qt.AlignCenter)
            self.table.setItem(row, 0, id_item)

            self.table.setItem(
                row, 1, QtWidgets.QTableWidgetItem(r.get("common_name", "N/A"))
            )
            self.table.setItem(
                row, 2, QtWidgets.QTableWidgetItem(r.get("formula", "N/A"))
            )
            self.table.setItem(
                row, 3, QtWidgets.QTableWidgetItem(r.get("spacegroup", "N/A"))
            )
            self.table.setItem(row, 4, QtWidgets.QTableWidgetItem(r.get("cell", "N/A")))

            doi = r.get("doi")
            if doi and doi not in ["N/A", "None"]:
                container = QtWidgets.QWidget()
                btn_layout = QtWidgets.QHBoxLayout(container)
                btn_layout.setContentsMargins(0, 0, 0, 0)
                btn_layout.setAlignment(QtCore.Qt.AlignCenter)

                btn_ref = self.create_icon_button("", "reference")
                btn_ref.setIconSize(QtCore.QSize(18, 18))
                btn_ref.setFixedSize(24, 24)
                btn_ref.setFlat(True)
                btn_ref.setCursor(QtCore.Qt.PointingHandCursor)
                btn_ref.setToolTip(f"Open DOI: {doi}")

                doi_url = f"https://doi.org/{doi}"
                btn_ref.clicked.connect(
                    lambda checked=False, url=doi_url: webbrowser.open(url)
                )

                btn_layout.addWidget(btn_ref)
                self.table.setCellWidget(row, 5, container)
            else:
                item_dash = QtWidgets.QTableWidgetItem("-")
                item_dash.setTextAlignment(QtCore.Qt.AlignCenter)
                self.table.setItem(row, 5, item_dash)

        self.table.setSortingEnabled(True)

    def load_selected(self) -> None:
        row = self.table.currentRow()
        if row < 0:
            return
        self.last_cod_id = self.table.item(row, 0).text()
        try:
            system = get_structure_by_cod_id(self.last_cod_id)
            if system:
                self.selected_system = system  # ← already an AtomicSystem
                self.accept()
            else:
                raise ValueError("Could not parse CIF.")
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Import Error", str(e))
