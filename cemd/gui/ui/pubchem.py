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
from typing import TYPE_CHECKING

from PySide6 import QtCore, QtWidgets

from ...core._io.sources.pubchem import (
    get_pubchem_details,
    get_structure,
    pubchem_sdq_search,
    pubchem_search_by_name,
)
from .. import _userdata
from .base_dialog import BaseBuilderDialog

if TYPE_CHECKING:
    from ..main_window import AtomViewerGUI


class PubChemBrowserDialog(BaseBuilderDialog):
    def __init__(self, parent: AtomViewerGUI = None) -> None:
        super().__init__(parent, "PubChem Database Explorer", 850)
        self.setMinimumHeight(600)

        self.cache_file = _userdata.cache_file("pubchem_cache.json")

        self.selected_system = None
        self.setup_ui()
        QtCore.QTimer.singleShot(100, self.load_last_search)

    def setup_ui(self) -> None:
        layout = self.main_layout

        # Filters
        search_group = QtWidgets.QGroupBox("Search in PubChem")
        search_layout = QtWidgets.QVBoxLayout(search_group)

        top_row = QtWidgets.QHBoxLayout()
        self.combo_type = QtWidgets.QComboBox()
        self.combo_type.addItems(["Compound Name or Formula", "PubChem CID"])
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

        # Table
        self.table = self.create_table(
            ["CID", "Name", "Formula", "Mol. Weight", "SMILES", "Ref."], selectable=True
        )

        self.table.setColumnWidth(0, 80)
        self.table.setColumnWidth(5, 50)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.doubleClicked.connect(self.load_selected)
        layout.addWidget(self.table)

        # Actions
        bottom_layout = QtWidgets.QHBoxLayout()
        self.btn_import = self.create_action_button(
            "Open selected structure", primary=True
        )
        self.btn_import.clicked.connect(self.load_selected)
        self.btn_import.setEnabled(False)

        bottom_layout.addStretch()
        bottom_layout.addWidget(self.btn_import)
        layout.addLayout(bottom_layout)

        self.table.itemSelectionChanged.connect(
            lambda: self.btn_import.setEnabled(True)
        )

        self.status_bar = QtWidgets.QStatusBar()
        self.status_bar.setSizeGripEnabled(False)
        self.status_bar.showMessage("Ready")
        layout.addWidget(self.status_bar)

    def run_search(self) -> None:
        query = self.search_input.text().strip()
        if not query:
            return

        mode = self.combo_type.currentText()
        self.status_bar.showMessage(f"⏱ Searching in PubChem for '{query}'...")
        # Force the refresh of the UI to display the message
        QtCore.QCoreApplication.processEvents()

        try:
            results = []
            if mode == "Compound Name or Formula":
                results = pubchem_search_by_name(query)
                if len(results) == 0:
                    results = pubchem_sdq_search(query)

            elif mode == "PubChem CID":
                # Direct search by ID (we keep your details logic)
                details = get_pubchem_details(query)
                if details:
                    results = [
                        {
                            "id": query,
                            "common_name": details.get("title", "N/A"),
                            "formula": details.get("formula", "N/A"),
                            "weight": details.get("weight", "N/A"),
                            "smiles": details.get("smiles", "N/A"),
                        }
                    ]

            # Display and feedback
            if results:
                self.display_results(results)
                self.save_last_search(query, mode, results)
                self.status_bar.showMessage(f"{len(results)} compounds found")
            else:
                self.table.setRowCount(0)
                self.status_bar.showMessage(f"No results found for '{query}'")

        except Exception as e:
            self.status_bar.showMessage("⚠ Error during search")
            QtWidgets.QMessageBox.critical(
                self, "Search Error", f"Search failed: {str(e)}"
            )

    def display_results(self, results: list) -> None:
        self.table.setRowCount(0)
        for r in results:
            row = self.table.rowCount()
            self.table.insertRow(row)

            self.table.setItem(row, 0, QtWidgets.QTableWidgetItem(str(r["id"])))
            self.table.setItem(row, 1, QtWidgets.QTableWidgetItem(r["common_name"]))
            self.table.setItem(row, 2, QtWidgets.QTableWidgetItem(r["formula"]))
            self.table.setItem(row, 3, QtWidgets.QTableWidgetItem(str(r["weight"])))
            self.table.setItem(row, 4, QtWidgets.QTableWidgetItem(r["smiles"]))

            # Web link button
            btn_web = self.create_icon_button("", "reference")  # Reuses COD icon
            url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{r['id']}"
            btn_web.clicked.connect(lambda chk=False, u=url: webbrowser.open(u))
            self.table.setCellWidget(row, 5, btn_web)

    def load_selected(self) -> None:
        row = self.table.currentRow()
        if row < 0:
            return

        cid = self.table.item(row, 0).text()
        smiles = self.table.item(row, 4).text()
        self.status_bar.showMessage(f"Fetching 3D coordinates for CID {cid}...")

        try:
            system = get_structure(cid, smiles)  # ← already an AtomicSystem
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                self, "3D Error", f"Failed to fetch structure for CID {cid}:\n{e}"
            )
            return

        if system:
            self.selected_system = system
            self.accept()
        else:
            QtWidgets.QMessageBox.warning(
                self,
                "3D Error",
                "No 3D conformer available for this compound on PubChem.",
            )

    def save_last_search(self, query, mode, results) -> None:
        try:
            with open(self.cache_file, "w") as f:
                json.dump({"query": query, "mode": mode, "results": results}, f)
        except:
            pass

    def load_last_search(self) -> None:
        if not os.path.exists(self.cache_file):
            return
        try:
            with open(self.cache_file) as f:
                data = json.load(f)
                self.search_input.setText(data["query"])
                self.combo_type.setCurrentText(data["mode"])
                self.display_results(data["results"])
        except:
            pass
