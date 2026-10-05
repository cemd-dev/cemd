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

import numpy as np
from PySide6 import QtCore, QtWidgets


class AtomicModel(QtCore.QAbstractTableModel):
    def __init__(self, system_obj):
        super().__init__()
        self.system_obj = system_obj
        self._last_known_len = len(system_obj.atoms)
        self.columns_order = ["type", "x", "y", "z", "charge"]

    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self.system_obj.atoms)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return len(self.columns_order) + 1

    def flags(self, index):
        if not index.isValid():
            return QtCore.Qt.ItemFlag.NoItemFlags
        fl = QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable
        if index.column() != 0:
            fl |= QtCore.Qt.ItemFlag.ItemIsEditable
        return fl

    def data(self, index, role=QtCore.Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None

        row = index.row()
        col = index.column()

        if role in [
            QtCore.Qt.ItemDataRole.DisplayRole,
            QtCore.Qt.ItemDataRole.EditRole,
        ]:
            if col == 0:
                return str(self.system_obj.atoms.index[row])
            col_name = self.columns_order[col - 1]
            val = self.system_obj.atoms.iloc[row].get(col_name, "")
            # Clean formatting for display
            return f"{val:.4f}" if isinstance(val, float) else str(val)

        if role == QtCore.Qt.ItemDataRole.TextAlignmentRole:
            return QtCore.Qt.AlignmentFlag.AlignCenter
        return None

    def setData(self, index, value, role=QtCore.Qt.ItemDataRole.EditRole):
        """C'est ici que la magie opère : on écrit dans le DataFrame original."""
        if role == QtCore.Qt.ItemDataRole.EditRole:
            col = index.column()
            if col == 0:
                return False

            row = index.row()
            col_name = self.columns_order[col - 1]

            try:
                if col_name in ["x", "y", "z", "charge"]:
                    new_value = float(value)
                    if not np.isfinite(new_value):
                        return False
                else:
                    new_value = str(value)

                self.system_obj.atoms.at[self.system_obj.atoms.index[row], col_name] = (
                    new_value
                )

                self.dataChanged.emit(index, index)

                return True

            except (ValueError, TypeError):
                return False
        return False

    def headerData(self, section, orientation, role):
        if (
            role == QtCore.Qt.ItemDataRole.DisplayRole
            and orientation == QtCore.Qt.Orientation.Horizontal
        ):
            if section == 0:
                return "ID"
            return self.columns_order[section - 1].capitalize()
        return None


class AtomTable(QtWidgets.QTableView):  # We change QTableWidget to QTableView
    def __init__(self):
        super().__init__()
        self.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection
        )
        self.horizontalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeMode.Stretch
        )
        self.verticalHeader().setVisible(True)  # Optional: to see the IDs on the left

    def update_data(self, system_obj) -> None:
        current_model = self.model()
        if current_model is not None:
            same_obj = getattr(current_model, "system_obj", None) is system_obj
            last_len = getattr(current_model, "_last_known_len", 0)
            current_len = len(system_obj.atoms)

            if same_obj and last_len == current_len:
                return

        # Clean up the old model
        if current_model is not None:
            try:
                current_model.dataChanged.disconnect()
            except (RuntimeError, TypeError):  # nothing connected yet
                pass

        # Create the model (Force the reconstruction)
        model = AtomicModel(system_obj)
        self.setModel(model)
