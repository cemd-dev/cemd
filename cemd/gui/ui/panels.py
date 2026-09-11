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

import json
import os

from PySide6 import QtCore, QtGui, QtWidgets
from PySide6Qlementine import ColorButton, ColorMode

from .gui_utils import create_icon_button


class SystemSummaryPanel(QtWidgets.QFrame):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(2)

        # Physical Properties
        self.lbl_stats = QtWidgets.QLabel()
        layout.addWidget(self.lbl_stats)

        # Separator
        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        layout.addWidget(line)

        # Box Parameters (HTML table)
        self.lbl_box = QtWidgets.QLabel()
        layout.addWidget(self.lbl_box)

    def update_info(self, system):
        """Updates the frame from an AtomicSystem object"""
        if not system:
            self.lbl_stats.setText("No system loaded")
            self.lbl_box.clear()
            return

        charge = system.total_charge
        is_neutral = abs(charge) < 1e-3
        charge_style = (
            "color: #d32f2f; font-weight: bold;"
            if not is_neutral
            else "color: #263238;"
        )

        # --- Stats part (Uses your @property) ---
        stats_html = f"""
        <div style="line-height: 140%;">
            <span style="font-size: 10pt;">
                <b>Total charge:</b> 
                <span style="{charge_style}">{charge:.4f} e</span> 
            </span><br>
            <b>Volume:</b> {system.volume / 1e3:.2f} nm³<br>
            <b>Density:</b> {system.density:.2f} g/cm³<br>
            <span style="color: #607d8b; font-size: 9pt; display: block; margin-top: 5px;">
                &lt;AtomicSystem with {system.num_atoms} atoms, {system.num_bonds} bonds, {system.num_angles} angles&gt;
            </span>
        </div>
        """
        self.lbl_stats.setText(stats_html)

        # ---Box part (Uses system.box) ---
        b = system.box  # [a, b, c, alpha, beta, gamma]
        box_html = f"""
        <b style="color: #455a64; font-size: 9pt;">Box Parameters</b>
        <table width="100%" style="margin-top: 2px; border-collapse: collapse; line-height: 100%;">
            <tr style="color: #90a4ae; font-size: 8pt;">
                <th align="left" style="padding: 0px;">a (Å)</th>
                <th align="left" style="padding: 0px;">b (Å)</th>
                <th align="left" style="padding: 0px;">c (Å)</th>
            </tr>
            <tr style="font-size: 9pt;">
                <td style="padding: 0px 0px 2px 0px;">{b[0]:.2f}</td>
                <td style="padding: 0px 0px 2px 0px;">{b[1]:.2f}</td>
                <td style="padding: 0px 0px 2px 0px;">{b[2]:.2f}</td>
            </tr>
            <tr style="color: #90a4ae; font-size: 8pt;">
                <th align="left" style="padding: 0px;">α (°)</th>
                <th align="left" style="padding: 0px;">β (°)</th>
                <th align="left" style="padding: 0px;">γ (°)</th>
            </tr>
            <tr style="font-size: 9pt;">
                <td style="padding: 0px;">{b[3]:.2f}</td>
                <td style="padding: 0px;">{b[4]:.2f}</td>
                <td style="padding: 0px;">{b[5]:.2f}</td>
            </tr>
        </table>
        """
        self.lbl_box.setText(box_html)


class BaseManagerPanel(QtWidgets.QGroupBox):
    """Parent class that manages the look and structure Scroll + Button"""

    def __init__(self, title, height):
        super().__init__(title)
        self.setFixedHeight(height)

        # Main layout
        self.main_layout = QtWidgets.QVBoxLayout(self)
        self.main_layout.setContentsMargins(8, 15, 8, 8)

        # White scroll area (the central rectangle)
        self.scroll = QtWidgets.QScrollArea()
        self.scroll.setWidgetResizable(True)

        self.scroll_content = QtWidgets.QWidget()
        self.list_layout = QtWidgets.QVBoxLayout(self.scroll_content)
        self.list_layout.setAlignment(QtCore.Qt.AlignTop)
        self.list_layout.setContentsMargins(5, 5, 5, 5)
        self.list_layout.setSpacing(2)

        self.scroll.setWidget(self.scroll_content)
        self.main_layout.addWidget(self.scroll)

    def add_action_button(self, text):
        """Adds the standard button at the bottom, outside the scroll"""
        btn = QtWidgets.QPushButton(text)
        self.main_layout.addWidget(btn)
        return btn


class FilterPanel(BaseManagerPanel):
    type_changed = QtCore.Signal(dict)

    def __init__(self):
        # We initialize the base with the title and the height
        super().__init__("Types View Manager", 400)

        self.checkboxes = {}
        self.radius_spinboxes = {}
        self.color_buttons = {}

        # ---GLOBAL SCALING CONTROL (Slider) ---
        self.scale_container = QtWidgets.QWidget()
        scale_layout = QtWidgets.QHBoxLayout(self.scale_container)
        scale_layout.setContentsMargins(5, 5, 5, 5)
        scale_layout.setSpacing(10)

        # Title label
        lbl_title = QtWidgets.QLabel("Atom Scale:")

        # The Slider
        # Qt sliders use integers, so we will map 10 -> 100 to have precision
        self.global_scale_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.global_scale_slider.setRange(10, 300)  # From 0.1x to 3.0x
        self.global_scale_slider.setValue(100)  # 1.0x by default
        self.global_scale_slider.valueChanged.connect(self._on_slider_move_only)
        self.global_scale_slider.sliderReleased.connect(self._emit_everything)

        # Value label (to see the number)
        self.lbl_scale_val = QtWidgets.QLabel("x1.0")
        self.lbl_scale_val.setFixedWidth(35)

        scale_layout.addWidget(lbl_title)
        scale_layout.addWidget(self.global_scale_slider)
        scale_layout.addWidget(self.lbl_scale_val)

        # Insertion above the main Reset button (the last element of the layout is the Reset button)
        # We insert at index count -1
        self.layout().insertWidget(self.layout().count() - 1, self.scale_container)

        # Action button inherited from BaseManagerPanel
        self.btn_reset_vdw = self.add_action_button("Reset")

    def _on_slider_move_only(self, value):
        """Updates the text AND starts rendering in real time"""
        scale_float = value / 100.0
        self.lbl_scale_val.setText(f"x{scale_float:.1f}")

    def get_scale_value(self):
        """Returns the float value of the scale"""
        return self.global_scale_slider.value() / 100.0

    def set_scale_value(self, value):
        """Updates the slider and label from a float value (ex: 1.0)."""
        if hasattr(self, "global_scale_slider"):
            self.global_scale_slider.blockSignals(True)
            self.global_scale_slider.setValue(int(value * 100))
            self.global_scale_slider.blockSignals(False)

        if hasattr(self, "lbl_scale_val"):
            self.lbl_scale_val.setText(f"x{value:.1f}")

    def refresh(self, data, color_map, radius_map, global_scale):
        """
        Reconstructs the content from the AtomicSystem (data) object.
        Adds atom count by type.
        """
        # Clean up the legacy layout
        count_removed = 0
        while self.list_layout.count():
            item = self.list_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        self.checkboxes = {}
        self.radius_spinboxes = {}
        self.color_buttons = {}
        main_win = self.window()

        # ---NEW: Calculation of the number of atoms by type ---
        # We count the occurrences in the atoms DataFrame
        counts = data.atoms["type"].value_counts()
        counts.index = counts.index.astype(str)

        # We retrieve the properties of the data object
        types = data.atom_types

        # Adding lines for each type
        for i, atype in enumerate(types):
            atype_str = str(atype)

            row_widget = QtWidgets.QWidget()
            row_layout = QtWidgets.QHBoxLayout(row_widget)
            row_layout.setContentsMargins(5, 2, 5, 2)
            row_layout.setSpacing(5)  # Reduced spacing between elements on the left

            # Checkbox (Tight left) ---
            cb = QtWidgets.QCheckBox(atype_str)
            cb.setChecked(True)
            cb.stateChanged.connect(self._emit_everything)
            # Prevent the checkbox from taking up the full width
            cb.setSizePolicy(
                QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred
            )
            self.checkboxes[atype_str] = cb

            # The number (Paste right after) ---
            n_atoms = counts.get(atype, 0)
            count_label = QtWidgets.QLabel(f"[{n_atoms}]")
            count_label.setSizePolicy(
                QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred
            )

            # Tooltip on both for added comfort
            tip = f"Total pour {atype_str}: {n_atoms} atomes"
            cb.setToolTip(tip)
            count_label.setToolTip(tip)

            # Color (Default Gray [2026-03-08]) ---
            color_btn = ColorButton(row_widget)
            color_btn.setFixedSize(16, 16)
            color_btn.setColorMode(ColorMode.RGB)
            current_color = color_map.get(atype_str, "#808080")
            color_btn.setColor(QtGui.QColor(current_color))
            color_btn.clicked.connect(
                lambda _, t=atype_str, b=color_btn: self.pick_color(t, b, main_win)
            )
            self.color_buttons[atype_str] = color_btn

            # Spinbox (VdW Radius) ---
            sp = QtWidgets.QDoubleSpinBox()
            sp.setRange(0.1, 5.0)
            sp.setSingleStep(0.1)
            sp.setSuffix(" Å")
            sp.setValue(radius_map.get(atype_str, 1.0))
            sp.setFixedWidth(70)
            sp.valueChanged.connect(self._emit_everything)
            self.radius_spinboxes[atype_str] = sp

            row_layout.addWidget(cb)
            row_layout.addWidget(count_label)
            row_layout.addStretch()
            row_layout.addWidget(color_btn)
            row_layout.addWidget(sp)

            self.list_layout.addWidget(row_widget)

        # Reconnecting the Reset button
        try:
            self.btn_reset_vdw.clicked.disconnect()
        except TypeError:
            pass
        self.btn_reset_vdw.clicked.connect(self.reset_to_config)

        self.global_scale_slider.blockSignals(True)
        # We multiply by 100 because the slider uses integers (1.0 -> 100)
        self.global_scale_slider.setValue(int(global_scale * 100))
        self.lbl_scale_val.setText(f"x{global_scale:.1f}")
        self.global_scale_slider.blockSignals(False)

    def pick_color(self, atype_str, btn, main_win):
        """Opens the color picker initialized to the current color"""

        current_config = main_win.global_config
        color_map = current_config.get("color_map", {})

        current_hex = color_map.get(atype_str, "#ffffff")

        color = QtWidgets.QColorDialog.getColor(
            QtGui.QColor(current_hex), main_win, f"Select color for {atype_str}"
        )

        if color.isValid():
            new_hex = color.name()

            color_map[atype_str] = new_hex

            active_tab = main_win.tabs.currentWidget()
            if active_tab:
                active_tab.refresh_tab_view(full_rebuild=False)

    def get_settings(self):
        """Returns the settings: real radii, visibility AND scale factor"""
        radii = {}
        visibility = {}
        colors_instat = {t: btn.color().name() for t, btn in self.color_buttons.items()}

        # Recover the raw values ​​of the spinboxes
        for atype, sp in self.radius_spinboxes.items():
            radii[atype] = sp.value()

        for atype, cb in self.checkboxes.items():
            visibility[atype] = cb.isChecked()

        # Add the scale factor to the dictionary
        return {
            "radii": radii,
            "visibility": visibility,
            "colors": colors_instat,
            "global_scale": self.global_scale_slider.value(),
        }

    def reset_to_config(self):
        """Resets radii and scale from default_config.json"""

        base_dir = os.path.dirname(os.path.realpath(__file__))
        default_path = os.path.join(base_dir, "default_config.json")

        default_radii = {}
        default_scale = 1.0

        if os.path.exists(default_path):
            try:
                with open(default_path, encoding="utf-8") as f:
                    data = json.load(f)
                    default_radii = data.get("radius_map", {})
                    default_scale = data.get("global_scale", 1.0)
            except Exception as e:
                print(f"Erreur lors du reset (lecture JSON): {e}")

        # Block the signals to avoid refreshing the plotter at each iteration
        self.blockSignals(True)

        self.global_scale_slider.setValue(int(default_scale * 100))
        if hasattr(self, "lbl_scale_val"):
            self.lbl_scale_val.setText(f"x{default_scale:.1f}")

        for atype, sp in self.radius_spinboxes.items():
            val = default_radii.get(atype, 1.0)
            sp.setValue(val)

        self.blockSignals(False)
        self._emit_everything()

    def _emit_everything(self, *args):
        """Sends the complete dictionary to the Plotter"""
        settings = self.get_settings()
        self.type_changed.emit(settings)


class BondManagerPanel(BaseManagerPanel):
    bond_settings_changed = QtCore.Signal(dict)

    def __init__(self):
        # Base initialization (Title, Height)
        super().__init__("Bonds View Manager", 400)

        self.available_types = []
        self.pair_rules = []

        # ---THICKNESS CONTROL (Slider) ---
        self.scale_container = QtWidgets.QWidget()
        scale_layout = QtWidgets.QHBoxLayout(self.scale_container)
        scale_layout.setContentsMargins(5, 5, 5, 5)

        lbl_title = QtWidgets.QLabel("Bond Radius:")

        self.bond_radius_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.bond_radius_slider.setRange(1, 50)  # From 0.01 to 0.5
        self.bond_radius_slider.setValue(10)  # 0.1 by default
        self.bond_radius_slider.valueChanged.connect(self._on_slider_move)
        # We emit the final signal when we release to avoid rowing
        self.bond_radius_slider.sliderReleased.connect(self.emit_changes)

        self.lbl_radius_val = QtWidgets.QLabel("0.10x")
        self.lbl_radius_val.setFixedWidth(35)

        scale_layout.addWidget(lbl_title)
        scale_layout.addWidget(self.bond_radius_slider)
        scale_layout.addWidget(self.lbl_radius_val)

        # Insertion above the “Add Bond” button
        self.main_layout.insertWidget(0, self.scale_container)

        # Action button via BaseManagerPanel
        self.btn_add_pair = self.add_action_button("Add bond")
        self.btn_add_pair.clicked.connect(lambda: self.add_pair_row())

    def _on_slider_move(self, value):
        """Updates the text next to the slider without launching heavy rendering"""
        scale_float = value / 100.0
        self.lbl_radius_val.setText(f"x{scale_float:.2f}")

    def get_bond_radius(self):
        """Returns the real value of the radius (ex: 0.15)"""
        return self.bond_radius_slider.value() / 100.0

    def refresh(self, system, bond_map):
        if system is not None and hasattr(system, "atom_types"):
            self.available_types = [str(t) for t in system.atom_types]
        else:
            self.available_types = []

        self._is_refreshing = True
        try:
            self.set_settings(bond_map)
        finally:
            self._is_refreshing = False

    def add_pair_row(self, t1=None, t2=None, dist_val=1.5):
        """Adds a linking rule line"""
        if not self.available_types:
            return
        if not self.available_types:
            return

        row_widget = QtWidgets.QWidget()
        row_widget.setFixedHeight(35)
        row_layout = QtWidgets.QHBoxLayout(row_widget)
        row_layout.setContentsMargins(5, 2, 5, 2)
        row_layout.setSpacing(6)

        # Combo 1 & 2
        c1 = QtWidgets.QComboBox(row_widget)
        c2 = QtWidgets.QComboBox(row_widget)

        for c in [c1, c2]:
            c.blockSignals(True)
            c.addItems(self.available_types)
            c.blockSignals(False)

        if t1 in self.available_types:
            c1.setCurrentText(t1)
        if t2 in self.available_types:
            c2.setCurrentText(t2)

        # Distance SpinBox
        dist = QtWidgets.QDoubleSpinBox()
        dist.setRange(0.1, 10.0)
        dist.setSingleStep(0.1)
        dist.blockSignals(True)
        dist.setValue(dist_val)
        dist.blockSignals(False)
        dist.setSuffix(" Å")
        dist.setFixedWidth(85)

        # Delete button
        btn_del = create_icon_button("", "trash")
        btn_del.setFixedWidth(30)

        # Layout
        for w in [c1, c2, dist, btn_del]:
            row_layout.addWidget(w)

        self.list_layout.addWidget(row_widget)
        row_widget.show()

        rule = {"c1": c1, "c2": c2, "dist": dist, "widget": row_widget}
        self.pair_rules.append(rule)

        # Connections
        c1.currentIndexChanged.connect(self.emit_changes)
        c2.currentIndexChanged.connect(self.emit_changes)
        dist.valueChanged.connect(self.emit_changes)
        btn_del.clicked.connect(lambda: self.remove_pair_row(rule))

    def remove_pair_row(self, rule):
        rule["widget"].deleteLater()
        if rule in self.pair_rules:
            self.pair_rules.remove(rule)
        self.emit_changes()

    def get_settings(self):
        """Reads the UI and returns the current bond_map"""
        settings = {}
        for r in self.pair_rules:
            try:
                t1 = r["c1"].currentText()
                t2 = r["c2"].currentText()
                if t1 and t2:
                    key = "-".join(sorted([t1, t2]))
                    settings[key] = r["dist"].value()
            except RuntimeError:
                continue
        settings["global_bond_radius"] = self.get_bond_radius()
        return settings

    def set_settings(self, bond_map):

        while self.list_layout.count():
            item = self.list_layout.takeAt(0)
            widget = item.widget()
            if widget:
                # widget.setParent(None)
                widget.deleteLater()
        self.pair_rules = []

        current_types = set(self.available_types)
        for pair_str, dist_value in bond_map.items():
            elements = pair_str.split("-")
            if len(elements) == 2:
                t1, t2 = elements[0], elements[1]
                if t1 in current_types and t2 in current_types:
                    self.add_pair_row(t1, t2, dist_value)

    def emit_changes(self, *args):
        # We do not transmit if the panel is emptying (signals blocked)
        if not getattr(self, "_is_refreshing", False):
            self.bond_settings_changed.emit(self.get_settings())
