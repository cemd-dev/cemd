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

"""How structures are drawn: one set of settings shared by every tab.

The colors, radii, bond rules and background used to live in four places
at once -- each tab's plotter, a cache on the main window, the side
panels' widgets and a dictionary rebuilt on demand -- kept in step by
copies running in every direction. They are one `DisplaySettings` now:
the main window owns it, every plotter reads that same object, the panels
write to it when the user changes something, and it is what gets saved.
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from cemd.gui import _userdata

#: Settings shipped with the package; the user's file only stores what
#: differs from them.
DEFAULT_CONFIG = Path(__file__).resolve().parent / "default_config.json"

#: Backgrounds the "Background color" action cycles through.
BG_CYCLE = ("black", "#1A1A1A", "white", "#B0C4DE")

#: Backgrounds light enough to need dark text and box edges.
LIGHT_BACKGROUNDS = ("white", "#B0C4DE")


@dataclass
class DisplaySettings:
    color_map: dict[str, str] = field(default_factory=dict)
    radius_map: dict[str, float] = field(default_factory=dict)
    #: Bond cutoffs in Å, keyed "A-B" with the two types sorted.
    bond_map: dict[str, float] = field(default_factory=dict)
    global_scale: float = 1.0
    global_bond_radius: float = 0.1
    bg_color: str = BG_CYCLE[0]

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> DisplaySettings:
        bg_color = data.get("bg_color", BG_CYCLE[0])
        return cls(
            color_map=dict(data.get("color_map", {})),
            radius_map=dict(data.get("radius_map", {})),
            bond_map=dict(data.get("bond_map", {})),
            global_scale=float(data.get("global_scale", 1.0)),
            global_bond_radius=float(data.get("global_bond_radius", 0.1)),
            bg_color=bg_color if bg_color in BG_CYCLE else BG_CYCLE[0],
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "color_map": dict(self.color_map),
            "radius_map": dict(self.radius_map),
            "bond_map": dict(self.bond_map),
            "global_scale": self.global_scale,
            "global_bond_radius": self.global_bond_radius,
            "bg_color": self.bg_color,
        }

    # ------------------------------------------------------------------
    # Background
    # ------------------------------------------------------------------

    @property
    def is_light_background(self) -> bool:
        return self.bg_color.lower() in (c.lower() for c in LIGHT_BACKGROUNDS)

    def next_background(self) -> str:
        """Move to the next background of the cycle and return it."""
        index = BG_CYCLE.index(self.bg_color) if self.bg_color in BG_CYCLE else -1
        self.bg_color = BG_CYCLE[(index + 1) % len(BG_CYCLE)]
        return self.bg_color

    # ------------------------------------------------------------------
    # Bond rules
    # ------------------------------------------------------------------

    def bond_rules_for(self, types: Iterable[str]) -> dict[str, float]:
        """The bond rules between two of `types`, the ones a system uses.

        Rules for other structures stay in `bond_map` but are left out:
        the longest cutoff sets the neighbour-search radius.
        """
        present = {str(t) for t in types}
        return {
            key: cutoff
            for key, cutoff in self.bond_map.items()
            if _pair(key) is not None and _pair(key) <= present
        }

    def set_bond_rules(self, types: Iterable[str], rules: dict[str, float]) -> None:
        """Replace the rules between two of `types` with `rules`.

        The bond panel only shows the rules for the active structure's
        types. Merging its rows into `bond_map` kept a rule that had been
        deleted, or whose types had been changed, and the next rebuild of
        the panel brought it back.
        """
        present = {str(t) for t in types}
        self.bond_map = {
            key: cutoff
            for key, cutoff in self.bond_map.items()
            if _pair(key) is None or not _pair(key) <= present
        }
        self.bond_map.update(rules)


def _pair(key: str) -> frozenset[str] | None:
    parts = key.split("-")
    return frozenset(parts) if len(parts) == 2 else None


# ----------------------------------------------------------------------
# Disk
# ----------------------------------------------------------------------


def _read_json(path: Path) -> dict[str, Any]:
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f) or {}
    except (OSError, ValueError) as e:
        print(f"Could not read {path}: {e}")
        return {}


def _deep_update(base: dict, update_with: dict) -> None:
    """Merge `update_with` into `base`, nested dictionaries included."""
    for key, value in update_with.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _deep_update(base[key], value)
        else:
            base[key] = value


def _diff(default: dict, current: dict) -> dict:
    """The entries of `current` that are new or differ from `default`."""
    diff = {}
    for key, value in current.items():
        if key not in default:
            diff[key] = value
        elif isinstance(value, dict):
            nested = _diff(default.get(key) or {}, value)
            if nested:
                diff[key] = nested
        elif value != default[key]:
            diff[key] = value
    return diff


def load_defaults() -> dict[str, Any]:
    """The settings shipped with the package, as a dictionary."""
    return _read_json(DEFAULT_CONFIG) if DEFAULT_CONFIG.exists() else {}


def load() -> DisplaySettings:
    """The shipped defaults, overridden by what the user saved."""
    merged = load_defaults()

    # Fall back to the copy earlier versions left inside the package, so
    # that an upgrade does not discard settings already saved.
    user_path = _userdata.config_file()
    if not user_path.exists():
        legacy = _userdata.legacy_config_file()
        if legacy.exists():
            user_path = legacy

    if user_path.exists():
        _deep_update(merged, _read_json(user_path))

    return DisplaySettings.from_dict(merged)


def save(settings: DisplaySettings) -> Path:
    """Write what differs from the shipped defaults to the user's file."""
    user_path = _userdata.config_file()
    with open(user_path, "w", encoding="utf-8") as f:
        json.dump(_diff(load_defaults(), settings.to_dict()), f, indent=4)
    return user_path
