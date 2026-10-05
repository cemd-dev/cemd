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

from .rules import NeighborCriterion, TopologyRule

CLAYFF_RULES: list[TopologyRule] = [
    # --- Silicate oxygens ---
    # Bridging Si-O-Si (exactly 2 Si)
    TopologyRule("type O", NeighborCriterion("type Si", 2.0, 2), new_type="Ob"),
    # Bridging Si-O-Al (exactly 1 Si and 1 Al)
    TopologyRule(
        "type O",
        [NeighborCriterion("type Si", 2.0, 1), NeighborCriterion("type Al", 2.0, 1)],
        new_type="Obs",
    ),
    # Non-bridging Si-OH (exactly 1 Si)
    TopologyRule("type O", NeighborCriterion("type Si", 2.0, 1), new_type="Osi"),
    # --- Aluminate oxygens ---
    # Non-bridging Al-OH (exactly 1 Al)
    TopologyRule("type O", NeighborCriterion("type Al", 2.0, 1), new_type="Oa"),
    # --- Water molecules (O bonded to exactly 2 H) ---
    TopologyRule(
        "type O Oa Osi",
        NeighborCriterion("type H", 1.2, 2, "Hw"),
        new_type="Ow",
        bonds=True,
        angles=True,
    ),
    # --- Si-O-H silanols ---
    TopologyRule(
        "type Osi",
        NeighborCriterion("type H", 1.2, 1, "Hsi"),
        new_type="Osih",
        bonds=True,
    ),
    # --- Al-O-H aluminols ---
    TopologyRule(
        "type Oa", NeighborCriterion("type H", 1.2, 1, "Ha"), new_type="Oah", bonds=True
    ),
    # --- Hydroxide (generic O-H) ---
    TopologyRule(
        "type O", NeighborCriterion("type H", 1.2, 1, "Hh"), new_type="Oh", bonds=True
    ),
    # --- Carbonate oxygens ---
    TopologyRule("type O", NeighborCriterion("type C", 1.6, 1), new_type="Oc"),
    # --- Sulfate oxygens ---
    TopologyRule("type O", NeighborCriterion("type S", 1.6, 1), new_type="Os"),
]


CSHFF_RULES = CLAYFF_RULES + [
    TopologyRule(
        center="type Ca",
        neighbors=NeighborCriterion(
            selection="type Ow",
            cutoff=3.2,
            count=1,
            exact_match=False,
        ),
        new_type="Cw",
    ),
]
