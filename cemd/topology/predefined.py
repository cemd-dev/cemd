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

from .rules import LayerRule, NeighborCriterion, TopologyRule

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


# Ca of the main layers (between the planes of paired tetrahedra) stay Ca,
# every other Ca (interlayer, pore, solution) becomes Cw. A neighbour rule
# ("Ca with a water O within 3.2 A -> Cw", used from 2026-09-01) mistyped
# 25-50 % of the Ca of a defective C-S-H: see LayerRule.
CSHFF_RULES = CLAYFF_RULES + [
    LayerRule("type Ca", "type Si", inside_type="Ca", outside_type="Cw"),
]


# INTERFACE FF silica (Emami et al. 2014, https://doi.org/10.1021/cm500365c),
# types and charges of `iff.toml`. From element types (Si, O, H, Na), it
# writes the Si-O and O-H bonds and every angle the model defines:
#   bulk Si-O-Si oxygen   o_silica_bulk      -0.55
#   silanol Si-O-H        o_silica_silanol   -0.675, h_silica +0.4
#   siloxide Si-O(-)      o_silica_siloxide  -0.9
#   Si                    si_silica          +1.1, or si_silica_siloxide
#                         +0.725 when it bears one siloxide group
#   counter-ion           na_silica          +1.0
# The Emami models never put two siloxide groups on one Si: such a Si (and
# any Si without four oxygens) is left untyped, for the caller to see.
# Water oxygens are farther than 2 A from any Si and stay "O".
_SI_O = 2.0  # Si-O cutoff, r0 = 1.68 A
_O_H = 1.2

IFF_SILICA_RULES: list[TopologyRule] = [
    TopologyRule(
        "type O",
        NeighborCriterion("type Si", _SI_O, 2),
        new_type="o_silica_bulk",
        bonds=True,
        angles=True,
    ),
    TopologyRule(
        "type O",
        [
            NeighborCriterion("type Si", _SI_O, 1),
            NeighborCriterion("type H", _O_H, 1, "h_silica"),
        ],
        new_type="o_silica_silanol",
        bonds=True,
        angles=True,
    ),
    TopologyRule(
        "type O",
        [NeighborCriterion("type Si", _SI_O, 1), NeighborCriterion("type H", _O_H, 0)],
        new_type="o_silica_siloxide",
        bonds=True,
    ),
    TopologyRule(
        "type Si",
        [
            NeighborCriterion("type o_silica_siloxide", _SI_O, 1),
            NeighborCriterion("type o_silica_bulk o_silica_silanol", _SI_O, 3),
        ],
        new_type="si_silica_siloxide",
        angles=True,
    ),
    TopologyRule(
        "type Si",
        NeighborCriterion("type o_silica_bulk o_silica_silanol", _SI_O, 4),
        new_type="si_silica",
        angles=True,
    ),
    TopologyRule("type Na", [], new_type="na_silica"),
]
