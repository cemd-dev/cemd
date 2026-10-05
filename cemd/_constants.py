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

AVOGADRO = 6.02214076e23

MASSES_DICT = {
    "H": 1.007947,
    "He": 4.0026022,
    "Li": 6.9412,
    "Be": 9.0121823,
    "B": 10.8117,
    "C": 12.01078,
    "N": 14.00672,
    "O": 15.99943,
    "F": 18.99840325,
    "Ne": 20.17976,
    "Na": 22.989769282,
    "Mg": 24.30506,
    "Al": 26.981538613,
    "Si": 28.08553,
    "P": 30.9737622,
    "S": 32.0655,
    "Cl": 35.4532,
    "Ar": 39.9481,
    "K": 39.09831,
    "Ca": 40.0784,
    "Sc": 44.9559126,
    "Ti": 47.8671,
    "V": 50.94151,
    "Cr": 51.99616,
    "Mn": 54.9380451,
    "Fe": 55.8452,
    "Co": 58.9331955,
    "Ni": 58.69344,
    "Cu": 63.5463,
    "Zn": 65.382,
    "Ga": 69.7231,
    "Ge": 72.6308,
    "As": 74.921602,
    "Se": 78.963,
    "Br": 79.9041,
    "Kr": 83.7982,
    "Rb": 85.46783,
    "Sr": 87.621,
    "Y": 88.905852,
    "Zr": 91.2242,
    "Nb": 92.906382,
    "Mo": 95.962,
    "Tc": 98.0,
    "Ru": 101.072,
    "Rh": 102.905502,
    "Pd": 106.421,
    "Ag": 107.86822,
    "Cd": 112.4118,
    "In": 114.8181,
    "Sn": 118.7107,
    "Sb": 121.7601,
    "Te": 127.603,
    "I": 126.904473,
    "Xe": 131.2936,
    "Cs": 132.90545196,
    "Ba": 137.3277,
    "La": 138.90547,
    "Ce": 140.116,
    "Pr": 140.90765,
    "Nd": 144.242,
    "Pm": 145.0,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.92535,
    "Dy": 162.5,
    "Ho": 164.93032,
    "Er": 167.259,
    "Tm": 168.93421,
    "Yb": 173.04,
    "Lu": 174.967,
    "Hf": 178.49,
    "Ta": 180.94788,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.084,
    "Au": 196.966569,
    "Hg": 200.59,
    "Tl": 204.3833,
    "Pb": 207.2,
    "Bi": 208.9804,
    "Th": 232.03806,
    "Pa": 231.03588,
    "U": 238.02891,
    "D": 2.01410178,
    "T": 3.01604928,
}

VDW_RADII = {
    "H": 1.10,
    "He": 1.40,
    "Li": 1.82,
    "Be": 1.53,
    "B": 1.92,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "Ne": 1.54,
    "Na": 2.27,
    "Mg": 1.73,
    "Al": 1.84,
    "Si": 2.10,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Ar": 1.88,
    "K": 2.75,
    "Ca": 2.31,
    "Sc": 2.11,
    "Ti": 1.95,
    "V": 1.91,
    "Cr": 1.85,
    "Mn": 1.79,
    "Fe": 1.63,
    "Co": 1.63,
    "Ni": 1.63,
    "Cu": 1.40,
    "Zn": 1.39,
    "Ga": 1.87,
    "Ge": 2.11,
    "As": 1.85,
    "Se": 1.90,
    "Br": 1.85,
    "Kr": 2.02,
    "Rb": 3.03,
    "Sr": 2.49,
    "Y": 2.32,
    "Zr": 2.23,
    "Nb": 2.18,
    "Mo": 2.17,
    "Tc": 2.16,
    "Ru": 2.13,
    "Rh": 2.10,
    "Pd": 1.63,
    "Ag": 1.72,
    "Cd": 1.58,
    "In": 1.93,
    "Sn": 2.17,
    "Sb": 2.06,
    "Te": 2.06,
    "I": 1.98,
    "Xe": 2.16,
    "Cs": 3.43,
    "Ba": 2.68,
    "La": 2.43,
    "Ce": 2.42,
    "Pr": 2.40,
    "Nd": 2.39,
    "Pm": 2.38,
    "Sm": 2.36,
    "Eu": 2.35,
    "Gd": 2.34,
    "Tb": 2.33,
    "Dy": 2.31,
    "Ho": 2.30,
    "Er": 2.29,
    "Tm": 2.27,
    "Yb": 2.26,
    "Lu": 2.24,
    "Hf": 2.23,
    "Ta": 2.22,
    "W": 2.18,
    "Re": 2.16,
    "Os": 2.16,
    "Ir": 2.13,
    "Pt": 1.75,
    "Au": 1.66,
    "Hg": 1.55,
    "Tl": 1.96,
    "Pb": 2.02,
    "Bi": 2.07,
    "Po": 1.97,
    "At": 2.02,
    "Rn": 2.20,
    "Fr": 3.48,
    "Ra": 2.83,
    "Ac": 2.47,
    "Th": 2.45,
    "Pa": 2.43,
    "U": 1.86,
    "Np": 2.41,
    "Pu": 2.39,
    "Am": 2.38,
    "Cm": 2.37,
    "Bk": 2.35,
    "Cf": 2.34,
    "Es": 2.33,
    "Fm": 2.31,
    "Md": 2.30,
    "No": 2.29,
    "Lr": 2.27,
    "Rr": 3.03,
}

CHARGES_DICT = {
    "H": 1,
    "Li": 1,
    "Na": 1,
    "K": 1,
    "Rb": 1,
    "Cs": 1,
    "Be": 2,
    "Mg": 2,
    "Ca": 2,
    "Sr": 2,
    "Ba": 2,
    "B": 3,
    "Al": 3,
    "Ga": 3,
    "In": 3,
    "C": 4,
    "Si": 4,
    "Ge": 4,
    "Sn": 2,
    "Pb": 2,
    "N": -3,
    "P": -3,
    "As": -3,
    "Sb": 3,
    "Bi": 3,
    "O": -2,
    "S": -2,
    "Se": -2,
    "Te": -2,
    "F": -1,
    "Cl": -1,
    "Br": -1,
    "I": -1,
    "Sc": 3,
    "Ti": 4,
    "V": 5,
    "Cr": 3,
    "Mn": 2,
    "Fe": 3,
    "Co": 2,
    "Ni": 2,
    "Cu": 2,
    "Zn": 2,
    "Y": 3,
    "Zr": 4,
    "Nb": 5,
    "Mo": 6,
    "Ag": 1,
    "Cd": 2,
    "W": 6,
    "Pt": 4,
    "Au": 3,
    "Hg": 2,
    "La": 3,
    "Ce": 4,
    "Nd": 3,
    "Gd": 3,
    "Yb": 3,
}

INV_MASSES = {float(v): k for k, v in MASSES_DICT.items()}
MASS_KEYS = np.array(list(INV_MASSES.keys()))

#: How far an atom type's mass may sit from a real element's before the type
#: is judged not to be that element. Every genuine element in the bundled
#: force fields matches to better than 0.05 amu, while a GROMOS united atom
#: -- a carbon carrying its apolar hydrogens, CH3 at 15.035 -- sits about
#: 1 amu away from anything real. The gap is wide.
ELEMENT_MASS_TOLERANCE = 0.2


def element_from_mass(mass: float, tolerance: float = ELEMENT_MASS_TOLERANCE) -> str | None:
    """Return the element symbol whose known mass is closest to `mass`, or
    None if even the closest one is farther than `tolerance` amu away."""
    distances = np.abs(MASS_KEYS - mass)
    index = distances.argmin()
    if distances[index] > tolerance:
        return None
    return str(INV_MASSES[MASS_KEYS[index]])


def two_letter_element_from_type(
    atom_type: str, mass: float, tolerance: float = ELEMENT_MASS_TOLERANCE
) -> str | None:
    """Read `atom_type`'s first two characters as an element symbol (e.g.
    GROMOS ``CL`` or CHARMM ``NA+`` as chlorine and sodium), and accept it
    only if that element's known mass agrees with `mass`.

    This disambiguates real two-letter elements -- Cl, Br, Fe, Si, Na... --
    from one-letter-plus-digit force field type names that happen to share
    the same first two capitals, such as GROMOS's ``CE1`` (carbon) or
    ``NR1`` (nitrogen), which the mass check rules out.
    """
    if len(atom_type) < 2:
        return None

    candidate = atom_type[:2].capitalize()
    if candidate in MASSES_DICT and abs(MASSES_DICT[candidate] - mass) <= tolerance:
        return candidate

    return None
