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

from dataclasses import dataclass, field

# Pair Coeffs, PairIJ Coeffs, Bond Coeffs, Angle Coeffs, Dihedral Coeffs, Improper Coeffs = force field sections

# BondBond Coeffs, BondAngle Coeffs, MiddleBondTorsion Coeffs, EndBondTorsion Coeffs, AngleTorsion Coeffs, AngleAngleTorsion Coeffs, BondBond Coeffs, AngleAngle Coeffs = class 2 force field sections


@dataclass
class ForceFieldModel:
    """Force field model metadata."""

    name: str
    description: str = ""
    ref: str = ""
    tags: list[str] = field(default_factory=list)
    # Kind of system the model is meant for: "minerals", "carbonates",
    # "ions" or "water" (groups the models in the documentation).
    category: str = ""
    # Combination rule for the LJ pairs the model does not list:
    # "arithmetic" (Lorentz-Berthelot), "geometric", or "none" when every
    # pair is given explicitly and unlisted pairs are Coulomb only.
    mixing: str = ""
    mixing_note: str = ""


@dataclass
class AtomType:
    """Atom type parameters."""

    element: str
    # None means the force field defines no per-type charge (GROMOS and the
    # CHARMM Interface parameter files carry the charges in the molecular
    # topology instead). It is distinct from a genuine 0.0, and callers must
    # not overwrite a system's real charges with it -- see
    # `ForceFieldMixin._update_masses_and_charges`.
    charge: float | None = None
    environment: str = ""
    ref: str = ""
    mass: float | None = None
    model: str | None = None


@dataclass
class LJParams:
    """Lennard-Jones 12-6 parameters."""

    epsilon: float
    sigma: float
    ref: str = ""
    model: str | None = None


@dataclass
class BuckinghamParams:
    """Buckingham potential parameters."""

    a: float
    rho: float
    c: float
    ref: str = ""
    model: str | None = None


@dataclass
class HarmonicBondParams:
    """Harmonic bond parameters."""

    k: float
    r0: float
    ref: str = ""
    model: str | None = None


@dataclass
class MorseBondParams:
    """Morse bond parameters."""

    r0: float
    D: float
    alpha: float
    ref: str = ""
    model: str | None = None


@dataclass
class Class2BondParams:
    """Class2 bond parameters."""

    r0: float
    k2: float
    k3: float
    k4: float
    ref: str = "https://docs.lammps.org/bond_class2.html"
    model: str | None = "COMPASS class2"


@dataclass
class HarmonicAngleParams:
    """Harmonic angle parameters."""

    k: float  # Kcal/(mol·rad²)
    theta0: float  # Degrees
    ref: str = ""
    model: str | None = None


@dataclass
class Class2AngleParams:
    """
    Parameters for the main angle term (E_a) of angle_style class2.

    The potential is: E_a = K2*(θ-θ0)^2 + K3*(θ-θ0)^3 + K4*(θ-θ0)^4
    """

    theta0: float  # Balance angle (in degrees)
    k2: float  # Quadratic coefficient
    k3: float  # Cubic coefficient
    k4: float  # Quartic coefficient
    ref: str = "https://docs.lammps.org/angle_class2.html"
    model: str | None = "COMPASS class2"


@dataclass
class Class2BondBondParams:
    """
    Parameters for the bond-bond term (E_bb) of angle_style class2.

    The potential is: ``E_bb = M (r_ij - r1) (r_jk - r2)``
    """

    m: float  # M (energy/distance^2)
    r1: float  # Equilibrium length of first bond
    r2: float  # Equilibrium length of the second bond
    ref: str = "https://docs.lammps.org/angle_class2.html"
    model: str | None = "COMPASS class2"


@dataclass
class Class2BondAngleParams:
    """
    Parameters for the angle-bond term (E_ba) of angle_style class2.

    The potential is: E_ba = N1*(r_ij -r1)*(θ-θ0) + N2*(r_jk -r2)*(θ-θ0)
    """

    n1: float  # N1 (energy/distance)
    n2: float  # N2 (energy/distance)
    r1: float  # Equilibrium length of first bond
    r2: float  # Equilibrium length of the second bond
    ref: str = "https://docs.lammps.org/angle_class2.html"
    model: str | None = "COMPASS class2"


@dataclass
class HarmonicDihedralParams:
    """Harmonic dihedral parameters."""

    k: float  # Kcal/(mol·rad²)
    d: int
    n: int
    ref: str = "https://docs.lammps.org/dihedral_harmonic.html"
    model: str | None = "Harmonic"

    def __post_init__(self):
        """Validate the dihedral parameters."""
        if self.d not in (1, -1):
            raise ValueError(f"d must be +1 or -1, got {self.d}")
        if self.n < 0:
            raise ValueError(f"n must be >= 0, got {self.n}")


@dataclass
class FourierTerm:
    """Represents an individual term (K, n, delta) in the Fourier sum."""

    k: float  # kcal/mol (can be negative)
    n: int  # int >= 0
    delta: float  # Degrees

    def __post_init__(self):
        """Validates the parameters of an individual term."""
        if self.n < 0:
            raise ValueError(f"n doit être >= 0, reçu : {self.n}")
        if not -360.0 <= self.delta <= 360.0:
            raise ValueError(f"delta doit être entre -360 et 360°, reçu : {self.delta}")


@dataclass
class FourierDihedralParams:
    """Parameters for Fourier-style dihedral (LAMMPS).

    The potential is defined as a sum of 'm' terms:
        E(phi) = Sum_{i=1}^{m} K_i *[1.0 + cos(n_i *phi -d_i)]
    """

    terms: list[FourierTerm]
    ref: str = "https://docs.lammps.org/dihedral_fourier.html"
    model: str = "Fourier"

    def __post_init__(self):
        """Validates that the dihedral contains at least one term (m >= 1)."""
        if not self.terms:
            raise ValueError(
                "Un dièdre Fourier doit contenir au moins un terme (m >= 1)."
            )

    @property
    def m(self) -> int:
        """Returns the total number of terms (the 'm' parameter of LAMMPS)."""
        return len(self.terms)


@dataclass
class Class2DihedralParams:
    """Class2 dihedral parameters."""

    k1: float
    k2: float
    k3: float
    phi1: float
    phi2: float
    phi3: float
    ref: str = "https://docs.lammps.org/dihedral_class2.html"
    model: str | None = "Class2"


@dataclass
class CHARMMDihedralParams:
    """Charm Dihydral Parameters."""

    k: float  # Kcal/(mol·rad²)
    d: int
    n: int
    w: float
    ref: str = "https://docs.lammps.org/dihedral_charmm.html"
    model: str | None = "CHARMM"

    def __post_init__(self):
        """Validate the dihedral parameters."""
        if self.d not in (1, -1):
            raise ValueError(f"d must be +1 or -1, got {self.d}")
        if self.n < 0:
            raise ValueError(f"n must be >= 0, got {self.n}")
        if self.w not in (0, 0.5, 1):
            raise ValueError(f"w must be 0, 0.5 or 1, got {self.w}")


@dataclass
class Class2AngleAngleTorsionParams:
    """Class2 angle-angle-torsion parameters."""

    m: float
    theta1: float
    theta2: float
    ref: str = "https://docs.lammps.org/dihedral_class2.html"
    model: str | None = "Class2"


@dataclass
class HarmonicImproperParams:
    """Harmonic improper parameters."""

    k: float
    chi0: float
    ref: str = "https://docs.lammps.org/improper_harmonic.html"
    model: str | None = "Harmonic"


@dataclass
class Class2ImproperParams:
    """Class2 improper parameters."""

    k: float
    chi0: float
    ref: str = "https://docs.lammps.org/improper_class2.html"
    model: str | None = "Class2"


@dataclass
class CVFFImproperParams:
    """CVFF improper parameters."""

    k: float  # Kcal/(mol·rad²)
    d: int
    n: int
    ref: str = "https://docs.lammps.org/improper_cvff.html"
    model: str | None = "CVFF"

    def __post_init__(self):
        """Validate the dihedral parameters."""
        if self.d not in (1, -1):
            raise ValueError(f"d must be +1 or -1, got {self.d}")
        if self.n < 0:
            raise ValueError(f"n must be >= 0, got {self.n}")


@dataclass
class Class2AngleAngleParams:
    """Class2 angle-angle parameters."""

    m1: float
    m2: float
    m3: float
    theta1: float
    theta2: float
    theta3: float
    ref: str = "https://docs.lammps.org/improper_class2.html"
    model: str | None = "Class2"


@dataclass
class DistanceImproperParams:
    """Distance improper parameters."""

    k2: float
    k4: float
    ref: str = "https://docs.lammps.org/improper_distance.html"
    model: str | None = "Distance"
