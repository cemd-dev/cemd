"""Tests for cemd.build.GlassBuilder."""

from __future__ import annotations

import pytest

from cemd.build import GlassBuilder

from conftest import requires_packmol

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_non_positive_density_raises():
    with pytest.raises(ValueError, match="density must be positive"):
        GlassBuilder(density=0.0, composition={"SiO2": 1})


def test_empty_composition_raises():
    with pytest.raises(ValueError, match="composition must not be empty"):
        GlassBuilder(density=2.3, composition={})


def test_invalid_formula_raises():
    with pytest.raises(ValueError, match="Invalid chemical formula"):
        GlassBuilder(density=2.3, composition={"NotAFormula!!": 1})


def test_non_positive_coefficient_raises():
    with pytest.raises(ValueError, match="must be positive"):
        GlassBuilder(density=2.3, composition={"SiO2": 0})


def test_charge_imbalance_warns():
    with pytest.warns(UserWarning, match="not charge neutral"):
        GlassBuilder(density=2.3, composition={"Na": 5, "O": 1})


def test_charge_neutral_composition_does_not_warn(recwarn):
    GlassBuilder(density=2.3, composition={"SiO2": 1})
    assert not any("charge neutral" in str(w.message) for w in recwarn.list)


# ---------------------------------------------------------------------------
# Composition math
# ---------------------------------------------------------------------------


def test_get_elemental_composition_from_oxides():
    blueprint = GlassBuilder(density=2.3, composition={"SiO2": 1})
    elemental = blueprint.get_elemental_composition()
    assert elemental["Si"] == pytest.approx(1.0)
    assert elemental["O"] == pytest.approx(2.0)


def test_get_mass_per_formula_unit_matches_manual_calc():
    blueprint = GlassBuilder(density=2.3, composition={"Si": 1, "O": 2})
    # 28.08553 (Si) + 2 * 15.99943 (O)
    assert blueprint.get_mass_per_formula_unit() == pytest.approx(60.084, rel=1e-3)


def test_get_num_formula_units_scales_with_volume():
    blueprint = GlassBuilder(density=2.3, composition={"SiO2": 1})
    small = blueprint.get_num_formula_units(10.0**3)
    large = blueprint.get_num_formula_units(20.0**3)
    assert large > small


def test_is_pure_and_get_components():
    blueprint = GlassBuilder(density=2.3, composition={"SiO2": 1})
    assert blueprint.is_pure()
    assert blueprint.get_components() == ["SiO2"]

    mixed = GlassBuilder(density=2.3, composition={"SiO2": 3, "Al2O3": 2})
    assert not mixed.is_pure()


def test_normalize_sums_to_one():
    blueprint = GlassBuilder(density=2.3, composition={"SiO2": 3, "Al2O3": 1})
    normalized = blueprint.normalize()
    assert sum(normalized.composition.values()) == pytest.approx(1.0)
    # normalize() returns a new instance, doesn't mutate the original.
    assert blueprint.composition == {"SiO2": 3, "Al2O3": 1}


def test_from_elements_and_from_oxides_aliases():
    by_elements = GlassBuilder.from_elements(density=2.3, elements={"Si": 1, "O": 2})
    by_oxides = GlassBuilder.from_oxides(density=2.3, oxides={"SiO2": 1})
    assert by_elements.get_elemental_composition() == pytest.approx(
        by_oxides.get_elemental_composition()
    )


# ---------------------------------------------------------------------------
# build() (requires packmol)
# ---------------------------------------------------------------------------


@requires_packmol
def test_build_glass_from_elements():
    blueprint = GlassBuilder(density=2.3, composition={"Si": 1, "Al": 1, "Na": 1, "O": 4})
    system = blueprint.build(box=[15.0, 15.0, 15.0])

    assert system.num_atoms > 0
    assert set(system.atom_types) == {"Si", "Al", "Na", "O"}


@requires_packmol
def test_build_glass_from_oxides():
    blueprint = GlassBuilder(
        density=2.3, composition={"SiO2": 3, "Al2O3": 2, "Na2O": 2}
    )
    system = blueprint.build(box=[15.0, 15.0, 15.0])

    assert system.num_atoms > 0
    assert set(system.atom_types) == {"Si", "Al", "Na", "O"}
    # Roughly 2 O per Si (SiO2 stoichiometry dominates), sanity bound only.
    assert system.get_count("O") > system.get_count("Si")
