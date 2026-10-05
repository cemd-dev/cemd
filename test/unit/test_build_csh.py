"""Tests for cemd.build.CSHBuilder and cemd.build.AFBuilder.

`CSHBuilder.build()`/`build_pycsh()` and `AFBuilder.build_aft()`/
`build_afm()` run a full Packmol-based construction pipeline that takes
minutes even for small supercells (confirmed manually); they are not
exercised here to keep the suite fast. Everything reachable without
running Packmol -- validation, `from_system`, `analyze`, `to_cash` -- is
covered with real (non-mocked) systems.
"""

from __future__ import annotations

import pandas as pd
import pytest

from cemd import AtomicSystem
from cemd.build import AFBuilder, CSHBuilder
from cemd.build.cement_hydrates._silicate_helpers import (
    calculate_csh_modifiers,
    neutralize_csh_charge,
)


def _make_silicate_chain() -> AtomicSystem:
    """Two corner-sharing SiO4 tetrahedra (a Q1-Q1 dimer)."""
    atoms = pd.DataFrame(
        {
            "type": ["Si", "Si", "O", "O", "O", "O", "O", "O"],
            "charge": [0.0] * 8,
            "x": [0.0, 3.2, 1.6, -1.6, 0.0, 0.0, 3.2, 3.2],
            "y": [0.0, 0.0, 0.0, 0.0, 1.6, -1.6, 1.6, -1.6],
            "z": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        },
        index=range(1, 9),
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
            "masses": {"Si": 28.085, "O": 15.999},
            "charges": {},
        }
    )
    system.add_bond([1, 3])  # bridging O between the two Si
    system.add_bond([2, 3])
    system.add_bond([1, 4])
    system.add_bond([1, 5])
    system.add_bond([1, 6])
    system.add_bond([2, 7])
    system.add_bond([2, 8])
    return system


# ---------------------------------------------------------------------------
# CSHBuilder: validation
# ---------------------------------------------------------------------------


def test_negative_cs_ratio_raises():
    with pytest.raises(ValueError, match="cs_ratio must be positive"):
        CSHBuilder(cs_ratio=-1.0, ws_ratio=1.0)


def test_negative_ws_ratio_raises():
    with pytest.raises(ValueError, match="ws_ratio must be >= 0"):
        CSHBuilder(cs_ratio=1.5, ws_ratio=-1.0)


def test_bare_constructor_allows_none_ratios():
    # `from_system` relies on CSHBuilder() with no ratios being valid.
    builder = CSHBuilder()
    assert builder.cs_ratio is None
    assert builder.ws_ratio is None


def test_build_rejects_unsupported_model():
    builder = CSHBuilder(cs_ratio=1.5, ws_ratio=1.0)
    with pytest.raises(ValueError, match="Unsupported model"):
        builder.build(model="not_a_real_model.cif")


def test_build_rejects_min_mcl_below_two():
    builder = CSHBuilder(cs_ratio=1.5, ws_ratio=1.0)
    with pytest.raises(ValueError, match="min_mcl must be >= 2"):
        builder.build(min_mcl=1.0)


# ---------------------------------------------------------------------------
# CSHBuilder: from_system / analyze / repr
# ---------------------------------------------------------------------------


def test_from_system_runs_initial_analysis():
    system = _make_silicate_chain()
    builder = CSHBuilder.from_system(system)

    assert builder._system is system
    assert builder._analysis is not None
    assert "error" not in builder._analysis
    # A Q1-Q1 dimer: each Si has exactly 1 bridging neighbor.
    assert builder._analysis["Qn_distribution"][1] == pytest.approx(100.0)


def test_analyze_without_system_raises():
    builder = CSHBuilder(cs_ratio=1.5, ws_ratio=1.0)
    with pytest.raises(RuntimeError, match="No system available"):
        builder.analyze()


def test_analyze_matches_from_system_analysis():
    system = _make_silicate_chain()
    builder = CSHBuilder.from_system(system)
    result = builder.analyze()
    assert result == builder._analysis


def test_analyze_silicates_raises_on_an_invalid_selection():
    # A malformed selection used to be swallowed and turned into `None`,
    # surfacing later as a misleading "types must be defined" error -- or,
    # for the hydrogens, as silent zero ratios.
    from MDAnalysis.exceptions import SelectionError

    from cemd.analysis.silicates import analyze_silicates

    with pytest.raises(SelectionError):
        analyze_silicates(_make_silicate_chain(), types_map={"si_types": "Si and and"})


def test_repr_includes_system_and_analysis():
    system = _make_silicate_chain()
    builder = CSHBuilder.from_system(system)
    text = repr(builder)
    assert "CSHBuilder" in text
    assert "System       : 8 atoms" in text
    assert "Qⁿ dist" in text


def test_repr_without_system_says_none():
    builder = CSHBuilder(cs_ratio=1.5, ws_ratio=1.0)
    assert "System       : None" in repr(builder)


# ---------------------------------------------------------------------------
# CSHBuilder: to_cash (Si -> Al substitution)
#
# Regression coverage: `CSHBuilder(system)` used to bind `system`
# positionally to the `cs_ratio` dataclass field instead of wrapping it
# (must use `.from_system(system)`), and `to_cash()` takes no positional
# `system` argument -- both crashed every GUI call site using this API
# (see cemd/gui/ui/build.py fixes).
# ---------------------------------------------------------------------------


def test_to_cash_with_no_eligible_bridging_pairs_is_a_graceful_noop():
    # Regression test: `to_cash` finds substitution candidates via
    # `_find_symmetric_bridging_pairs`, which needs real interlayer-pore
    # geometry; a minimal isolated Si-O-Si dimer (no such geometry) used
    # to make `substitute_si_by_al` return an empty id list, which then
    # crashed deep inside `set_type2atoms` with a confusing pyarrow
    # IndexError instead of leaving the system untouched.
    system = _make_silicate_chain()
    builder = CSHBuilder.from_system(system)

    with pytest.warns(UserWarning, match="No indices given"):
        result = builder.to_cash(as_ratio=1.0)

    assert isinstance(result, AtomicSystem)
    assert result.get_count("Si") == 2
    assert result.get_count("Al") == 0


def test_to_cash_without_system_raises():
    builder = CSHBuilder(cs_ratio=1.5, ws_ratio=1.0)
    with pytest.raises(RuntimeError, match="No system available"):
        builder.to_cash(as_ratio=0.1)


def test_to_cash_no_silicon_raises():
    atoms = pd.DataFrame(
        {"type": ["O"], "charge": [0.0], "x": [0.0], "y": [0.0], "z": [0.0]},
        index=[1],
    )
    system = AtomicSystem(
        {
            "atoms": atoms,
            "box": [10.0, 10.0, 10.0, 90.0, 90.0, 90.0],
            "masses": {"O": 15.999},
            "charges": {},
        }
    )
    builder = CSHBuilder.from_system(system)
    with pytest.raises(ValueError, match="No silicon atoms"):
        builder.to_cash(as_ratio=0.1)


# ---------------------------------------------------------------------------
# AFBuilder: validation only (build_aft/build_afm run Packmol)
# ---------------------------------------------------------------------------


def test_afbuilder_requires_ws_ratio():
    with pytest.raises(ValueError, match="ws_ratio must be >= 0"):
        AFBuilder(ws_ratio=None)


def test_afbuilder_negative_ws_ratio_raises():
    with pytest.raises(ValueError, match="ws_ratio must be >= 0"):
        AFBuilder(ws_ratio=-5.0)


def test_afbuilder_repr():
    builder = AFBuilder(ws_ratio=10.0)
    assert "ws_ratio   : 10.000" in repr(builder)


def test_afbuilder_build_rejects_unknown_structure_type():
    builder = AFBuilder(ws_ratio=10.0)
    with pytest.raises(ValueError, match="Unknown structure_type"):
        builder.build(structure_type="not_a_real_type")


def test_afbuilder_build_af_rejects_bad_supercell_shape():
    builder = AFBuilder(ws_ratio=10.0)
    with pytest.raises(ValueError, match="supercell must have 3 elements"):
        builder.build_aft(supercell=[1, 1])


# ---------------------------------------------------------------------------
# calculate_csh_modifiers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cs_ratio, expected_n_ca_added, expected_n_si",
    [
        (1.4, 12, 180),
        (1.5, 30, 180),
        (1.6, 48, 180),
        (1.8, 84, 180),
    ],
)
def test_calculate_csh_modifiers_hits_the_target_ratio_exactly(
    cs_ratio, expected_n_ca_added, expected_n_si
):
    # Regression test: the calcium count was floored straight off a float
    # product, and 1.4 * 180 is 251.99999999999997 in binary -- so one
    # calcium was dropped and the built system came out at Ca/Si = 1.3944
    # instead of the 1.4 asked for.
    n_si_removed, n_ca_added, _ = calculate_csh_modifiers(240, 240, cs_ratio, 3.0)
    n_si = 240 - n_si_removed

    assert n_si == expected_n_si
    assert n_ca_added == expected_n_ca_added
    assert (240 + n_ca_added) / n_si == pytest.approx(cs_ratio)


def test_calculate_csh_modifiers_returns_an_integer_calcium_count():
    # np.floor returns a float, which propagated into the water count and
    # made the build log read "15.0 Ca2+".
    _, n_ca_added, _ = calculate_csh_modifiers(240, 240, 1.5, 3.0)
    assert isinstance(n_ca_added, int)


def test_calculate_csh_modifiers_adds_no_calcium_below_the_vacancy_cap():
    # Below Ca/Si ~ 1.33 the target is reached by removing bridging
    # silicates alone, so H2O/Si comes back exactly as requested.
    _, n_ca_added, _ = calculate_csh_modifiers(240, 240, 1.2, 3.0)
    assert n_ca_added == 0


# ---------------------------------------------------------------------------
# neutralize_csh_charge
# ---------------------------------------------------------------------------


def _make_csh_like_system(n_si: int, n_ca: int, ws_ratio: float) -> AtomicSystem:
    """A silica framework that is charge-neutral on its own (4*n_si ==
    2*n_o), plus `n_ca` inserted Ca2+ and the water CSHBuilder.build would
    add for them (`round(n_si * ws_ratio) + n_ca` molecules, one extra per
    Ca2+). Water is named O1/H1/H2, like the real h2o.lt structure, not
    bare O/H -- that distinction is exactly what the selection bug in
    `neutralize_csh_charge` used to miss.
    """
    n_o_framework = 2 * n_si
    n_water = round(n_si * ws_ratio) + n_ca

    rows = [{"type": "Si", "charge": 0.0} for _ in range(n_si)]
    rows += [{"type": "O", "charge": 0.0} for _ in range(n_o_framework)]
    rows += [{"type": "Ca", "charge": 0.0} for _ in range(n_ca)]
    for _ in range(n_water):
        rows.append({"type": "O1", "charge": 0.0})
        rows.append({"type": "H1", "charge": 0.0})
        rows.append({"type": "H2", "charge": 0.0})

    n = len(rows)
    atoms = pd.DataFrame(rows, index=range(1, n + 1))
    atoms["x"] = [float(i) for i in range(n)]
    atoms["y"] = 0.0
    atoms["z"] = 0.0

    return AtomicSystem(
        {
            "atoms": atoms,
            "box": [200.0, 200.0, 200.0, 90.0, 90.0, 90.0],
            "masses": {
                "Si": 28.085,
                "O": 15.999,
                "Ca": 40.078,
                "O1": 15.999,
                "H1": 1.008,
                "H2": 1.008,
            },
            "charges": {},
        }
    )


def test_neutralize_csh_charge_neutralizes_and_preserves_h2o_si(tmp_path):
    # Regression test: `neutralize_csh_charge` used to select `name H` /
    # `name O`, which matches nothing on h2o.lt's O1/H1/H2 water, so no
    # hydrogen was stripped and the box came out charged. With the fixed
    # element-based selection, the one extra water molecule added per
    # inserted Ca2+ loses exactly the two hydrogens needed to neutralize
    # that Ca2+, so the built system ends up both neutral and back at the
    # requested H2O/Si.
    n_si, n_ca, ws_ratio = 4, 2, 1.5
    system = _make_csh_like_system(n_si, n_ca, ws_ratio)
    pdb_path = tmp_path / "csh.pdb"
    system.write(str(pdb_path))

    neutralize_csh_charge(str(pdb_path))

    result = AtomicSystem.from_file(str(pdb_path))
    counts = result.atoms["type"].value_counts()
    n_h = counts.get("H1", 0) + counts.get("H2", 0)
    n_o = counts.get("O", 0) + counts.get("O1", 0)
    n_si_out = counts.get("Si", 0)
    n_ca_out = counts.get("Ca", 0)

    total_charge = 2 * n_ca_out + n_h + 4 * n_si_out - 2 * n_o
    assert total_charge == 0
    assert (n_h / 2) / n_si_out == pytest.approx(ws_ratio)
