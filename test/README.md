# test/

Only `unit/` is tracked in git. Everything else here is local working
data, kept out of the repository by `.gitignore` because it is either
large or not reproducible from the source.

## `unit/` — the automated suite

262 tests, run on every push by `.github/workflows/tests.yml`.

```bash
pytest
```

from the repository root. `pyproject.toml` restricts collection to
`test/unit`, so that command picks up nothing else.

| file | covers |
|------|--------|
| `test_core.py` | construction, properties, box, masses, charges, elements |
| `test_edit.py` | `EditMixin` — adding, removing, moving, replicating |
| `test_io.py` | readers and writers, and the MDAnalysis/pymatgen converters |
| `test_topology.py`, `test_topology_rules.py` | typing rules and connectivity |
| `test_forcefield.py`, `test_forcefield_database.py` | ff keys, parameters, the bundled database |
| `test_build_*.py` | one file per builder, plus `test_builders.py` for what they share |
| `test_analysis.py` | RDF, density, MSD and diffusion |
| `test_gui_smoke.py` | the interface opens and its actions are wired |

Two things the tests need from outside Python: **Packmol** on `PATH` (the
tests that pack a system skip themselves without it, so a run can look
green while never exercising the builders) and, for `test_gui_smoke.py`,
the `[gui]` extra — it skips when PySide6 or pyvistaqt is missing.

### `unit/data/` — fixtures

Six small files, 56 kB in total: `calcite.cif`, `calcite_ortho.data`,
`h2o.pdb`, `h2o.lt`, `ho.sdf` and `caffeine.lt`. They are in the
repository because a test that cannot find its input is a test that does
not run — several of these are duplicated from `cemd/build/_structures/`
for exactly that reason, so the suite does not depend on package
internals staying where they are.

Most tests need no file at all: they build what they check, from a box of
uniformly distributed atoms whose density is known to Brownian walkers
whose diffusion coefficient is known. Asserting against an analytic value
rather than a stored result is what catches a normalisation being wrong,
and it is how the density-map factor-of-forty error surfaced.

## `legacy/` — old working files, not tracked

Structures and data from earlier work: C₃S and portlandite slabs,
tobermorite variants, C-S-H at several Ca/Si ratios, and `test.dcd`, an
orphaned 297-atom trajectory with no matching topology. Nothing here is
used by the suite.

## `tutorials/` — inputs for the documented walkthroughs, not tracked

- `csh/` — `csh_{12,14,16,18}_reacted.data`, C-S-H relaxed with ReaxFF at
  Ca/Si 1.2 to 1.8. These are the reference structures the Ca/Si and
  H₂O/Si ratios were validated against.
- `slab_caco3/` — includes `caffeine.lt`, the ATB topology used by the
  solvation tutorial (copied into `unit/data/` so the test suite does not
  depend on this directory).
- `caffeine/`

## `traj_solution/` — a reference trajectory, not tracked

A 4344-atom box (1433 SPC/E waters, sodium, sulfate, hydroxide) over
10 001 frames, 100 fs apart. Half a gigabyte, hence untracked. It is what
the analysis module was validated against: the water O–O peak lands at
2.73 Å and the diffusion coefficient at 3.11 × 10⁻⁹ m²/s, both within
reach of the published SPC/E values.

The frame interval is **100 fs**, and it cannot be read from the file:
`universe.trajectory.dt` reports 4.888821 ps, a factor of 49 out. Passing
that to `msd()` yields a diffusion coefficient wrong by the same factor,
with a perfectly straight curve and nothing to suggest a problem.
