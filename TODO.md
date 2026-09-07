# TODO

## Before the first release

Rehearsed on TestPyPI on 2026-09-07: the workflow published 0.1.0 through
trusted publishing, and the package installed and ran from that upload --
builders, Packmol, the LAMMPS round trip and the interface, nine checks.

- [ ] Declare the trusted publisher on **PyPI** (environment `pypi`, the
      name the workflow uses for a release; TestPyPI used `testpypi`).
- [ ] Enable the Zenodo integration **before** tagging: Zenodo only sees
      releases published after it is switched on.
- [ ] Create the release `v0.1.0`. The workflow refuses to publish if the
      tag and `pyproject.toml` disagree, so bump both together or neither.
- [ ] Put the Zenodo concept DOI -- the one that always resolves to the
      latest version -- in `CITATION.cff`, and add the DOI badge.
- [ ] Check `date-released` in `CITATION.cff` still matches the day the
      release actually goes out. It currently says 2026-09-08.

## analysis

The weakest part of the package. `compute_rdf`, `density_profile`,
`density_map`, `msd` and `diffusion_coefficient` are now checked against
values known analytically or from the literature; the rest is not checked
at all.

- [ ] `tcf.py`, `velocities.py` and `util.py` have no tests and no
      validation. Nothing says whether they are right.
- [ ] Rework the module around `MDAnalysis.analysis.base.AnalysisBase`:
      a standard `run(start, stop, step)`, a progress bar, a `results`
      namespace, and one pass over the trajectory instead of one per atom
      type. Accumulate 1D and 2D histograms — a 3D grid costs 7.6 GB per
      atom type on a 100 Å box at the default `bin_size=0.1`.
- [ ] Consider delegating `compute_rdf`, `msd` and `density_profile` to
      `InterRDF`, `EinsteinMSD` and `LinearDensity`. They do the same work,
      are tested by a large community, and would leave us maintaining only
      what is genuinely ours: the silicate analysis, the 2D map, the
      profiles resolved along an axis, and `tcf`.

## Robustness

- [ ] `SurfaceBuilder` builds from `system._pmg_struct`, the pymatgen
      structure cached when the file was read, and that cache is
      invalidated by some operations but not others. Delete every carbon
      from a calcite and ask for a (10-14) surface: the carbon comes back,
      silently. Type the same system first and the cache *is* cleared, so
      it falls through to `to_pmg()` -- which builds `species` from
      `atoms["type"]` and dies on `Can't parse Element or Species from
      'Oc'`. The fix is both halves: rebuild from `elements` rather than
      the force-field types, and drop the cache on every mutation. Until
      then, generate the surface before editing or typing the system,
      which is what the tutorials do.

- [ ] 15 bare `except:` clauses swallow every error. One of them hid the
      `H2O/Si = 0` bug for as long as it existed: an invalid selection
      returned `None` instead of raising.
- [ ] Several analysis functions index `universe.dimensions` without
      checking it exists, so a trajectory carrying no box dies on
      `TypeError: 'NoneType' object is not subscriptable`.
- [ ] `MASSES_DICT` stops at barium, so any heavier element resolves to the
      nearest one that is present. `elements` now refuses rather than
      guessing, but the table itself is still short.
- [ ] Element guessing from a type name mishandles two-letter symbols in
      capitals: `gromos.CL` reads as carbon while carrying chlorine's mass,
      `iff_charmm.NA+` as nitrogen. Trying the two-letter form first and
      validating it against the mass would fix `CL`, `BR`, `FE`, `AR`,
      `SE`, `SI`, `HE`, `NE` and `NA+` without breaking the united atoms.

## Provenance

- [ ] A system does not remember where it came from. `from_cod(9016705)`
      lets a script name its source, but the object it returns cannot
      report one, so a model read back from a `.data` file has lost the
      link to the entry it was built from. A `source` field carried
      through `copy` and `_replace_internals` would close it, and would
      also give `CITATION.cff` something to point at for the structures.

## GUI

- [ ] `test_gui_smoke.py` checks that the window opens and the actions are
      wired. Nothing exercises a dialog or a build.
- [ ] Undo covers the operations that change a structure in place. It does
      not cover the type and connectivity managers.
- [ ] Level of detail while rotating, for large systems.

## Housekeeping

- [ ] 43 ruff findings, mostly pre-existing. Several are false positives on
      Qt conventions (`closeEvent` must stay camelCase). Worth triaging so a
      lint gate can be added to CI.
- [ ] `docs/_build/` is no longer tracked, but `docs/api/generated/` still
      is. Those files are produced by autosummary at build time.

## Done

- Test suite for `AtomicSystem` and the builders — 263 tests, run in CI.
