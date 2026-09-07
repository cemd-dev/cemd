# TODO

## Before the first release

- [ ] Declare the trusted publishers on PyPI and TestPyPI, rehearse on
      TestPyPI, then publish.
- [ ] Enable the Zenodo integration **before** tagging: Zenodo only sees
      releases published after it is switched on.
- [ ] Update `date-released` in `CITATION.cff` to the tag's date, and put
      the concept DOI in once it exists.
- [ ] Replace the `git+https://` install line with `pip install cemd` in the
      README and the installation guide, and add the PyPI version badge --
      it reads "package not found" until the project exists.

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
