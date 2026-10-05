# TODO

## analysis

`compute_rdf`, `density_profile`, `density_map`, `msd` and
`diffusion_coefficient` are checked against values known analytically or
from the literature; `tcf`, `velocities` and `util` against synthetic
trajectories and a brute-force reference on `traj_solution`.

- [ ] `end=-1` means "every frame" in `density_profile` but drops the last
      frame in `density_map` (it is a plain `[start:-1]` slice), and
      `compute_rdf` takes `skip` where the others take `start`/`end`. A
      shared `start`/`stop`/`step` would settle it. Rebuilding the module on
      `AnalysisBase` was weighed and left: the speed-up and the
      reproducibility fix (`density_profile` counted frames in dask threads
      sharing one reader) came from reading the trajectory once, in order,
      which needed no framework.

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

## Cement hydrates

- [ ] `CSHBuilder.build`'s bridging-silicate vacancies are charge-neutralized
      globally: `neutralize_csh_charge` (`_silicate_helpers.py`) strips H at
      random from whatever water ended up in the box, unrelated to which
      vacancy produced the imbalance. Reviewers suggest compensating each
      vacancy locally with a Ca2+ instead, bridging the two dangling
      non-bridging O left by the removed tetrahedron, which matches
      EXAFS/PDF Ca coordination data and existing C-S-H molecular models
      (Pellenq et al. 2009, Qomi et al. 2014) better than a diffuse H
      removal, and should matter for Ca/Si-dependent stiffness (Ca-O
      crosslinks between chains vs silanol H-bonds). `_find_symmetric_-
      bridging_pairs` already tracks the paired vacancy positions per pore,
      so the site data needed for placement exists; today this local
      compensation only happens once `min_mcl` blocks further vacancies
      (via `nca_to_add`), everything below that threshold goes through the
      random-H path. A fraction parameter (silanol-cap vs Ca-bridge per
      vacancy) would let the mix be calibrated against measured Q^n
      distribution and density instead of picking one mechanism outright.

      Resolved 2026-09-22: `tolerance` is now computed per layer in
      `fill_csh_interlayers`, from that layer's own packing density
      (`_estimate_packing_tolerance` in `_interlayer_helpers.py`), instead
      of a single hard-coded 2.0. It scales down linearly, past the bulk
      liquid-water density (0.03346 molecule/A^3), from 2.0 A down to a
      1.5 A floor. Measured (not estimated) on the Ca/Si 1.2-1.7 sweep:
      the cliff is sharp and sits almost exactly at bulk water density —
      below it, tolerance=2.0 already converges in seconds; above it, the
      build time explodes (up to ~250 s for a single interlayer) and then
      code 173 appears once density is ~1.7x bulk. A single-layer sweep at
      Ca/Si 1.7 (density ~1.7x bulk) went from 252.6 s + warning at
      tolerance=2.0 to 8.3 s + no warning at 1.8. Full build,
      Ca/Si 1.2-1.7: 570.9 s with 2 warnings before, 28.3 s with 0
      warnings after for the 1.7 case alone; all six models together now
      build in 3 min 7 s with zero Packmol warnings (was: the 1.7 case
      alone took 9 min 33 s).

      Not done: the O-O / O-H distances were only checked right after
      Packmol (nearest-neighbour distance 0.97 A, consistent with an
      intramolecular O-H bond, no pathological overlap), not after the
      ReaxFF minimisation in `react.lmp` as the original plan asked.
      `tolerance` was also not exposed as a `build` parameter -- the fix
      picks it automatically per layer instead, so scripts don't need to
      choose one. The other options considered (smaller margin, multi-pass
      packing, failing loudly / capping `ws_ratio` past a density
      threshold) were not implemented; the adaptive tolerance alone
      removed the warning and the slowdown on the tested range.

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
- [ ] Level of detail while rotating, for large systems.
- [ ] Remove `_TooltipColorFixer` (`main_window.py`) once Qlementine ships
      its fix. Today Qlementine paints the tooltip background dark, but Qt
      resets the tooltip palette to the system default (black text) after
      startup, so tooltips are black on dark grey. The filter re-applies
      the theme colors to each tooltip as it is shown. When Qlementine is
      fixed, drop the class and its `installEventFilter` line in `main()`,
      bump the `PySide6-Qlementine` pin, and check that tooltips are still
      readable.

## Housekeeping

- [ ] 43 ruff findings, mostly pre-existing. Several are false positives on
      Qt conventions (`closeEvent` must stay camelCase). Worth triaging so a
      lint gate can be added to CI.
