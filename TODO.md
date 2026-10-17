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

- [ ] Check the O-O and O-H distances in the C-S-H interlayers after the
      ReaxFF minimisation in `react.lmp`. They were only checked right after
      Packmol (nearest neighbour 0.97 A, no overlap), when the packing
      tolerance became adaptive.

## Performance

- [ ] Writing a LAMMPS data file formats each section with
      `DataFrame.to_string`: 2.7 s for 120 000 atoms, most of it in pandas'
      column-width pass. Formatting the columns with numpy would be several
      times faster, but the output has to stay byte-identical (LAMMPS does
      not care, diffs between versions of a model do).

## Provenance

- [ ] A system does not remember where it came from. `from_cod(9016705)`
      lets a script name its source, but the object it returns cannot
      report one, so a model read back from a `.data` file has lost the
      link to the entry it was built from. A `source` field carried
      through `copy` and `_replace_internals` would close it, and would
      also give `CITATION.cff` something to point at for the structures.

## GUI

- [ ] `test_gui_smoke.py` covers the window, the side panels, the atom
      table, undo, unsaved changes and the background runner. Nothing
      exercises a build dialog or the 3D view, which needs a GL context.
- [ ] Level of detail while rotating, for large systems.
- [ ] Packmol builds cannot be cancelled. `run_packmol` starts Packmol with
      `subprocess.run`, so nothing outside it can stop the process; a
      cancel button needs it to hand the process out, or to take an event
      to wait on.
- [ ] The surface search in `SurfaceDialog.on_search` (pymatgen slab
      generation) still runs on the GUI thread. `run_in_background`
      (`logic/build.py`) is the tool for it.
- [ ] Remove `_TooltipColorFixer` (`main_window.py`) once Qlementine's fix
      covers the GNOME platform themes. Qlementine paints the tooltip
      background dark (#404040), but with the `gtk3` or `xdgdesktopportal`
      Qt platform theme -- the default under GNOME -- the tooltip text
      stays black, at a 2:1 contrast. Qlementine 1.5.1 ("Fix tooltip
      foreground color") fixes it only with the `generic` platform theme
      and offscreen. Re-test with each new release: drop the filter and its
      `installEventFilter` line in `main()`, then show a tooltip under
      GNOME with `QT_QPA_PLATFORMTHEME` unset and check the text is light.

## Housekeeping

- [ ] 44 ruff findings (`cemd` and `test`). Several are false positives on
      Qt conventions (`closeEvent` must stay camelCase). Worth triaging so a
      lint gate can be added to CI.
- [ ] The `cemd` conda environment has pyvista 0.49 and vtk 9.7 instead of
      the pinned 0.47.2 and 9.6.2, and the GUI fails there as soon as a
      structure opens (`SetInteractor` in `AtomicPlotter.add_axes`).
      `environment.yml` should install the pinned versions; `cemd_ui` has
      them.
