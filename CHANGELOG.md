# Changelog

## 0.3.0 — 2026-11-13

### Breaking changes

- Force-field models renamed after their paper: `cshff_duque` -> `duque2018`,
  `sperinck` -> `sadat2016`, `wang2008` -> `wang2009`.
- Cross pairs missing from the database are mixed by the rule each model
  declares (`mixing` in its TOML), not always Lorentz-Berthelot; pairs of
  pair-specific models (`raiteri2015`, `solcon`, `freitas2015`, ...) left
  unlisted stay Coulomb only. `apply_pair_mixing_rules()` defaults to
  `rule=None` accordingly; pass `rule=` to keep the 0.2.0 behaviour.
- Builders return a new system and leave their input unchanged:
  `Splitter.split()` no longer cuts in place, `CSHBuilder.from_system()` copies.
- `density_profile`, `density_map`, `compute_rdf`, `minmax_position` and
  `write_dcd` take `start`/`stop`/`step` as a Python slice; `end` and `skip`
  are gone.
- `Splitter.repair_report` no longer has `skipped`; `msd_profile` returns
  `None` as third result unless `coords_var=True`.

### Results to redo

Results from earlier versions are wrong in these cases:

- **C-S-H typed with `CSHFF_RULES` since 2026-09-01**: Ca/Cw was decided by
  the nearest water and mistyped 26-50 % of the Cw. Ca/Cw are now typed by
  position in the layer (`LayerRule`). Retype.
- **Systems typed with `cshff2014`**: the o/ob charges were swapped. Reassign
  the charges.
- **Triclinic or monoclinic CIFs**: atoms were placed in the wrong Cartesian
  frame (up to a quarter of a cell in kaolinite). Read them again.
- **`msd_profile` with several blocks**: only the first block was used.
  Recompute.
- **`density_map` / `minmax_position` with the default range** dropped the
  last frame.

### Fixed

- Force-field database checked against the source papers: `spc_fw` (water O
  never loaded), `raiteri2015` (O-H bond stored as Ow-Ow), `clayff` (Mg sigma,
  missing `ohs` LJ), `iff` (silanol/aluminol terms, swapped silica charges),
  `sadat2016` (SPC/Fw twice too stiff), `tip4p_ew` (TIP4P charges),
  `cshff2014` (charges), `solcon` (replaced by Armstrong et al. 2023).
- Pair parameters of types from two different models are no longer taken from
  an unrelated model by short name (`clayff.cao`-`clayff.o*` got raiteri2015's
  Ca-Ow); a cross term comes only from a model both types belong to.
- C-A-S-H: `to_cash()` finds bridging Si locally, so it works on relaxed,
  buckled sheets; it no longer creates Al-O-Al across the interlayer, spreads
  Al at random over the replicas, and accepts a C-S-H straight from `build()`.
- `analyze_silicates`: mean chain length follows Richardson's formula with
  bridging Al (new `"Q2(1Al)"` entry).
- `set_type2atoms` / `substitute_atoms` update the labels of the bonds and
  angles around retyped atoms, and keep the new type's mass and charge (a Si
  turned Al had the mass of H).
- PDB: bond types from CONECT were shifted by one (every Packmol output);
  files over 9 999 atoms are read and written (hybrid-36).
- LAMMPS data: mixed bond styles written and read in hybrid form;
  "PairIJ Coeffs" read correctly; mixed-length coefficient lines read.
- `SurfaceBuilder` handles partially occupied sites (majority species), no
  longer cuts a stale cached crystal after edits, and refuses an edited
  crystal or a slab.
- GUI builder: H2O/Si was 0 for a C-S-H straight from the builder.

### Added

- Force-field models `clayff_moh` (Cygan 2021), `afff` (AFt/AFm/C3A, Honorio
  2020), `freitas2015` (clinker, CAS glasses), `madrid2019` (Zeron 2019),
  `joung2008` (alkali and halide ions), `duque2018`; `iff` gains silica,
  calcium sulfates, portlandite and tobermorite 14 A. Morse bonds
  (`[bond.morse]`); `category` and `mixing` in every `[model]`.
- `CSHBuilder.to_cash(al_si=..., compensation="OH", seed=...)`.
- `AtomicSystem.substitute_atoms()`: random isomorphic substitution,
  respecting Loewenstein's rule.
- `analyze_coordination()`: cation coordination distribution, to compare with
  27Al NMR.
- `LayerRule` for typing by position in a lamellar solid.
- `Splitter(repair=True)` rebuilds severed bonds by hydrolysis (silanols,
  Ca-OH) with typed caps, and reports what it did; `Splitter.describe_cut()`
  and `best_in_scan()`.
- `msd_profile(edges=...)` for explicit, uneven bins.
- Packmol runs can be cancelled (`run_packmol(cancel=...)`,
  `cancel_packmol_on`).
- Docs: tutorials "From C-S-H to C-A-S-H" and "Kaolinite with ClayFF-MOH",
  a visualisation page, API pages for `cemd.topology` and `cemd.forcefield`;
  force-field tables generated from the database.

### Performance

- `set_topology` on 24 000 atoms: 70 s -> 1 s (5 s on 120 000).
- `set_types` 100x, `guess_connections` 7x, PDB writing 10x and LAMMPS data
  writing 9x faster.
- `msd_profile` memory no longer grows with the trajectory (3 TB -> a few
  hundred MB on 300 000 frames), and it runs about 20x faster.

### GUI

- Split dialog: live 3D preview of the cut and the bonds it severs, a scan of
  every position that opens on the best cut, and the three cut modes (plane,
  repaired, molecules kept whole).
- Long builds (solution, glass, liquid layer, droplet, channel, surface
  search) run in the background with a Cancel button that stops Packmol.
- Systems over 100 000 atoms stay fluid while rotating (about 3x faster);
  faster atom table and selection.
- Unsaved tabs are marked and saving is offered on close or quit; atom
  translation and deletion can be undone; closing a tab frees its memory.
- Fixes: table edits now reach the 3D view, glass tabs no longer save to a
  stray file, settings are no longer reset on close, deleted bond rules stay
  deleted, slider clicks and wheel apply, and more.
- `PySide6-Qlementine` 1.5.1.0 (stable) replaces 1.4.3.0rc7.

## 0.2.0 — 2026-10-05

### Breaking changes

- Force-field entries are keyed by the database file name: `clayff.st`,
  `cshff2014.si`, `spc.ospc`. The model prefix is case-insensitive, so
  `ClayFF.st` and `SPC.ospc` still work. Two kinds of type name changed:
  ClayFF's `h_star` / `o_star` are now `h*` / `o*`, and the CSHFF types
  dropped their suffix (`si_cshff2014` → `si`).
- `AtomicSystem.write(path, **kwargs)` passes its keyword arguments straight
  to the writer chosen by the extension. An argument that writer does not
  know raises `TypeError`, and the PDB writer takes none.
- `dask` is no longer a dependency.

### Added

- `AtomicSystem.neutralize_charge()` spreads the excess charge over every
  atom, or over chosen types with relative weights. `set_charges` and
  `set_ff_from_database` warn when the system is left charged.
- `Splitter(keep_molecules=True)` moves whole covalent units instead of a
  half-space, so opening a pore in a disordered solid breaks no covalent
  bond. `coordinate=None` picks the cut that breaks the fewest bonds and
  then the fewest ionic contacts (`ionic_dict`, `best_coordinate`,
  `count_broken_ionic`).
- `com(..., unwrap=True)` and `shift2com(..., unwrap=True)` make the
  reference group whole across the periodic boundaries before taking its
  centre of mass. Use it for a slab the boundary cuts through. Without it,
  that slab's centre lands in the middle of the pore.
- `.lmp` is accepted as a LAMMPS data file extension by `write`.
- The element mass table now runs up to bismuth, plus Th, Pa and U.

### Fixed

- A force-field key written with the display name (`SPC.ospc`,
  `CSHFF2014.ca`) silently took the mass and charge of another model that
  defines a type with the same short name (IFF-CVFF, ClayFF).
- `density_profile` could return different profiles for two identical calls:
  the frames were counted by threads sharing one trajectory reader. It now
  reads the trajectory once, for all types together, and is faster.
- `density_map` no longer keeps every position in memory until the end.
- Analyses run on a trajectory without a periodic box raise a clear
  `ValueError` instead of a `TypeError`, or of an `inf` RDF.
- `analyze_silicates` reported H2O/Si = 0: its hydrogen selection missed
  water hydrogens, and a bare `except` hid the error. The other bare
  `except` clauses in the package are gone too.
- LAMMPS data files store charges with 10 decimals. With 6, rounding left a
  neutralised system about 1e-3 e off neutral.
- CIF/pymatgen import: a reindexed cell (a, b, c reordered) mirrored or
  distorted the structure. The lattice and the fractional coordinates are
  now permuted together.
- `tcf.bondcorr` counted each atom paired with itself when both types are
  the same, so an Ow–Ow curve plateaued about 0.16 too high instead of
  decaying. It also ignored `nblocks`, and a `corrlength` too long for the
  trajectory raised `IndexError` instead of `ValueError`. The box is now
  read at every frame, which matters in NPT.
- `velocity_profile` returned noise: speeds were paired with the wrong
  atoms, and a list of types selected nothing. It now returns the mean speed
  `|v|` per z bin, and accepts a string or a list of types.
- `mean_pos` returned x, y, z as object columns and failed when called twice
  on the same universe. `minmax_position` raised `UnboundLocalError` for an
  unknown axis instead of `ValueError`.
- `AtomicSystem.copy()` lost the `.box` attribute.
- `set_ff_from_database` dropped stale bond and angle keys when atoms are
  re-typed, and no longer warns about missing class2 cross-terms for class1
  force fields such as ClayFF.
- Element guessing from GROMOS and CHARMM type names: `CL`, `BR`, `FE`, `SI`
  and `NA+` were read as one-letter elements.
- C-S-H builder: the Packmol tolerance adapts to each interlayer's packing
  density. A full Ca/Si 1.2–1.7 sweep goes from about 10 min with Packmol
  warnings to 3 min without. The periodic seam is moved into a dense
  silicate sheet.
- ClayFF/CSHFF typing rules: the Si–O and Al–O cutoffs for bridging and
  non-bridging oxygens went from 1.85 to 2.0 Å, so stretched bonds are
  still detected.

### GUI

- Undo now covers the type and connectivity managers, neutralisation and
  type guessing.
- Tooltips are readable with the dark theme.
