# Changelog

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
