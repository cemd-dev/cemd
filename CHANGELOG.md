# Changelog

## Unreleased

### Fixed

- Reading a non-orthogonal crystal (CIF or pymatgen Structure) without
  an axis reindexing placed the atoms in pymatgen's Cartesian frame (c along
  z) inside a box oriented a along x. The fractional coordinates were off by
  up to a quarter of a cell in triclinic kaolinite (COD 9009234, read with
  `refine=False`): half the Al came out five-fold and a quarter of the Si
  three-fold, from bonds across the boundary. Positions are now always taken
  from the fractional coordinates in the box cemd builds (`lattice2vectors`).
  The reindexed branch used `Lattice.from_parameters`, which has the same
  c-along-z convention, and is fixed too. Monoclinic cells with beta close to
  90 degrees were only slightly off. Systems read from a triclinic or
  monoclinic CIF before this fix should be read again.

- `CSHFF_RULES` typed as Cw every Ca with a water O within 3.2 Å
  (since 2026-09-01; positional typing before that). In a defective C-S-H,
  the apical ligand of an intralayer Ca is often an interlayer water at
  2.7 Å, and many interlayer Ca are solvated by OH only: in the ReaxFF-relaxed
  C-S-H of Ca/Si 1.1 to 1.8, 26 to 50 % of the Cw were intralayer Ca, and up
  to a third of the interlayer Ca stayed Ca. Ca/Cw are now typed by position
  with the new `LayerRule`: Ca between the two planes of paired Si of a main
  layer, Cw elsewhere (interlayer, pore, solution). Models built since
  2026-09-01 with `CSHFF_RULES` should be retyped.
- Pair parameters of two types from different force-field models were
  looked up in the model of the first type with the short name of the
  second, then by short name in every model of the database: `cshff2014.ob`
  with `clayff.oh` got CSH-FF's "ob-oh" (a silanol O), and `clayff.cao` with
  `clayff.o*` got the "Ca-Ow" of raiteri2015, a calcite-water model. An
  explicit cross term is now only taken from the model both types belong to,
  and any other pair is mixed from the two self terms. Short names are only
  searched when the assignment names no model.

- The slab generator failed on crystals with mixed or partial site
  occupancy, such as the Ca0.99Mg0.01 sites of a magnesian calcite (COD
  9017359): "species property only supports ordered structures". Each
  disordered site now takes the species with the largest occupancy, both
  for the atoms of the system and for the structure the slab is cut from,
  and reading such a structure warns. The atoms used to take the first
  species listed on the site, whichever had the larger occupancy.
- `SurfaceBuilder` cut the structure cached when the file was read, even
  after the system was edited: delete every carbon from a calcite and the
  (10-14) slab had them back. The cache is now checked against the cell,
  positions and masses of the system before use; typing, which moves
  nothing, keeps it. A crystal edited since it was read, or a slab, is no
  bulk to cut, and `SurfaceBuilder` now refuses it with a `ValueError`
  instead of cutting whatever the system had become. A system read from a
  format that carries no crystal (`.data`, `.xyz`) is still converted as
  it is; `to_pmg` built species from the atom types and failed on
  force-field types ("Can't parse Element or Species from 'Oc'"); it now
  uses the elements.
- Packmol builds of more than 9 999 atoms failed with "invalid literal for
  int()". From 10 000 on, the PDB serial number fills its five columns and
  meets the record name (`HETATM10001`), and the PDB reader split lines on
  whitespace. It now reads the standard columns, and the serial numbers
  Packmol writes past 99 999 in hexadecimal. The PDB writer uses
  hybrid-36 past 99 999, so files of any size read back.
- `msd_profile` kept the MSD of the first block only, and binned it by the
  mean positions of the last block: the profile came from one window of
  `corrlength` frames, with each atom assigned to the bin it occupied at
  the end of the trajectory. Every block now contributes, each atom binned
  by its own mean position over the block. Profiles computed with more
  than one block differ from those of 0.2.0 and earlier, and should be
  recomputed.
- `CSHBuilder.to_cash` put Al on paired silicates of relaxed C-S-H.
  Bridging Si were found by forcing the Si heights into 8 planes and
  keeping the planes smaller than the densest one; once a relaxation
  buckles the sheets, dense and bridging planes merge (on ReaxFF-relaxed
  models, 320 to 556 Si reported instead of 80 to 232), or no plane is
  left and the search crashes. Bridging Si are now identified locally: a
  paired Si has its dimer partner at the same height about 3 Å away, a
  bridging Si sits about 2 Å above the paired Si around it and still has
  a Si at Si-O-Si distance (a monomer from a hydrolysed dimer that drifted
  toward the interlayer is not a bridging site). Pores are
  grouped from the side each bridging Si points to, so pairs face each
  other across the interlayer, never across the CaO sheet, and a bridging
  Si whose partner was removed is kept as a single instead of being
  skipped. Above Ca/Si 1.2 every remaining bridging Si is such a single.
- `to_cash` substituted both bridging Si facing each other across the
  interlayer, which tobermorite 11 A links through an O (Q3): Al-O-Al
  bonds, against Loewenstein's rule (43 of them in a Ca/Si 1.1 model).
  A site is now skipped when an Al already sits at Si-O-Si distance.
  Within a pore, sites were taken in index order, so the Al gathered in
  one replica of a replicated cell; they are now drawn at random, pores
  still being filled in turn.
- A "PairIJ Coeffs" section of a LAMMPS data file was read as "Pair
  Coeffs" (`i j eps sigma` taken as `i eps sigma`, three numbers read as a
  Buckingham potential), under the bare type instead of the "A-B" key the
  rest of the code uses: a file read then written lost every pair
  coefficient. "Pair Coeffs" lines are now keyed "A-A". The writer also
  tested the first two characters of each key instead of its two types
  to choose between the two sections.
- `set_type2atoms` never stored the mass or charge of a new type: it
  tested it against the `masses`/`charges` properties, which already list
  it. The next `set_types` then gave the type a mass of 1.008, so a Si
  turned Al was written with the mass of hydrogen. `set_types` now carries
  the effective masses and charges over, atomic masses included.

### Added

- `AtomicSystem.substitute_atoms()`: replace randomly chosen atoms of one
  type by another, for isomorphic substitution (Al for Si in a mica, Mg for
  Al in a smectite). By default no two substituted atoms share a bridging
  oxygen (Loewenstein's rule), with periodic boundaries; `candidates`
  restricts the choice to a set of sites, one sheet at a time for instance.
  A CIF with mixed occupancy is read with the majority species, so this is
  how the minority one is put back.
- `cemd.analysis.analyze_coordination()`: distribution of the coordination
  number of a cation (Al by default) by oxygen, for an AtomicSystem or frame
  by frame along a trajectory. The 4-, 5- and 6-fold Al fractions compare
  directly with 27Al MAS NMR.

- Force-field model `cshff_duque` (`cshff_duque.toml`): the CSH-FF + ClayFF
  hybrid of Duque-Redondo et al. (doi:10.1039/c8cp00654g) with the charges of Ihaddadene et al.
  (2026), self LJ terms only (every cross term mixed). Ca 1.43, Cw 2.0,
  Si 1.72, Hsi 0.29, Ob -1.14 and Osi = Osih -1.0 before neutralization,
  ClayFF hydroxide, flexible SPC water, Na, Cl, Cs.
- `LayerRule` (`cemd.topology`): retypes atoms inside / outside the main
  layers of a lamellar solid, found from the density of a sheet atom along
  z. Usable in any rule list given to `set_topology`.
- `Splitter.describe_cut(coordinate)`: the side mask, the severed bonds
  and the number of severed ionic contacts of a cut, from a single
  placement of the atoms. With `keep_molecules`, that placement is most
  of the cost, which `find_broken_bonds` and `count_broken_ionic` each
  paid again. It also lists the severed bonds `repair` cannot rebuild.
  `Splitter.best_in_scan(scan)` picks the position `best_coordinate`
  would from a `scan_broken_bonds` table already made.
- `msd_profile(edges=...)`: explicit bin edges, possibly uneven, in place
  of the `delta` bins laid from the lowest coordinate. Aligning the bins
  on each interface of a slit pore puts the bins of both faces at the
  same distances to their surface.
- `CSHBuilder.to_cash(al_si=..., compensation="OH", seed=...)`. `al_si` is
  the Al/Si of the final C-A-S-H (Si counted after the substitution);
  `as_ratio`, a fraction of the initial Si, is kept. `compensation="OH"`
  balances each Si4+ -> Al3+ by removing the nearest hydroxide ion (an O
  with one H and no Si/Al, such as an interlayer Ca-OH), as Duque-Redondo
  et al. (2021) do: Ca and molecular water are untouched, so Ca/(Si+Al)
  keeps its C-S-H value. `seed` makes the site selection reproducible.
- Packmol runs can be cancelled. `run_packmol(cancel=event)`, or any build
  inside `with cancel_packmol_on(event):` (`SolutionBuilder.build`,
  `Splitter.split`, `CSHBuilder.build`, ...), kills Packmol as soon as the
  `threading.Event` is set and raises `PackmolCancelledError`. Packmol is
  now started without a shell, so a missing binary raises
  `FileNotFoundError` instead of "Packmol failed with return code 127".

### Performance

Measured on water boxes with bonds and angles.

- `set_topology("clayff")` / `("cshff")`: 70 s → 1.0 s on 24 000 atoms,
  5 s on 120 000. Each center searched the whole system for its
  neighbors; the time grew with the square of the size. The typing and
  connections are unchanged, checked against the previous engine on dense
  random mixtures.
- Renaming types (`set_types`): 2 s → 21 ms on 24 000 atoms.
- `guess_connections`: 2.3 s → 0.33 s on 120 000 atoms.
- Writing a PDB file: 1.7 s → 0.17 s on 24 000 atoms; the output is
  unchanged.
- Writing a LAMMPS data file: 7.4 s → 0.8 s on 130 000 atoms. The
  sections are formatted without `DataFrame.to_string`, whose output they
  reproduce byte for byte, scientific notation included; the cases left
  out (NaN, infinite values, strings over 50 characters) still go through
  it.
- `DihedralRule` matched every pair of each kind against all the others;
  it now indexes them by their shared atom.
- `msd_profile` accumulates per-bin sums block by block instead of
  keeping the MSD of every atom of every block: its memory no longer grows
  with the length of the trajectory. On 3 216 water oxygens over 300 000
  frames (blocks of 300 frames every 5), 0.2.0 needed about 3 TB and
  stopped out of memory; it now holds a few hundred MB. Overlapping blocks
  share one unwrapping of the trajectory, which makes it about 20 times
  faster (112 s → 5.7 s on 1 000 atoms, 3 000 frames).

### Changed

- `Splitter(repair=True)` rebuilds the bonds a cut severs by hydrolysis,
  explicit ones and `bonds_dict` contacts (the Si-O of ClayFF) alike. A
  hydrogen now goes with its O, so no O-H is cut. Of every other bond to
  an O, the atom that loses its O gets an OH back and the O gets an H:
  Si-O-Si gives two silanols, Si-O...Ca a silanol and a Ca-OH. The caps,
  and the atoms whose covalent environment changed, take the type and
  charge of an atom of the system in the same environment (Osih/Hsi,
  Oh/Hh, Ow/Hw under CSHFF), and explicit bonds and angles where the
  system has them between the same types, under its own labels. Caps
  that would land on an atom tilt away from the bond instead of being
  skipped. A bond without an O (C-C) cannot be hydrolysed: the split
  raises `ValueError` before anything moves. `repair` is now a field of
  `Splitter`, so the counts and scans see the bonds the repair rebuilds;
  `split(repair=...)` still sets it. `repair_report` gives `broken`,
  `capped`, `crowded`, `retyped`, `untyped`, `bonds_added`,
  `angles_added` and `charge_change`; `skipped` is gone. On a C-S-H cut
  through a sheet (319 Si-O severed), the old repair skipped 21 caps,
  left 32 atom pairs closer than 0.9 Å and typed every cap `O`/`H` with
  no charge; the new one leaves no overlap, and 2 caps crowded and
  untyped, reported.
- `density_profile`, `density_map`, `compute_rdf`, `minmax_position` and
  `write_dcd` take the same frame range, `start`/`stop`/`step`, read as a
  Python slice: `stop=None` (the default) runs to the last frame. `end` and
  `compute_rdf`'s `skip` are gone, and passing them raises `TypeError`.
  `end=-1` meant every frame in `density_profile` but dropped the last one
  in `density_map` and `minmax_position`, which reached the slice as
  `[start:-1]`; maps computed with the default range differ from those of
  0.2.0 and earlier by that frame. An empty range raises `ValueError`
  instead of dividing by zero.
- `msd_profile` returns `None` for the coordinate variation, its third
  result, unless called with `coords_var=True`: that one still keeps every
  position of every block in memory.

### GUI

- The Split dialog previews the cut. It shows the structure in a 3D view
  of its own, with the cut plane and, in red, the bonds it would sever,
  counted live with the ionic contacts as the coordinate moves (spinbox,
  slider, or a click on the plot). A scan of every cut position along the
  axis is plotted under the slider; the dialog opens on the cut that
  breaks the fewest bonds, then the fewest ionic contacts, and Best goes
  back to it. The view looks at the cut side on. The coordinate is
  limited to the span of the atoms along the axis.
- The Split dialog offers three cuts: "Plane", "Plane, repaired"
  (`repair=True`) and "Keep molecules whole" (`keep_molecules=True`).
  Repaired, the preview counts the bonds to rebuild and warns about any
  without an O, and a summary of `repair_report` (charge change, crowded
  and untyped caps) follows the split. Kept whole, no bond breaks, and the
  preview halos in orange the atoms dragged across the plane, and counts
  them. Plain cuts through a C-S-H used to be the only choice in the GUI,
  and left the severed Si-O bonds broken.
- `PySide6-Qlementine` 1.5.1.0, a stable release, replaces 1.4.3.0rc7.
- Solution, glass, liquid-layer, droplet and channel builds run in the
  background behind a progress dialog. The window used to freeze for the
  whole Packmol run. The dialog has a Cancel button that stops Packmol.
- Closing the C-(A)-S-H builder during a build stops Packmol; it used to
  run on in the background.
- The surface search in the slab generator runs in the background too.
- Escape and the title-bar button no longer close a progress dialog that
  has no Cancel button. Closing it gave back the window with the build
  still reading the structure.
- Systems of more than 100 000 atoms are drawn from a subset of their
  atoms, and without bonds, while the view is rotated, panned or zoomed,
  and whole once it stops. On an integrated GPU, 500 000 atoms rotate at
  about 18 frames per second instead of 6. Types of 1000 atoms or fewer
  are always drawn whole.
- A glass opened in a tab named after its density; the name was taken for
  a file path, so Save wrote a file called `Glass_<density>gcm3` into the
  working directory.
- A tab with unsaved changes carries a dot. Closing it offers to save
  first, and quitting lists the structures that would be lost. A saved tab
  closes without asking.
- A coordinate or type typed in the atom table reaches the 3D view. The
  edit was written to the table but failed silently before the view, the
  masses and the bonds were updated.
- Display settings (colors, radii, bond rules, scales, background) are one
  object shared by every tab, instead of copies kept in step between each
  plotter, the side panels and the main window.
- The saved atom scale and bond radius are no longer reset to 1.0 and 0.1
  when the window is closed without opening a structure.
- A bond rule deleted from the bond panel, or one whose types were
  changed, no longer comes back on the next rebuild of the panel.
- On the light blue background, the axis labels are dark like the box
  edges; they used to be white.
- Translating atoms can be undone. A failed liquid-layer, droplet or
  channel build no longer leaves an empty step in the undo history.
- Picking an atom color opens one color dialog instead of two.
- The atom-scale and bond-radius sliders apply on every move: a click on
  the groove, the mouse wheel and the arrow keys used to change the slider
  without changing the view. A burst of changes refreshes the view once.
- Changing atom radii or the scale no longer recomputes the bonds.
- Deleting atoms can be undone, and no longer leaves the deleted atoms
  selected.
- Closing a tab frees its 3D view and its structure; they used to stay in
  memory until the application quit.
- The atom table scrolls faster on large systems, and selecting a block of
  atoms in the 3D view is about 15× faster on 60 000 atoms.
- The per-type atom counts in the tooltips were 0 when the types were
  integers.

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
