.. _tutorial_csh_nanopore:

======================================================
Tutorial: A C-S-H Nanopore Filled with a NaCl Solution
======================================================

This tutorial builds a calcium silicate hydrate (C-S-H) matrix, opens a
nanopore inside it, fills that pore with a sodium chloride solution and
assigns a force field. It is the starting point of the molecular simulations
of chloride transport and binding in cement paste.

The interesting step is *where* to cut. Under CSHFF (as under ClayFF) the
silicate framework carries no explicit bond, so nothing in the topology
stops you from slicing straight through a silicate sheet. CEMD can find
those contacts geometrically, which turns the choice of a cut plane into
something you can measure rather than guess.

.. contents::
   :local:
   :depth: 2

Prerequisites
==============

.. code-block:: python

   from cemd.build import CSHBuilder, SolutionBuilder, Splitter

Packmol must be installed (see :doc:`../installation`).

Step 1: Build the C-S-H matrix
===============================

:class:`~cemd.build.CSHBuilder` starts from a tobermorite model and turns
it into C-S-H: it removes bridging silicates to reach the requested Ca/Si
ratio, then fills the interlayers with water and charge-balancing Ca\
:sup:`2+`.

.. code-block:: python

   builder = CSHBuilder(cs_ratio=1.5, ws_ratio=1.0)
   csh = builder.build(supercell=[4, 1, 1], model="tob11a_merlino.cif")

   print(csh)
   print(csh.box)

.. code-block:: text

   Adding 105 H2O and 15 Ca2+ in layer 1
   Adding 105 H2O and 15 Ca2+ in layer 2
   <AtomicSystem with 1620 atoms, 360 bonds>
   [26.9279995  30.83762169 22.68000031 90.         90.         90.        ]

The interlayers are packed with Packmol, so the exact water positions --
and therefore some of the counts later in this tutorial -- change slightly
from one run to the next.

Step 2: Check the chemistry
============================

:meth:`~cemd.build.CSHBuilder.analyze` reports the stoichiometry and the
polymerization of the silicate chains:

.. code-block:: python

   for key, value in builder.analyze().items():
       print(f"{key}: {value}")

.. code-block:: text

   Ca/(Si+Al): 1.5
   Ca/Si: 1.5
   Al/Si: 0.0
   H2O/(Si+Al): 1.0
   H2O/Si: 1.0
   MCL: 3.0
   Qn_distribution: [ 0.         66.66666667 33.33333333  0.          0.        ]
   Q2(1Al): 0.0

Both ratios are the ones requested. ``MCL: 3.0`` with two thirds of the
silicon in Q\ :sup:`1` and one third in Q\ :sup:`2` is the signature of
dreierketten chains cut down to dimers and pentamers -- the expected result
of removing bridging tetrahedra at Ca/Si = 1.5. ``H2O/Si`` counts the
hydroxides as half a water each (H\ :sub:`2`\ O = H / 2).

Step 3: Assign a topology
==========================

The builder returns crystallographic labels. Reset them to elements, then
let the CSHFF rules assign force-field types:

.. code-block:: python

   csh.set_types_from_elements()
   csh.set_topology("cshff")

   print(csh.atom_types)
   csh.summary()

.. code-block:: text

   ['Ca', 'Cw', 'Hh', 'Hsi', 'Hw', 'O', 'Ob', 'Oh', 'Osi', 'Osih', 'Ow', 'Si']

.. code-block:: text

   <AtomicSystem with 1620 atoms, 361 bonds, 150 angles>

   Box
    a (Å)  b (Å)  c (Å)  α (°)  β (°)  γ (°)
    26.93  30.84  22.68  90.00  90.00  90.00

   Atoms
   type  number      mass  charge
     Ca     160 40.078400     0.0
     Cw     110 40.078400     0.0
     Hh      59  1.007947     0.0
    Hsi       1  1.007947     0.0
     Hw     300  1.007947     0.0
      O       1 15.999430     0.0
     Ob     120 15.999430     0.0
     Oh      59 15.999430     0.0
    Osi     479 15.999430     0.0
   Osih       1 15.999430     0.0
     Ow     150 15.999430     0.0
     Si     180 28.085530     0.0

   Bonds
       type  number
      Hh-Oh      59
      Hsi-O       1
   Hsi-Osih       1
      Hw-Ow     300

   Angles
       type  number
   Hw-Ow-Hw     150

   Total charge: 0.000 e
   Volume: 18.83 nm³
   Density: 2.57 g/cm³

The rules separate structural calcium (``Ca``) from interlayer calcium
(``Cw``), bridging oxygens (``Ob``) from the rest of the silicate oxygens
(``Osi``), silanols (``Osih``/``Hsi``) from the hydroxides of the interlayer
(``Oh``/``Hh``) and from water. Note what the topology contains: **every bond
is an O-H**. Not a single Si-O bond is listed -- that framework is described
by non-bonded interactions in CSHFF.

.. note::

   One silanol oxygen stays ``O``: in this structure its Si-O distance is
   2.24 Å, beyond the 2.0 Å of the typing rules. It is typed as a silanol in
   Step 8. Check ``atom_types`` for such leftovers after every typing.

Step 4: Choose where to cut
============================

This is why the cut plane cannot be picked from the bond list. Cutting at
an arbitrary height looks harmless to the topology and wrecks the
structure anyway:

.. code-block:: python

   naive = Splitter(csh, coordinate=5.0, axis="z", gap_size=25.0)
   print(naive.count_broken_bonds())

.. code-block:: text

   80

Eighty severed Si-O contacts -- none of which appear as a bond.
:meth:`~cemd.build.Splitter.scan_broken_bonds` sweeps the whole cell and
reports the damage as a function of the cut coordinate:

.. code-block:: python

   scanner = Splitter(csh, coordinate=0.0, axis="z", gap_size=25.0)
   scan = scanner.scan_broken_bonds(start=0.0, stop=22.68, step=0.5)
   print(scan.to_string(index=False))

.. code-block:: text

    coordinate  n_broken  n_ionic
           0.0        41       36
           0.5        80       14
           ...
           1.5        80       14
           2.0         0      251
           2.5         0      160
           3.0         0      160
           3.5         0      160
           4.0         0      205
           4.5        80       17
           ...
          12.5        80        6
          13.0        80        6
          13.5         0      160
          14.0         0      160
          14.5         0      160
          15.0         0      160
          15.5         0      254
          16.0        80       14
           ...

The profile maps the material: plateaus at zero are the interlayers, the
values of 80 are the silicate sheets, anything in between clips the chains
partially. ``n_ionic`` counts the Ca-O contacts the plane separates: an
interlayer cut always splits the interlayer calcium from one of the two
walls (160 here), and the edges of the window (2.0, 4.0, 13.5, 15.5) cut
more of them.

.. code-block:: python

   clean = scan[scan["n_broken"] == 0]["coordinate"].tolist()
   print(clean)

.. code-block:: text

   [2.0, 2.5, 3.0, 3.5, 4.0, 13.5, 14.0, 14.5, 15.0, 15.5]

Pick a plane in the middle of one of those windows -- ``z = 14.5`` sits in
the middle of the second one, away from the cell boundary.

Step 5: Define the pore solution
=================================

.. code-block:: python

   pore = SolutionBuilder(density=1.0, counts={"Na": 8, "Cl": 8})
   print(pore)

.. code-block:: text

   <SolutionBuilder>

   ┌─ Composition
   │   density:  1.00 g/cm³
   │   type:     counts
   │             Na: 8
   │             Cl: 8

Only the ions are given explicitly: the builder works out how many water
molecules are needed to reach 1.0 g/cm³ once it knows the volume of the
pore. Here it gives about 0.8 mol/kg of NaCl, in the range of the chloride
concentrations studied for marine exposure and de-icing salts. Polyatomic
species (``HO``, ``CO3``, ``SO4``) can be named the same way, and any other
solute supplied as an :class:`~cemd.core.atomic_system.AtomicSystem` through
``structures=``, as in :ref:`the caffeine tutorial <tutorial_caffeine_solvation>`.

.. note::

   The real pore solution of a Portland cement is an alkaline hydroxide
   solution, mostly KOH with some NaOH, at pH 13-14. NaCl is the usual
   choice of molecular simulations because they study chloride ingress, and
   because Na\ :sup:`+` and Cl\ :sup:`-` have well-tested parameters
   compatible with CSH-FF and ClayFF. A dissolved OH\ :sup:`-` has neither:
   the hydroxyls of these force fields describe OH bound in a solid, and a
   classical model cannot exchange the proton with water.

Step 6: Open the pore and fill it
==================================

:class:`~cemd.build.Splitter` cuts the system, translates the upper part
to open a 25 Å gap, and packs the solution inside it -- ``padding`` keeps
the fluid off the two pore walls:

.. code-block:: python

   splitter = Splitter(csh, coordinate=14.5, axis="z", gap_size=25.0)
   print(splitter.count_broken_bonds())

   result = splitter.add_solution(pore, padding=2.0).split()

.. code-block:: text

   0

Zero broken contacts: the pore opens along an interlayer, so no silicate
chain is cut. ``split()`` returns a new system and leaves ``csh`` as it was,
ready for another cut.

Step 7: Inspect
================

.. code-block:: python

   print(result.box)
   result.summary()

.. code-block:: text

   [26.9279995  30.83762169 47.68000031 90.         90.         90.        ]

.. code-block:: text

   <AtomicSystem with 3307 atoms, 1475 bonds, 707 angles>

   Box
    a (Å)  b (Å)  c (Å)  α (°)  β (°)  γ (°)
    26.93  30.84  47.68  90.00  90.00  90.00

   Atoms
   type  number      mass  charge
     Ca     160 40.078400     0.0
     Cl       8 35.453200     0.0
     Cw     110 40.078400     0.0
     Hh      59  1.007947     0.0
    Hsi       1  1.007947     0.0
     Hw    1414  1.007947     0.0
     Na       8 22.989769     0.0
      O       1 15.999430     0.0
     Ob     120 15.999430     0.0
     Oh      59 15.999430     0.0
    Osi     479 15.999430     0.0
   Osih       1 15.999430     0.0
     Ow     707 15.999430     0.0
     Si     180 28.085530     0.0

   Bonds
       type  number
      Hh-Oh      59
      Hsi-O       1
   Hsi-Osih       1
      Hw-Ow    1414

   Angles
       type  number
   Hw-Ow-Hw     707

   Total charge: 0.000 e
   Volume: 39.59 nm³
   Density: 1.67 g/cm³

The cell grew by exactly the gap size (22.68 + 25 = 47.68 Å). The pore
holds the 8 Na\ :sup:`+`, the 8 Cl\ :sup:`-` and 557 new water molecules
(707 - 150), already typed ``Ow``/``Hw``, while the original interlayer
water is still in place.

.. image:: /_static/images/tuto_csh_nanopore.png
   :alt: C-S-H nanopore filled with NaCl solution, seen along a
   :align: center
   :height: 520px

The view is along *a*: the C-S-H walls above and below, the pore in the
middle with Na\ :sup:`+` (purple) and Cl\ :sup:`-` (cyan). The 2 Å
``padding`` leaves a thin empty layer on each wall, which closes during
equilibration.

Step 8: Assign the force field
==============================

``duque2018`` is the C-S-H force field of Duque-Redondo et al. (CSH-FF with
ClayFF for the interlayer calcium, the hydroxides, the water and the ions),
used for Cs\ :sup:`+` in C-S-H by Duque-Redondo et al. (2018) and for
Na\ :sup:`+` and Cl\ :sup:`-` in C-S-H gel pores by Ihaddadene et al. (2026):

.. code-block:: python

   assignments = {
       "Ca": "duque2018.ca", "Cw": "duque2018.cw", "Si": "duque2018.si",
       "Ob": "duque2018.ob", "Osi": "duque2018.osi",
       "Osih": "duque2018.osih", "O": "duque2018.osih", "Hsi": "duque2018.hsi",
       "Oh": "duque2018.oh", "Hh": "duque2018.hh",
       "Ow": "duque2018.ow", "Hw": "duque2018.hw",
       "Na": "duque2018.na", "Cl": "duque2018.cl",
   }
   result.set_ff_from_database(assignments)
   print(f"Total charge: {result.total_charge:.3f} e")

.. code-block:: text

   Total charge: 109.915 e

The silicate oxygen charges of this model are given *before
neutralization*: the excess is spread over them, as in the original work:

.. code-block:: python

   result.neutralize_charge({"Osi": 1.0, "Ob": 1.0, "Osih": 1.0, "O": 1.0})
   print(f"Total charge: {result.total_charge:.3f} e")
   result.write("csh_nanopore.data")

.. code-block:: text

   Total charge: -0.000 e

Every silicate oxygen received the same -0.183 e. The offset grows with the
interlayer calcium (``Cw`` carries the full +2 e in this model): check that
the resulting charges stay reasonable for your Ca/Si.

Variant: cutting through a sheet
=================================

Sometimes the cut plane is imposed -- a given pore width, a given
orientation -- and no clean window is available. ``repair=True`` then
rebuilds what the cut breaks, by hydrolysis:

.. code-block:: python

   rough = Splitter(csh, coordinate=5.0, axis="z", gap_size=25.0, repair=True)
   print(rough.count_broken_bonds())

   repaired = rough.split()
   print(rough.repair_report)
   print(repaired)

.. code-block:: text

   80
   {'broken': 80, 'capped': 240, 'crowded': 0, 'retyped': 80, 'untyped': 0, 'bonds_added': 80, 'angles_added': 0, 'charge_change': 0.0}
   <AtomicSystem with 1860 atoms, 441 bonds, 150 angles>

Each of the 80 severed Si-O bonds gave a silanol on each side: the Si got
an OH back where its O was, the O an H, three atoms per bond. The new
atoms are typed like the silanols already in the system, and their O-H
get an explicit bond wherever the system has one.

.. note::

   The caps are placed along the broken bond: a starting geometry to
   relax. Their charges come from the system; check ``charge_change``
   and neutralize before running.
