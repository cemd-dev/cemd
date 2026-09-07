.. _tutorial_csh_nanopore:

==================================================================
Tutorial: A C-S-H Nanopore Filled with an Alkaline Pore Solution
==================================================================

This tutorial builds a calcium silicate hydrate (C-S-H) matrix, opens a
nanopore inside it, and fills that pore with the alkaline solution that
sits in real cement paste -- here a sodium hydroxide solution.

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

   Adding 105.0 H2O and 15.0 Ca2+ in layer 1
   Adding 105.0 H2O and 15.0 Ca2+ in layer 2
   <AtomicSystem with 1680 atoms, 420 bonds>
   [26.93 30.84 22.68 90. 90. 90.]

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
   H2O/(Si+Al): 1.1666666666666667
   H2O/Si: 1.1666666666666667
   MCL: 3.0
   Qn_distribution: [ 0.  66.67  33.33  0.  0.]

The Ca/Si ratio is exactly the one requested. ``MCL: 3.0`` with two thirds
of the silicon in Q\ :sup:`1` and one third in Q\ :sup:`2` is the signature
of dreierketten chains cut down to trimers -- the expected result of
removing bridging tetrahedra at Ca/Si = 1.5.

.. note::

   ``H2O/Si`` comes out slightly above the requested ``ws_ratio`` here.
   Up to Ca/Si ≈ 1.33 the target is reached by removing bridging
   silicates alone and both ratios come back exactly as asked -- at
   ``cs_ratio=1.2, ws_ratio=1.0`` this cell gives Ca/Si = 1.2000 and
   H2O/Si = 1.0000. Beyond that, ``min_mcl`` (3.0 by default) forbids
   further vacancies, so the extra calcium has to go into the interlayer,
   and each of those Ca\ :sup:`2+` brings a water molecule with it.

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

   ['Ca', 'Cw', 'Hw', 'Ob', 'Osi', 'Ow', 'Si']

.. code-block:: text

   <AtomicSystem with 1680 atoms, 420 bonds, 210 angles>

   Box
    a (Å)  b (Å)  c (Å)  α (°)  β (°)  γ (°)
    26.93  30.84  22.68  90.00  90.00  90.00

   Atoms
   type  number      mass  charge
     Ca     123 40.078400     0.0
     Cw     147 40.078400     0.0
     Hw     420  1.007947     0.0
     Ob     120 15.999430     0.0
    Osi     480 15.999430     0.0
     Ow     210 15.999430     0.0
     Si     180 28.085530     0.0

   Bonds
    type  number
   Hw-Ow     420

   Angles
       type  number
   Hw-Ow-Hw     210

   Total charge: 0.000 e
   Volume: 18.83 nm³
   Density: 2.58 g/cm³

The rules separate structural calcium (``Ca``) from interlayer calcium
(``Cw``), and bridging oxygens (``Ob``) from the rest of the silicate
oxygens (``Osi``). Note what the topology contains: **420 bonds, all of
them O-H inside water molecules**. Not a single Si-O bond is listed --
that framework is described by non-bonded interactions in CSHFF.

Step 4: Choose where to cut
============================

This is why the cut plane cannot be picked from the bond list. Cutting at
an arbitrary height looks harmless to the topology and wrecks the
structure anyway:

.. code-block:: python

   naive = Splitter(csh, coordinate=3.0, axis="z", gap_size=25.0)
   print(naive.count_broken_bonds())

.. code-block:: text

   74

Seventy-four severed Si-O contacts -- none of which appear as a bond.
:meth:`~cemd.build.Splitter.scan_broken_bonds` sweeps the whole cell and
reports the damage as a function of the cut coordinate:

.. code-block:: python

   scanner = Splitter(csh, coordinate=0.0, axis="z", gap_size=25.0)
   scan = scanner.scan_broken_bonds(start=0.0, stop=22.68, step=0.5)
   print(scan.to_string(index=False))

.. code-block:: text

    coordinate  n_broken
           0.0         0
           0.5         8
           1.0         8
           ...
           2.5        39
           3.0        74
           3.5        74
           4.0        73
           4.5        34
           5.0         0
           5.5         0
           6.0         0
           6.5         0
           7.0        73
           ...
          14.5        74
          15.0        73
          15.5        73
          16.0         0
          16.5         0
          17.0         0
          17.5         0
          18.0         0
          18.5        73
           ...

The profile maps the material: plateaus at zero are the interlayers, the
peaks around 73-74 are the silicate sheets. Anything in between clips the
chains partially.

.. code-block:: python

   clean = scan[scan["n_broken"] == 0]["coordinate"].tolist()
   print(clean)

.. code-block:: text

   [0.0, 5.0, 5.5, 6.0, 6.5, 16.0, 16.5, 17.0, 17.5, 18.0]

Pick a plane in the middle of one of those windows -- ``z = 17.0`` sits in
the widest one, away from the cell boundary.

Step 5: Define the pore solution
=================================

Cement pore solution is strongly alkaline. ``HO`` is one of the polyatomic
species shipped with CEMD, so it can be named directly -- alongside
``CO3`` and ``SO4``, and any element:

.. code-block:: python

   pore = SolutionBuilder(density=1.0, counts={"Na": 8, "HO": 8})
   print(pore)

.. code-block:: text

   <SolutionBuilder>

   ┌─ Composition
   │   density:  1.00 g/cm³
   │   type:     counts
   │             Na: 8
   │             HO: 8

Only the ions are given explicitly: the builder works out how many water
molecules are needed to reach 1.0 g/cm³ once it knows the volume of the
pore. Any other solute can be supplied as an
:class:`~cemd.core.atomic_system.AtomicSystem` through ``structures=``, as
in :ref:`the caffeine tutorial <tutorial_caffeine_solvation>`.

Step 6: Open the pore and fill it
==================================

:class:`~cemd.build.Splitter` cuts the system, translates the upper part
to open a 25 Å gap, and packs the solution inside it -- ``padding`` keeps
the fluid off the two pore walls:

.. code-block:: python

   splitter = Splitter(csh, coordinate=17.0, axis="z", gap_size=25.0)
   print(splitter.count_broken_bonds())

   result = splitter.add_solution(pore, padding=2.0).split()

.. code-block:: text

   0

Zero broken contacts: the pore opens along an interlayer, so no silicate
chain is cut.

Step 7: Inspect and export
===========================

.. code-block:: python

   print(result.box)
   result.summary()

.. code-block:: text

   [26.93 30.84 47.68 90. 90. 90.]

.. code-block:: text

   <AtomicSystem with 3399 atoms, 1558 bonds, 775 angles>

   Box
    a (Å)  b (Å)  c (Å)  α (°)  β (°)  γ (°)
    26.93  30.84  47.68  90.00  90.00  90.00

   Atoms
   type  number      mass  charge
     Ca     123 40.078400     0.0
     Cw     147 40.078400     0.0
      H       8  1.007947     0.0
     Hw    1550  1.007947     0.0
     Na       8 22.989769     0.0
      O       8 15.999430     0.0
     Ob     120 15.999430     0.0
    Osi     480 15.999430     0.0
     Ow     775 15.999430     0.0
     Si     180 28.085530     0.0

   Bonds
    type  number
     H-O       8
   Hw-Ow    1550

   Angles
       type  number
   Hw-Ow-Hw     775

   Total charge: 0.000 e
   Volume: 39.59 nm³
   Density: 1.67 g/cm³

The cell grew by exactly the gap size (22.68 + 25 = 47.68 Å). The pore
holds the 8 Na\ :sup:`+`, the 8 hydroxides (``H-O``, 8 bonds) and 565 new
water molecules, while the original interlayer water is still in place.

The hydroxide and sodium arrive with the generic types ``H``, ``O`` and
``Na``: run :meth:`~cemd.core.atomic_system.AtomicSystem.set_topology` and
:meth:`~cemd.core.atomic_system.AtomicSystem.set_ff_from_database` again on
the result to type and parameterize them.

.. code-block:: python

   result.view()
   result.write("csh_nanopore.data")

Variant: cutting through a sheet
=================================

Sometimes the cut plane is imposed -- a given pore width, a given
orientation -- and no clean window is available. ``repair=True`` then caps
whatever the cut leaves under-coordinated:

.. code-block:: python

   rough = Splitter(csh, coordinate=12.5, axis="z", gap_size=25.0)
   print(rough.count_broken_bonds())

   repaired = rough.split(repair=True)
   print(rough.repair_report)
   print(repaired)

.. code-block:: text

   19
   {'broken': 19, 'capped': 25, 'skipped': 0}
   <AtomicSystem with 1705 atoms, 404 bonds, 190 angles>

Nineteen severed contacts produced 25 caps: exposed oxygens are
protonated into hydroxyls, exposed cations get a hydroxyl back along the
direction their partner used to occupy. ``skipped`` would count caps
landing on top of an existing atom.

.. note::

   The capping atoms are placed collinearly with the broken contact and
   carry no charge -- this is a starting geometry to relax, not a
   parameterized structure. Re-run
   :meth:`~cemd.core.atomic_system.AtomicSystem.set_topology` and
   :meth:`~cemd.core.atomic_system.AtomicSystem.set_ff_from_database`
   afterwards.
