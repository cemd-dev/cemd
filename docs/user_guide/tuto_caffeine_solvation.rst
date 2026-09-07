.. _tutorial_caffeine_solvation:

============================================================================
Tutorial: Solvating Caffeine in Saline Water, Force Field Included
============================================================================

This tutorial solvates 3 caffeine molecules in a sodium chloride solution.
The molecule is read from a **moltemplate** file produced by the
`Automated Topology Builder <https://atb.uq.edu.au>`_, so it arrives with
its topology *and* its force-field assignment already attached: bonds,
angles, dihedrals, impropers, GROMOS atom-type keys and ATB partial
charges. Nothing has to be guessed afterwards.

.. contents::
   :local:
   :depth: 2

Prerequisites
==============

.. code-block:: python

   from cemd import AtomicSystem
   from cemd.build import SolutionBuilder

Step 1: Read the molecule and its topology
===========================================

.. code-block:: python

   caffeine = AtomicSystem.from_file("caffeine.lt")
   print(caffeine)
   print(caffeine.atom_types)

.. code-block:: text

   <AtomicSystem with 24 atoms, 25 bonds, 43 angles, 54 dihedrals, 8 impropers>
   ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'H1', 'H10', 'H2', 'H3',
    'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'N1', 'N2', 'N3', 'N4', 'O1', 'O2']

The connectivity is complete straight away -- compare this with a
structure fetched from PubChem, which carries coordinates and bonds but no
angles or dihedrals and would need
:meth:`~cemd.core.atomic_system.AtomicSystem.guess_connections`.

Each atom is its own type (``C1``, ``C2``, ... rather than a single ``C``),
because the ATB assigns a distinct partial charge to every position:

.. code-block:: python

   print(f"{caffeine.total_charge:+.4f} e")
   print(caffeine.charges["C1"], caffeine.charges["N1"], caffeine.charges["O1"])

.. code-block:: text

   +0.0000 e
   0.536 -0.358 -0.583

Step 2: Inspect what the file already provides
===============================================

The moltemplate file also names, for every atom and every interaction, the
force-field entry it should use -- the *ff keys*. It does not carry the
parameters themselves:

.. code-block:: python

   print(len(caffeine.ff_keys.atom), len(caffeine.ff_keys.bond),
         len(caffeine.ff_keys.angle), len(caffeine.ff_keys.dihedral))
   print(list(caffeine.ff_keys.atom.items())[:3])

   print(len(caffeine.ff_params.bond), len(caffeine.ff_params.angle))

.. code-block:: text

   24 25 43 54
   [('C1', 'gromos.CAro'), ('C2', 'gromos.COpt'), ('C3', 'gromos.CAro')]
   0 0

Keys but no parameters: that gap is what
:meth:`~cemd.core.atomic_system.AtomicSystem.set_ff_from_database` fills in
at step 5, once the whole system is assembled.

Step 3: Describe the solution
==============================

Ask for 3 caffeine molecules and 10 Na\ :sup:`+` / 10 Cl\ :sup:`-` at a
target density of 1.0 g/cm³. Species that are not plain elements are
supplied through ``structures``:

.. code-block:: python

   blueprint = SolutionBuilder(
       density=1.0,
       counts={"caffeine": 3, "Na": 10, "Cl": 10},
       structures={"caffeine": caffeine},
   )
   print(blueprint)

.. code-block:: text

   <SolutionBuilder>

   ┌─ Composition
   │   density:  1.00 g/cm³
   │   type:     counts
   │             caffeine: 3
   │             Na: 10
   │             Cl: 10
   │
   ├─ Custom structures
   │   caffeine: 24 atoms

Step 4: Build the box
======================

The builder works out how much water is needed to reach the target density
once the solutes are accounted for, then packs everything with Packmol:

.. code-block:: python

   system = blueprint.build(box=[30.0, 30.0, 30.0])
   print(system)
   print(system.get_count("Ow"), "water molecules")

.. code-block:: text

   <AtomicSystem with 2606 atoms, 1751 bonds, 967 angles, 162 dihedrals, 24 impropers>
   838 water molecules

The caffeine topology was replicated for each copy: 3 × 54 = 162
dihedrals and 3 × 8 = 24 impropers survive the packing, and so do the ATB
charges.

.. code-block:: python

   print(system.charges["C1"], system.charges["Ow"], system.charges["Hw"])

.. code-block:: text

   0.536 -0.892 0.446

Step 5: Give the ions their charge
===================================

Monoatomic species are generated on the fly by the builder, so they come
out **neutral** -- the builder has no way to know which oxidation state you
mean. Set them explicitly:

.. code-block:: python

   system.set_charges({"Na": 1.0, "Cl": -1.0})
   print(f"{system.total_charge:+.4f} e")

.. code-block:: text

   +0.0000 e

.. warning::

   Check this before running anything. A silently neutral Na\ :sup:`+` is
   easy to miss and produces a physically meaningless trajectory.

Step 6: Resolve the force-field parameters
===========================================

With the ff keys already in place, the database lookup needs no arguments:

.. code-block:: python

   system.set_ff_from_database()

   print(len(system.ff_params.bond), len(system.ff_params.angle),
         len(system.ff_params.dihedral), len(system.ff_params.improper))
   print(list(system.ff_params.bond.items())[:1])

.. code-block:: text

   26 43 42 8
   [('C1-C3', HarmonicBondParams(k=500.678776290631, r0=1.38, ref='', model='gromos'))]

Every bond, angle, dihedral and improper type in the box is parameterized,
and the partial charges are untouched -- GROMOS defines its charges per
atom in the topology rather than per atom type, so the lookup leaves them
alone.

Step 7: Export
===============

.. code-block:: python

   system.view()
   system.write("caffeine_nacl.data")

.. image:: /_static/images/solvated_caffeine.png
   :alt: Caffeine solvated in an NaCl solution
   :align: center
   :width: 400px

.. note::

   **Starting from PubChem instead.**
   :meth:`~cemd.core.atomic_system.AtomicSystem.from_pubchem` fetches a 3D
   structure by name, which is the quickest route when you only need
   geometry. It returns atoms and bonds only, so you then have to call
   :meth:`~cemd.core.atomic_system.AtomicSystem.guess_connections` for the
   angles and dihedrals and assign the force-field keys yourself. The
   moltemplate route used above skips both steps.
