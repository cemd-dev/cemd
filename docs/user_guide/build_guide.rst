Building Systems
================

The **build** module purpose is to automate the construction of complex atomistic systems, from simple crystal bulks to intricate solid-liquid interfaces.

.. contents:: Table of Contents
   :local:
   :depth: 2

Workflow Overview
-----------------

Building a simulation-ready system generally follows these three steps:

1. **Blueprint**: Define composition and target density with a builder
   (:class:`~cemd.build.SolutionBuilder`, :class:`~cemd.build.GlassBuilder`,
   ...). A blueprint describes *what* to build, not *how*.
2. **Construction**: Call the blueprint's ``build()`` method (or an
   ``AtomicSystem`` method such as
   :meth:`~cemd.core.atomic_system.AtomicSystem.add_liquid_layer`) to
   generate the system, usually via `Packmol <https://m3g.github.io/packmol/>`__.
3. **Refinement**: Assign force-field parameters, apply topology rules, or
   repair broken bonds as needed.

.. note::

   Every builder that calls Packmol requires the ``packmol`` binary to be
   installed and available on your ``$PATH`` (see :doc:`../installation`).

Construction Examples
---------------------

Building a Liquid Solution
^^^^^^^^^^^^^^^^^^^^^^^^^^

:class:`~cemd.build.SolutionBuilder` is a blueprint for a solution: give it
a target density and either explicit molecule/ion counts or target molar
concentrations, and it automatically balances the rest with water.

.. code-block:: python

   from cemd.build import SolutionBuilder

   # By explicit counts
   blueprint = SolutionBuilder(density=1.04, counts={"Na": 16, "Cl": 16})

   # Or by target molar concentration -- converted to counts for the box
   # volume when built
   blueprint = SolutionBuilder(density=1.04, molarities={"Na": 1.0, "Cl": 1.0})

   # Build a standalone box
   system = blueprint.build(box=[30, 30, 30])
   print(system)

.. code-block:: none

   <AtomicSystem with 2693 atoms, 1774 bonds, 887 angles>

   Box
   a (Å)  b (Å)  c (Å)  α (°)  β (°)  γ (°)
      30     30     30     90     90     90

   Atoms
   type  number      mass  charge
     Cl      16 35.453200   0.000
     Hw    1774  1.007947   0.446
     Na      16 22.989769   0.000
     Ow     887 15.999430  -0.892

   Bonds
    type  number
   Hw-Ow    1774

   Angles
       type  number
   Hw-Ow-Hw     887

   Total charge: 0.000 e
   Volume: 27.00 nm³
   Density: 1.04 g/cm³

.. note::

   The water template's own O/H charges come along with it (here from its
   built-in ``.lt`` file). ``Na``/``Cl`` come from a generated single-atom
   template with no assigned charge (0.0) -- assign real force-field
   charges afterward with
   :meth:`~cemd.core.atomic_system.AtomicSystem.set_ff_from_database` or
   :meth:`~cemd.core.atomic_system.AtomicSystem.set_charges` (see the
   :doc:`force-field guide <ff_database>`).

Pure water, or a hemispherical droplet instead of a box, are also
available:

.. code-block:: python

   # Pure water
   water = SolutionBuilder.from_water(density=1.0)
   box_of_water = water.build(box=[20, 20, 20])

   # A standalone hemispherical droplet
   droplet = water.build_hemisphere(radius=15.0, axis="z")

Building a Glass Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^

:class:`~cemd.build.GlassBuilder` builds an amorphous structure from a
composition given as elements, oxides, or a mix of both:

.. code-block:: python

   from cemd.build import GlassBuilder

   # By oxides
   blueprint = GlassBuilder(
       density=2.3,
       composition={"SiO2": 3, "Al2O3": 2, "Na2O": 2},
   )
   glass = blueprint.build(box=[25, 25, 25])

   print(glass.density)  # close to 2.3 g/cm³

   # Equivalently, by elements
   blueprint = GlassBuilder.from_elements(
       density=2.3,
       elements={"Si": 3, "Al": 4, "Na": 4, "O": 14},
   )

If the composition is not charge-neutral, ``GlassBuilder`` raises a
``UserWarning`` when the blueprint is created (it still builds the
structure).

Building C-S-H Models
^^^^^^^^^^^^^^^^^^^^^

For cementitious materials, :class:`~cemd.build.CSHBuilder` generates
C-S-H structures from a target Calcium/Silicate (Ca/Si) and Water/Silicate
(H2O/Si) ratio, starting from a tobermorite crystal:

.. code-block:: python

   from cemd.build import CSHBuilder

   builder = CSHBuilder(cs_ratio=1.5, ws_ratio=1.2)
   csh_system = builder.build(
       supercell=[3, 5, 1],           # replication of the base cell
       model="tob11a_hamid.cif",      # or "tob11a_merlino.cif"
       min_mcl=3.0,                   # minimum mean chain length
   )

   # The builder keeps an analysis of the built structure
   analysis = builder.analyze()
   print(f"Real Ca/Si: {analysis['Ca/Si']:.3f}")
   print(f"MCL: {analysis['MCL']:.2f}")
   print(f"Qⁿ distribution: {analysis['Qn_distribution']}")

``CSHBuilder`` can also wrap an existing C-S-H ``AtomicSystem`` (e.g. one
loaded from a file) to analyze or convert it without rebuilding it:

.. code-block:: python

   system = AtomicSystem.from_file("csh_structure.data")
   builder = CSHBuilder.from_system(system)
   analysis = builder.analyze()

C-A-S-H Models (Aluminum Substitution)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Convert a C-S-H structure to C-A-S-H by substituting a fraction of the
bridging silicon atoms with aluminum:

.. code-block:: python

   cash_system = builder.to_cash(as_ratio=0.10)  # target Al/Si ratio

.. note::

   Substitution targets *bridging* silicon atoms specifically (matching
   real Al-for-Si substitution chemistry), so it requires enough of them
   to be present -- a very short or non-periodic silicate chain may have
   none, in which case the system is returned unchanged.

Building AFt/AFm Structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^

:class:`~cemd.build.AFBuilder` builds AFt (ettringite) or AFm structures
by hydrating a sulfate-containing framework:

.. code-block:: python

   from cemd.build import AFBuilder

   # AFt (ettringite)
   builder = AFBuilder(ws_ratio=10.0)
   aft_system = builder.build_aft(supercell=[2, 2, 1])

   # AFm
   builder = AFBuilder(ws_ratio=8.0)
   afm_system = builder.build_afm(supercell=[3, 1, 1])

Creating Surfaces
-----------------

:class:`~cemd.build.SurfaceBuilder` cuts crystalline surfaces from a bulk
structure along a chosen crystallographic plane and adds a vacuum layer.

.. code-block:: python

   from cemd import AtomicSystem
   from cemd.build import SurfaceBuilder

   # Load a bulk structure
   bulk = AtomicSystem.from_file("tobermorite.cif")

   # Generate surfaces cut along the (1, 0, 4) plane
   builder = SurfaceBuilder(bulk)
   surfaces, shifts, dipoles, broken_bonds = builder.build(
       miller_indices=[1, 0, 4],
       min_slab_size=25.0,   # Å
       min_vacuum_size=15.0,  # Å
   )

   # Select the first generated surface
   my_surface = surfaces[0]

   print(f"Generated {len(surfaces)} surfaces")
   print(f"Broken bonds: {broken_bonds}")

Cut the slab before editing the system. A crystal whose atoms have been
moved, added or removed since it was read, or a slab, is refused with a
``ValueError``: it is no longer the bulk the Miller indices refer to.
Typing (:meth:`~cemd.core.atomic_system.AtomicSystem.set_topology`) moves
nothing and is accepted.

Understanding the Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **Miller Indices**: Define the cutting plane orientation (``[1,0,0]``, ``[1,1,0]``, ``[1,1,1]``, etc.)
- **min_slab_size**: Minimum thickness of the material (too thin = bulk properties not representative)
- **min_vacuum_size**: "Gap" in PBC to prevent surface-surface interactions
- **Broken Bonds**: Number of bonds broken during the cut (diagnostic metric); ``build()`` automatically tries increasing tolerances up to ``max_broken_bonds`` until it finds a valid cut

Creating Solid/Liquid Interfaces
--------------------------------

Once you have a solid surface (and optionally a :class:`~cemd.build.SolutionBuilder`
blueprint), add liquid directly onto it with
:class:`~cemd.core.atomic_system.AtomicSystem` methods:

.. code-block:: python

   from cemd.build import SolutionBuilder

   blueprint = SolutionBuilder(density=1.0, counts={"Na": 10, "Cl": 10})

   # Add a flat liquid layer on top of the surface (in-place, returns self)
   my_surface.add_liquid_layer(
       blueprint,
       thickness=20.0,  # Å
       distance=2.0,    # gap between solid and liquid
       vacuum=10.0,     # vacuum above the liquid
   )

.. code-block:: python

   # Or a hemispherical droplet instead of a flat layer
   droplet_blueprint = SolutionBuilder(density=1.0, counts={"Ca": 5, "Cl": 10})
   my_surface.add_droplet(
       droplet_blueprint,
       radius=15.0,
       distance=2.0,
   )

Adding Custom Structures
------------------------

To place any existing ``AtomicSystem`` on top of a surface (aligned by
center of mass, at a fixed distance), use
:meth:`~cemd.core.atomic_system.AtomicSystem.add_structure`:

.. code-block:: python

   from cemd import AtomicSystem

   custom = AtomicSystem.from_file("custom_molecule.pdb")

   my_surface.add_structure(
       custom,
       distance=3.0,
       axis="z",
       vacuum=15.0,
   )

Splitting Systems
------------------

:class:`~cemd.build.Splitter` cuts a system open along an axis, creating a
gap, and can optionally fill that gap with a solution in the same call
(fluent API):

.. code-block:: python

   from cemd.build import Splitter, SolutionBuilder

   # Simple split: cut at z=15 and open a 20 A gap there
   split_system = Splitter(system, coordinate=15.0, axis="z", gap_size=20.0).split()

   # Split and fill the gap with a solution
   blueprint = SolutionBuilder.from_water(density=1.0)
   filled_system = (
       Splitter(system, coordinate=15.0, axis="z", gap_size=30.0)
       .add_solution(blueprint, padding=2.0)
       .split()
   )

Any bond stretched across the newly opened gap (and any angle/dihedral/
improper built on it) is automatically removed.

Choosing where to cut
^^^^^^^^^^^^^^^^^^^^^

A cut can also sever contacts that carry *no* explicit bond -- the Si-O
framework under ClayFF/CSHFF is treated as non-bonded, so nothing in the
topology marks it, yet slicing through it leaves broken tetrahedra.
:meth:`~cemd.build.Splitter.count_broken_bonds` and
:meth:`~cemd.build.Splitter.scan_broken_bonds` find those contacts
geometrically, from a cutoff table, so you can pick a cut plane that goes
through an interlayer rather than through a silicate sheet:

.. code-block:: python

   splitter = Splitter(csh_system, coordinate=0.0, axis="z", gap_size=20.0)

   scan = splitter.scan_broken_bonds(step=0.5)
   print(scan.head())

.. code-block:: none

      coordinate  n_broken
   0    0.003344         0
   1    0.503344        15
   2    1.003344        19
   3    1.503344        24
   4    2.003344        37

.. code-block:: python

   # Cut where the structure suffers least
   best = scan.loc[scan["n_broken"].idxmin(), "coordinate"]
   result = splitter.with_coordinate(best).split()

The cutoff table (``bonds_dict``) accepts **atom types first, elements as
a fallback** -- the same layered convention as
:class:`~cemd.build.SurfaceBuilder`. Types matter because elements alone
cannot separate structure-bearing oxygens from ones that merely sit
nearby:

.. code-block:: python

   # Default table keys on elements, so water O-H counts as breakable too
   splitter.count_broken_bonds()

   # Framework only: water is free to sit anywhere across the cut
   splitter.count_broken_bonds(bonds_dict={("Si", "O"): 1.8})

Repairing what the cut breaks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When a cut through the framework is unavoidable, ``repair=True`` rebuilds
what it breaks, whether the bond is explicit in the topology or only a
contact of ``bonds_dict`` (the Si-O of ClayFF and CSHFF):

* a hydrogen always goes to the side of its O, so no O-H is cut;
* every other severed bond to an O is hydrolysed. The atom that loses its
  O gets an OH back where the O was, and the O gets an H: a bridging
  Si-O-Si gives two silanols, a non-bridging Si-O...Ca a silanol and a
  Ca-OH, and Si-OH + H2O a silanol and a water;
* a severed bond without an O (C-C, Si-Si...) has no such repair, and
  the split refuses it with a ``ValueError`` before anything moves. Move
  the cut, or use ``keep_molecules=True``.

Set on the splitter, ``repair`` also applies to the counts and the scans,
which then report the bonds the repair will rebuild.

.. code-block:: python

   splitter = Splitter(csh_system, coordinate=12.2, axis="z", gap_size=20.0,
                       repair=True)
   result = splitter.split()

   print(splitter.repair_report)

.. code-block:: none

   {'broken': 319, 'capped': 957, 'crowded': 2, 'retyped': 306, 'untyped': 2,
    'bonds_added': 637, 'angles_added': 29, 'charge_change': -48.27}

The caps go along the broken bond, tilted away when another atom sits too
close; ``crowded`` counts those that found no clear spot. Each atom added,
and each atom whose covalent environment changed (an ``Ob`` that took an
H, the H of an ``Osih`` turned water), takes the type and charge of an
atom of the system in the same environment: ``Osih``/``Hsi`` for a
silanol, ``Oh``/``Hh`` for a Ca-OH, ``Ow``/``Hw`` for a water under CSHFF.
``untyped`` counts the atoms with no such model, typed by their element
(``O``, ``H``) and, if new, uncharged: type them before assigning the
force field. Explicit bonds and angles are added wherever the system
already has them between the same types; Si-O stays implicit.

.. note::

   The charges copied from the system do not add up to zero in general
   (``charge_change``): call
   :meth:`~cemd.core.atomic_system.AtomicSystem.neutralize_charge`
   afterwards, as after
   :meth:`~cemd.core.atomic_system.AtomicSystem.set_ff_from_database`.
