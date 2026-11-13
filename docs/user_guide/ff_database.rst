Force Field Database
====================

.. raw:: html

   <iframe src="../_static/ff_viewer.html" style="width:100%; height:650px; border:none; border-radius:12px;"></iframe>

cemd ships a database of force fields for cementitious and geological
materials, their pore solutions and water. Each model is a TOML file in
``cemd/forcefield/db/`` (see :ref:`ff-models`), loaded together with the
original INTERFACE force field files (IFF v1.5, CHARMM and CVFF forms) and the
GROMOS 54A7 library used for organic molecules.

Every parameter comes with the reference it was taken from, and every file has
been checked against its source paper (tables of the article or of its
supporting information). The database covers:

- **Anhydrous clinker phases**: C\ :sub:`3`\ S, C\ :sub:`3`\ A (``iff``),
  CaO, C\ :sub:`2`\ S, C\ :sub:`3`\ A and CaO-Al\ :sub:`2`\ O\ :sub:`3`-SiO\ :sub:`2`
  glasses (``freitas2015``), silicate melts (``guillot2007``).
- **Hydrates**: C-S-H / tobermorite (``cshff2014``, ``duque2018``, ``iff``),
  portlandite and layered hydroxides (``clayff``, ``clayff_moh``, ``iff``),
  ettringite, monosulfoaluminate (``afff``), gypsum, bassanite and anhydrite
  (``iff``), geopolymers (``sadat2016``).
- **Carbonates**: calcite, aragonite and other alkaline-earth carbonates
  (``raiteri2015``, ``solcon``, ``wang2009``).
- **Clays, silica and oxides**: ``clayff``, ``clayff_moh``, ``iff``.
- **Ions and water**: ClayFF ions, Joung-Cheatham alkali/halide ions for each
  water model (``joung2008``), scaled-charge ions and sulfate (``madrid2019``),
  and the water models SPC, SPC/E, SPC/Fw, TIP3P, TIP4P-Ew and TIP4P/2005.

Supported functional forms: Lennard-Jones 12-6 and Buckingham pair terms;
harmonic, Morse and class2 bonds; harmonic and class2 angles; harmonic and
distance impropers, with the class2 cross terms.

Programmatic Access
-------------------

You can also access the database programmatically using methods available on any :any:`AtomicSystem` object.

Interactive Explorer
^^^^^^^^^^^^^^^^^^^^

The interactive explorer allows you to browse and select force field parameters for each atom type in your system:

.. code-block:: python

   from cemd import AtomicSystem

   # Load your system
   system = AtomicSystem.from_file("my_structure.data")

   # Launch the interactive explorer
   system.explore_ff_database()

.. code-block:: none

   --- Select Forcefield Parameters for type O [117/147] ---
   ↑↓ Move   [Enter] Assign   [p] Open DOI   [q] Quit

     TYPE                 MODEL           | ENVIRONMENT
   --------------------------------------------------------------------------------------
     oc1                  iff             | silicate ion in c3s
     oc2                  iff             | oxide ion in c3s
     oc3                  iff             | aluminate ring in c3a
   ➜ oc4                  iff             | aluminate apical in c3a
     oc5_c3a              iff             | superficial hydroxide ion in c3a
     oc5_c3s              iff             | superficial hydroxide ion in c3s
     oc5a                 iff             | aloh group in c3a and aft
     oc6                  iff             | silanol group in c3s
     o_silica_bulk        iff             | oxygen atom in silica (bulk)
     o_silica_silanol     iff             | oxygen atom in silica (silanol)
     o_sulfate            iff             | oxygen in sulfate ion (calcium sulfates)
     ow_sulfate           iff             | water oxygen in calcium sulfate hydrates (fle
     ox                   iff             | hydroxide oxygen in portlandite
     oc12                 iff             | bridging oxygen in tobermorite (14 A; -0.56 i
     oc13                 iff             | apical silicate oxygen in tobermorite
     oc22                 iff             | apical interlayer silicate oxygen in tobermor
     ocx                  iff             | interlayer bridging oxygen in tobermorite 11
     oc14                 iff             | silanol oxygen in tobermorite
     os                   madrid2019      | oxygen in sulfate
     ow                   madrid2019      | water oxygen (TIP4P/2005, charge on the M sit

   ↕ 114-133 of 147

Manual Assignment
^^^^^^^^^^^^^^^^^

If you know which parameters you want to use, you can assign them directly:

.. code-block:: python

   from cemd import AtomicSystem

   system = AtomicSystem.from_file("my_structure.data")

   # Map your system's atom types to database types. A bare short name
   # (e.g. "ospc") works only if it is unambiguous across every loaded
   # model; several models define the same short name, so prefer the
   # fully-qualified "Model.type" form shown here.
   assignments = {
       'O': 'SPC.ospc',      # Water oxygen (SPC model)
       'H': 'SPC.hspc',      # Water hydrogen (SPC model)
       'Si': 'ClayFF.st',    # Tetrahedral silicon (ClayFF)
       'Ca': 'ClayFF.ca',    # Aqueous calcium (ClayFF)
   }

   # Apply the force field parameters
   system.set_ff_from_database(assignments)

After assignment, the system will have:

- **Masses** and **charges** for each atom type
- **Pair** parameters (Lennard-Jones or Buckingham) for the self and cross
  interactions the model defines
- **Bond**, **angle** and **improper** parameters (if available in the database)

Charges are those of the model. Several models (CSH-FF, IFF tobermorite)
balance their charges structure by structure: check ``system.total_charge``
and use :meth:`~cemd.core.atomic_system.AtomicSystem.neutralize_charge` if the
system is not neutral.

Working with Different Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can mix parameters from different force field models:

.. code-block:: python

   # Mix ClayFF and SPC parameters. A model is named after its file
   # (ForceFieldDatabase().get_model_names() lists them), in any case:
   # "ClayFF.st" and "clayff.st" are the same entry. A type the model does
   # not define is reported and skipped, never looked up in another model.
   assignments = {
       'O': 'clayff.o*',         # Water oxygen from ClayFF
       'H': 'clayff.h*',         # Water hydrogen from ClayFF
       'Si': 'clayff.st',        # Silicon from ClayFF
       'Ca': 'clayff.ca',        # Calcium from ClayFF
       'Ow': 'spc.ospc',         # Water oxygen from SPC
       'Hw': 'spc.hspc',         # Water hydrogen from SPC
   }

   system.set_ff_from_database(assignments)

.. note::
   When mixing models, ensure that the parameters are compatible, especially for cross-interactions between different models.

Checking Applied Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^

After assignment, you can inspect the parameters that were applied:

.. code-block:: python

   # Check masses and charges
   print(f"Masses: {dict(system.masses)}")
   print(f"Charges: {dict(system.charges)}")

   # Force-field parameters live under `ff_params`, keyed by the
   # (canonical, hyphen-joined) interaction type -- e.g. "Ow-Ow" for a
   # pair, "Hw-Ow-Hw" for an angle.
   print(f"Pair parameters: {dict(system.ff_params.pair)}")
   print(f"Bond parameters: {dict(system.ff_params.bond)}")
   print(f"Angle parameters: {dict(system.ff_params.angle)}")

Applying Mixing Rules
^^^^^^^^^^^^^^^^^^^^^

:meth:`~cemd.core.atomic_system.AtomicSystem.set_ff_from_database` already
fills the cross pairs the database does not list, with the mixing rule each
model declares (table below): arithmetic (Lorentz-Berthelot: sigma
arithmetic, epsilon geometric) for most models, geometric for ``wang2009``,
and none for the pair-specific models, whose unlisted pairs stay Coulomb
only. Two types from models with different rules are mixed arithmetically,
with a warning.

To recompute cross pairs, for instance after setting self terms by hand:

.. code-block:: python

   # Rule of each type's model (default)
   system.apply_pair_mixing_rules()

   # Force a rule for every pair
   system.apply_pair_mixing_rules(rule='geometric')

   # Recompute pairs that already have parameters
   system.apply_pair_mixing_rules(overwrite=True)

Mixing only produces Lennard-Jones pairs; a Buckingham pair the model does
not list is left to Coulomb.

.. include:: _generated/ff_mixing.inc

Full Example: Setting Up a Water System
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here is a complete example of setting up a water system with force field parameters:

.. code-block:: python

   from cemd import AtomicSystem

   # Create a water box
   system = AtomicSystem.from_file("waterbox.data")

   # Check current atom types
   print(f"Atom types: {system.atom_types}")
   print(f"Elements: {dict(system.elements)}")

   # Assign SPC water parameters. Prefer the fully-qualified "model.type"
   # key ("SPC.ospc") over the bare short name ("ospc") -- several models
   # define a type with the same short name (e.g. "ospc" also exists under
   # "iff_cvff"), and an unqualified short name resolves to whichever one
   # happens to be loaded first, silently picking the wrong model.
   assignments = {
       'O': 'SPC.ospc',  # Water oxygen
       'H': 'SPC.hspc',  # Water hydrogen
   }

   system.set_ff_from_database(assignments)

   # Verify the assignment
   print(f"Masses: {dict(system.masses)}")
   print(f"Charges: {dict(system.charges)}")
   print(f"LJ parameters: {dict(system.ff_params.pair)}")
   print(f"Bond parameters: {dict(system.ff_params.bond)}")
   print(f"Angle parameters: {dict(system.ff_params.angle)}")

   # Export to LAMMPS data file
   system.write("water_spc.data")

.. code-block:: none

   Atom types: ['H', 'O']
   Elements: {'H': 'H', 'O': 'O'}
   Masses: {'H': 1.007947, 'O': 15.99943}
   Charges: {'H': 0.41, 'O': -0.82}
   LJ parameters: {'O-O': LJParams(epsilon=0.15535, sigma=3.166, ...), 'H-O': LJParams(epsilon=0.0, sigma=1.583, ...), 'H-H': LJParams(epsilon=0, sigma=0, ...)}
   Bond parameters: {'H-O': HarmonicBondParams(k=554.1349, r0=1.0, ...)}
   Angle parameters: {'H-O-H': HarmonicAngleParams(k=45.7696, theta0=109.47, ...)}

Full Example: Setting Up a C-S-H System
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a more complex example with cementitious materials:

.. code-block:: python

   from cemd import AtomicSystem

   # Load a C-S-H structure and type it with the CSH-FF rules:
   # Ca / Cw (intra- / interlayer calcium), Ob / Osi (bridging /
   # non-bridging silicate oxygen), Ow / Hw (water)
   system = AtomicSystem.from_file("csh_structure.data")
   system.set_types_from_elements()
   system.set_topology("cshff")

   # Assign CSH-FF parameters (Shahsavari et al. 2011; LJ of Qomi et al. 2014)
   assignments = {
       'Si': 'cshff2014.si',     # Silicon
       'Osi': 'cshff2014.o',     # Non-bridging silicate oxygen (CSH-FF "O")
       'Ob': 'cshff2014.ob',     # Bridging oxygen
       'Ca': 'cshff2014.ca',     # Intralayer calcium
       'Cw': 'cshff2014.cw',     # Interlayer calcium
       'Ow': 'cshff2014.ow',     # Water oxygen
       'Hw': 'cshff2014.hw',     # Water hydrogen
   }

   system.set_ff_from_database(assignments)

   # CSH-FF charges are balanced model by model: spread the excess over the
   # silicate oxygens
   print(f"Total charge: {system.total_charge:.3f} e")
   system.neutralize_charge({'Osi': 1.0, 'Ob': 1.0})

   # Export
   system.write("csh_cshff2014.data")

.. _ff-models:

Available Models
----------------

The database includes the following models, grouped by the kind of system
they describe. The key is the name used in ``model.type`` (``clayff.st``).

.. include:: _generated/ff_models.inc

Notes on specific models
^^^^^^^^^^^^^^^^^^^^^^^^

- ``clayff`` follows Cygan et al. (2004); ``clayff_moh`` follows the 2021
  tables of Cygan et al., which differ for epsilon of ``feo``/``lio`` and
  sigma of ``cah``. In ``clayff_moh`` the O-H bond is a Morse potential
  (dioctahedral set on ``oh``, trioctahedral set on ``ohs``) and the M-O-H
  angles use the bulk values (110 degrees for surface hydroxyls of Al and Mg).
- ``cshff2014`` carries the base CSH-FF charges; Qomi et al. (2014) adjust the
  oxygen charges of each model for neutrality (see the C-S-H example above).
- ``iff`` holds the 12-6 parameters of several IFF papers: C\ :sub:`3`\ S
  (Mishra 2013, Pustovgar 2017), C\ :sub:`3`\ A (Mishra 2014), silica
  (Emami 2014), calcium sulfates (Mishra 2021), portlandite and tobermorite
  14 Å (Mishra 2017). The tobermorite 11 Å differs only by a few charges,
  given in the file.
- ``solcon``: the carbonate intramolecular terms are not published by
  Armstrong et al. (2023); the harmonic terms of ``raiteri2015`` are used.
- ``joung2008`` has one type per ion and water model (``na_spce``,
  ``na_tip3p``, ``na_tip4pew``...): pick the one matching your water.
- ``madrid2019``, ``tip4p2005`` and ``tip4p_ew`` are rigid four-site models:
  the oxygen carries the M-site charge, bonds and angles have k = 0 and only
  fix the geometry. Use LAMMPS ``lj/cut/tip4p/long`` with d\ :sub:`OM` =
  0.1546 Å (TIP4P/2005) or 0.125 Å (TIP4P-Ew) and SHAKE.
- Not included, because cemd has no core-shell model or n-m potentials yet:
  Cement FF (Galmarini 2013) and ERICA FF (Valavi et al. 2022), as well as the
  reactive ReaxFF parameter sets.

Writing the LAMMPS input
^^^^^^^^^^^^^^^^^^^^^^^^

:meth:`~cemd.core.atomic_system.AtomicSystem.write` produces the LAMMPS data
file with its coefficient sections; the styles go in your input script. When a
bonded section mixes styles, as harmonic water bonds with the Morse hydroxyls
of ``clayff_moh``, the data file uses LAMMPS's hybrid form:

.. code-block:: none

   Bond Coeffs  # hybrid

   1   harmonic     554.1349     1.0000   # h*-o*
   2   morse        132.2491     2.1815     0.9450   # ho-oh

and the input script must declare ``bond_style hybrid harmonic morse``. A
pair section mixing Lennard-Jones and Buckingham terms (``raiteri2015``,
``solcon``, ``sadat2016``) likewise needs ``pair_style hybrid/overlay`` with
``lj/cut/coul/long`` and ``buck/coul/long``.

Energies are in kcal/mol and distances in Å (LAMMPS ``units real``).

.. _ff-add-model:

Adding a model
--------------

A model is one TOML file in ``cemd/forcefield/db/``; its name (the file stem)
is the prefix of every type (``clayff.st``). Name it after the force field
when it has one (``clayff``, ``solcon``, ``tip4p2005``), otherwise after the
first author and the year in lowercase (``raiteri2015``, ``freitas2015``). A
combination of existing force fields is named after its author and year
(``duque2018``).

.. code-block:: toml

   [model]
   name = "ClayFF"
   ref = "https://doi.org/10.1021/jp0363287"
   description = "Force field for clay minerals and aqueous solutions"
   category = "minerals"
   mixing = "arithmetic"

   [atom.st]
   element = "Si"
   charge = 2.1
   environment = "tetrahedral silicon"
   ref = "https://doi.org/10.1021/jp0363287"

   [lj."st-st"]               # E = 4 eps [(s/r)^12 - (s/r)^6]
   epsilon = 1.8405e-06       # kcal/mol
   sigma = 3.30202700891535   # A

   [buckingham."Ca-Oc"]       # E = A exp(-r/rho) - C/r^6
   A = 72914.9
   rho = 0.271511
   C = 0.0

   [bond.harmonic."ho-oh"]    # E = K (r - r0)^2
   k = 554.1349
   r0 = 1.0

   [bond.morse."ho-oh"]       # E = D [1 - exp(-alpha (r - r0))]^2
   D = 132.2491
   alpha = 2.1815
   r0 = 0.945

   [angle.harmonic."h*-o*-h*"]  # E = K (theta - theta0)^2
   k = 45.7696
   theta0 = 109.47

Other sections: ``[bond.class2]``, ``[angle.class2]``,
``[improper.harmonic]``, ``[improper.distance]``, ``[bondbond.class2]`` and
``[bondangle.class2]``. Keys join the short type names with ``-``; both
orders of a bond or angle key are looked up.

Conversions to keep in mind when copying a table:

- **Harmonic constants**: cemd and LAMMPS use E = K (r - r0)\ :sup:`2`. A
  paper that writes E = (k/2)(r - r0)\ :sup:`2` (SPC/Fw, GROMACS) needs
  K = k/2.
- **LJ distance**: INTERFACE and ClayFF tabulate the distance of the minimum
  (sigma0 or R0); sigma = R0 / 2\ :sup:`1/6`.
- **Units**: 1 eV = 23.060548 kcal/mol, 1 kcal = 4.184 kJ.

Besides ``name``, ``ref`` and ``description``, the ``[model]`` section
declares:

.. code-block:: toml

   category = "minerals"      # minerals, carbonates, ions or water
   mixing = "arithmetic"      # arithmetic, geometric or none
   mixing_note = "..."        # optional, shown in the mixing-rule table

The test suite checks that every type named in a pair, bond or angle key is
declared under ``[atom]`` of the same file, and that ``category`` and
``mixing`` take one of the values above, so a misspelt type fails the tests
instead of silently dropping a term. The viewer, the table of
:ref:`ff-models` and the mixing-rule table are generated from the TOML files
when the documentation is built: a new file appears on this page without
editing it.
