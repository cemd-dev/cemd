Force Field Database
====================

.. raw:: html

   <iframe src="../_static/ff_viewer.html" style="width:100%; height:650px; border:none; border-radius:12px;"></iframe>

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

   --- Select Forcefield Parameters for type O [1/32] ---
   ↑↓ Move   [Enter] Assign   [p] Open DOI   [q] Quit

   TYPE                 MODEL           | ENVIRONMENT
   --------------------------------------------------------------------------------------
   osm                  guillot2007     | oxygen in silica melt
   oc1                  iff             | silicate ion in c3s
   oc2                  iff             | oxide ion in c3s
   oc3                  iff             | aluminate ring in c3a
   oc4                  iff             | aluminate apical in c3a
   oc5_c3a              iff             | superficial hydroxide ion in c3a
   oc5_c3s              iff             | superficial hydroxide ion in c3s
   oc5a                 iff             | aloh group in c3a and aft
   oc6                  iff             | silanol group in c3s
   o_silica_bulk        iff             | oxygen atom in silica (bulk)
   o_silica_silanol     iff             | oxygen atom in silica (silanol)
   ogp                  sperinck        | geopolymer/metakaolin
   owgp                 sperinck        | fw/spc water for geopolymer/metakaolin
   o_star               clayff          | water oxygen
 ➜ ob                   clayff          | bridging oxygen in csh
   obos                 clayff          | bridging oxygen with octahedral substitution
   obss                 clayff          | bridging oxygen with double substitution
   obts                 clayff          | bridging oxygen with tetrahedral substitution
   oh                   clayff          | hydroxyl oxygen
   ohs                  clayff          | hydroxyl oxygen with substitution

   ↕ 1-20 of 32

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
- **Lennard-Jones 12-6** parameters for all pairs (self and cross-interactions)
- **Bond** and **angle** parameters (if available in the database)

Working with Different Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can mix parameters from different force field models:

.. code-block:: python

   # Mix ClayFF and SPC parameters. Model names are case-sensitive as
   # loaded (check with ForceFieldDatabase().get_model_names()) --
   # "clayff.o_star" would *not* raise an error, it would silently fall
   # back to an arbitrary model that also defines an "o_star" type.
   assignments = {
       'O': 'ClayFF.o_star',     # Water oxygen from ClayFF
       'H': 'ClayFF.h*',     # Water hydrogen from ClayFF
       'Si': 'ClayFF.st',        # Silicon from ClayFF
       'Ca': 'ClayFF.ca',        # Calcium from ClayFF
       'Ow': 'SPC.ospc',         # Water oxygen from SPC
       'Hw': 'SPC.hspc',         # Water hydrogen from SPC
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

If cross-interactions are not explicitly defined in the database, you can apply mixing rules:

.. code-block:: python

   # Apply arithmetic mixing rules
   system.apply_pair_mixing_rules(rule='arithmetic')

   # Or geometric mixing rules
   system.apply_pair_mixing_rules(rule='geometric')

   # Overwrite existing cross-interactions
   system.apply_pair_mixing_rules(rule='arithmetic', overwrite=True)

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

   # Load a C-S-H structure
   system = AtomicSystem.from_file("csh_structure.data")

   # Assign CSHFF2014 parameters
   assignments = {
       'Si': 'CSHFF2014.si_cshff2014',    # Silicon in C-S-H
       'O': 'CSHFF2014.o_cshff2014',      # Oxygen in C-S-H
       'Ca': 'CSHFF2014.ca_cshff2014',    # Calcium in C-S-H
       'H': 'CSHFF2014.h_cshff2014',      # Hydrogen in C-S-H
       'Ow': 'CSHFF2014.ow_cshff2014',    # Water oxygen in C-S-H
       'Hw': 'CSHFF2014.hw_cshff2014',    # Water hydrogen in C-S-H
   }

   system.set_ff_from_database(assignments)

   # Check the results
   print(f"Total charge: {system.total_charge:.3f} e")
   print(f"Number of atom types: {system.num_atom_types}")

   # Export
   system.write("csh_cshff2014.data")

Available Models
----------------

The database currently includes the following force field models:

.. list-table::
   :header-rows: 1
   :align: center
   :widths: 20 20 60

   * - Model
     - Reference
     - Description
   * - ClayFF
     - `10.1021/jp0363287 <https://doi.org/10.1021/jp0363287>`_
     - Force field for clay minerals and aqueous solutions
   * - CSHFF2014
     - `10.1038/ncomms5960 <https://doi.org/10.1038/ncomms5960>`_
     - Force field for calcium-silicate-hydrate
   * - IFF
     - `10.1039/C4DT00438H <https://doi.org/10.1039/C4DT00438H>`_
     - Interface Force Field for cement minerals
   * - SPC
     - `10.1021/j100308a038 <https://doi.org/10.1021/j100308a038>`_
     - Simple Point Charge water model
   * - TIP3P
     - `10.1063/1.1808117 <https://doi.org/10.1063/1.1808117>`_
     - Transferable Intermolecular Potential 3-point water model
   * - TIP4P-EW
     - `10.1063/1.1683075 <https://doi.org/10.1063/1.1683075>`_
     - Transferable Intermolecular Potential 4-point water model (Ewald)
   * - TIP4P-2005
     - `10.1063/1.2121687 <https://doi.org/10.1063/1.2121687>`_
     - Transferable Intermolecular Potential 4-point water model (2005)
   * - Guillot & Sator 2007
     - `10.1016/j.gca.2006.11.015 <https://doi.org/10.1016/j.gca.2006.11.015>`_
     - Force field for silica melts
   * - Sperinck 2016
     - `10.1063/1.4964301 <https://doi.org/10.1063/1.4964301>`_
     - Force field for geopolymer/metakaolin
   * - SolCon
     - `10.1098/rsta.2022.0250 <https://doi.org/10.1098/rsta.2022.0250>`_
     - Force field for calcium carbonate (solubility-consistent)
   * - Wang 2008
     - `10.2138/am.2009.2939 <https://doi.org/10.2138/am.2009.2939>`_
     - Force field for vaterite