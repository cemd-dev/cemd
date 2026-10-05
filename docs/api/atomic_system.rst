.. currentmodule:: cemd.core.atomic_system

AtomicSystem
============

.. autoclass:: AtomicSystem
   :members: atoms, bonds, angles, dihedrals, impropers, velocities

Properties
----------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   AtomicSystem.box
   AtomicSystem.volume
   AtomicSystem.density
   AtomicSystem.total_charge
   AtomicSystem.total_mass
   AtomicSystem.masses
   AtomicSystem.charges
   AtomicSystem.elements
   AtomicSystem.ff_keys
   AtomicSystem.ff_params
   AtomicSystem.num_atoms
   AtomicSystem.num_bonds
   AtomicSystem.num_angles
   AtomicSystem.num_dihedrals
   AtomicSystem.num_impropers
   AtomicSystem.num_atom_types
   AtomicSystem.num_bond_types
   AtomicSystem.num_angle_types
   AtomicSystem.num_dihedral_types
   AtomicSystem.num_improper_types
   AtomicSystem.atom_types
   AtomicSystem.bond_types
   AtomicSystem.angle_types
   AtomicSystem.dihedral_types
   AtomicSystem.improper_types

Inspection
----------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   AtomicSystem.summary
   AtomicSystem.get_count
   AtomicSystem.get_center_of_mass
   AtomicSystem.copy
   AtomicSystem.view

Input/Output
------------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   AtomicSystem.from_file
   AtomicSystem.from_smiles
   AtomicSystem.from_cod
   AtomicSystem.from_pubchem
   AtomicSystem.from_mda
   AtomicSystem.from_pmg
   AtomicSystem.to_mda
   AtomicSystem.to_pmg
   AtomicSystem.write

Topology
--------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   AtomicSystem.set_types
   AtomicSystem.set_type2atoms
   AtomicSystem.set_types_from_elements
   AtomicSystem.set_topology
   AtomicSystem.remove_all_connections
   AtomicSystem.remove_connection_types
   AtomicSystem.keep_connection_types
   AtomicSystem.guess_connections
   AtomicSystem.guess_angles
   AtomicSystem.guess_dihedrals
   AtomicSystem.guess_impropers
   AtomicSystem.add_bond
   AtomicSystem.add_angle
   AtomicSystem.add_dihedral
   AtomicSystem.add_improper

Force Field
-----------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   AtomicSystem.set_masses
   AtomicSystem.set_charges
   AtomicSystem.neutralize_charge
   AtomicSystem.set_ff_from_database
   AtomicSystem.explore_ff_database
   AtomicSystem.set_pair_params
   AtomicSystem.set_bond_params
   AtomicSystem.set_angle_params
   AtomicSystem.set_dihedral_params
   AtomicSystem.set_improper_params
   AtomicSystem.apply_pair_mixing_rules
   AtomicSystem.set_ff_keys
   AtomicSystem.set_atom_ff_keys
   AtomicSystem.set_bond_ff_keys
   AtomicSystem.set_angle_ff_keys
   AtomicSystem.set_dihedral_ff_keys
   AtomicSystem.set_improper_ff_keys

Modifiers
---------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   AtomicSystem.add_structure
   AtomicSystem.add_liquid_layer
   AtomicSystem.add_droplet
   AtomicSystem.add_atom
   AtomicSystem.add_atoms
   AtomicSystem.protonate_atom
   AtomicSystem.protonate_atoms
   AtomicSystem.remove_atom
   AtomicSystem.remove_atoms
   AtomicSystem.set_atom_position
   AtomicSystem.set_box
   AtomicSystem.orthogonalize
   AtomicSystem.unskew
   AtomicSystem.replicate
   AtomicSystem.wrap
   AtomicSystem.center_on_com

