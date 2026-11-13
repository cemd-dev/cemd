Topology Module
===============

The topology module types the atoms of a system and creates its bonds, angles
and dihedrals from neighbour rules, applied with
:meth:`~cemd.core.atomic_system.AtomicSystem.set_topology`.

Rules
-----

.. autosummary::
   :toctree: generated/

   cemd.topology.TopologyRule
   cemd.topology.NeighborCriterion
   cemd.topology.LayerRule
   cemd.topology.DihedralRule

Predefined rule sets
--------------------

``set_topology("clayff")`` and ``set_topology("cshff")`` apply these lists.

.. autodata:: cemd.topology.CLAYFF_RULES
   :no-value:

.. autodata:: cemd.topology.CSHFF_RULES
   :no-value:
