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

``set_topology("clayff")``, ``set_topology("cshff")`` and
``set_topology("iff_silica")`` apply these lists.

.. autodata:: cemd.topology.CLAYFF_RULES
   :no-value:

.. autodata:: cemd.topology.CSHFF_RULES
   :no-value:

.. autodata:: cemd.topology.IFF_SILICA_RULES
   :no-value:
