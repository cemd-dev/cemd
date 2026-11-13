Force Field Module
==================

The force-field database and the parameter classes it returns. The models
themselves are described in :doc:`../user_guide/ff_database`.

Database
--------

.. autosummary::
   :toctree: generated/

   cemd.forcefield.ForceFieldDatabase
   cemd.forcefield.ForceFieldModel
   cemd.forcefield.AtomType

Pair parameters
---------------

.. autosummary::
   :toctree: generated/

   cemd.forcefield.LJParams
   cemd.forcefield.BuckinghamParams

Bonded parameters
-----------------

.. autosummary::
   :toctree: generated/

   cemd.forcefield.HarmonicBondParams
   cemd.forcefield.MorseBondParams
   cemd.forcefield.Class2BondParams
   cemd.forcefield.HarmonicAngleParams
   cemd.forcefield.Class2AngleParams
   cemd.forcefield.Class2BondBondParams
   cemd.forcefield.Class2BondAngleParams
   cemd.forcefield.HarmonicImproperParams
   cemd.forcefield.DistanceImproperParams
