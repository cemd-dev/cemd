Analysis Module
===============

The analysis package provides tools for computing structural,
dynamical, and interfacial properties from molecular dynamics
simulations.

Distribution Functions
----------------------

.. autosummary::
   :toctree: generated/

   cemd.analysis.compute_rdf

Silicate Network
----------------

.. autosummary::
   :toctree: generated/

   cemd.analysis.analyze_silicates
   cemd.analysis.analyze_coordination

Diffusion
---------

.. autosummary::
   :toctree: generated/

   cemd.analysis.msd
   cemd.analysis.msd_profile
   cemd.analysis.diffusion_coefficient
   cemd.analysis.diffusion.diffusion_coefficient_profile
   cemd.analysis.diffusion.plot_msd
   cemd.analysis.diffusion.plot_diffusion_profile

Bond Correlation
----------------

.. autosummary::
   :toctree: generated/

   cemd.analysis.tcf.bondcorr
   cemd.analysis.tcf.lifetime

Interface Properties
--------------------

.. autosummary::
   :toctree: generated/

   cemd.analysis.density_profile
   cemd.analysis.density_map
   cemd.analysis.electrostatic_potential
   cemd.analysis.density.find_interfaces_coordinates
   cemd.analysis.density.shift_profile