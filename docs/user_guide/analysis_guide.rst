Analyzing MD Trajectories
=========================

The **analysis** module provides tools to extract physical properties from molecular dynamics simulations. It is designed to work seamlessly with `MDAnalysis` universes, allowing you to compute structural, dynamical, and electrostatic properties.

.. contents:: Table of Contents
   :local:
   :depth: 2

Loading Trajectories
--------------------

Before analysis, load your trajectory into an MDAnalysis Universe:

.. code-block:: python

   import MDAnalysis as mda
   
   # From LAMMPS data + trajectory
   universe = mda.Universe("system.data", "production.dcd")
   
   # From AtomicSystem + trajectory
   from cemd import AtomicSystem
   system = AtomicSystem.from_file("system.data")
   universe = system.to_mda()
   # Add trajectory
   mda.Universe(universe, "production.dcd")
   
   # From other formats
   universe = mda.Universe("topology.psf", "trajectory.dcd")
   universe = mda.Universe("structure.gro", "trajectory.xtc")

Structural Analysis
-------------------

Radial Distribution Function (RDF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The RDF, :math:`g(r)`, describes the probability of finding a particle at a distance :math:`r` from another.

.. code-block:: python

   from cemd.analysis import compute_rdf

   # Basic RDF between Si and O
   rdf_df = compute_rdf(
       universe,
       type1="Si",
       type2="O",
       cutoff=10.0,  # Maximum distance (Å)
       dr=0.1        # Bin width (Å)
   )
   
   # RDF between water oxygen atoms
   rdf_oo = compute_rdf(universe, type1="Ow", type2="Ow", cutoff=8.0)
   
   # RDF with all atoms (type1="all")
   rdf_all = compute_rdf(universe, type1="all", type2="Ca", cutoff=12.0)

The resulting DataFrame contains:
- `g_r`: Radial distribution function
- `G_r`: Total correlation function :math:`G(r) = 4πrρ(g(r)-1)`
- `n_r`: Running coordination number :math:`n(r)`

.. code-block:: python

   # Access the data
   print(rdf_df.head())
   print(f"First coordination shell at: {rdf_df.index[rdf_df['g_r'].argmax()]:.2f} Å")

Silicate Network Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^

For silicate materials, you can analyze the polymerization state:

.. code-block:: python

   from cemd.analysis.silicates import analyze_silicates

   # Analyze the silicate network
   results = analyze_silicates(
       universe,
       types_map={
           "si_types": "Si",
           "o_types": "O Ob Osi",
           "al_types": "Al",
           "ca_types": "Ca"
       },
       cutoff=1.85  # Si-O bond cutoff
   )
   
   print(f"Ca/(Si+Al): {results['Ca/(Si+Al)']:.3f}")
   print(f"Al/Si: {results['Al/Si']:.3f}")
   print(f"MCL (Mean Chain Length): {results['MCL']:.2f}")
   print(f"Q^n distribution: {results['Qn_distribution']}")

Dynamical Analysis
------------------

Mean Squared Displacement (MSD)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The MSD is used to calculate the diffusion coefficient :math:`D` via Einstein's relation:

.. math::
   D = \\frac{1}{6} \\lim_{\\tau \\to \\infty} \\frac{\\langle [\\mathbf{r}(t_0 + \\tau) - \\mathbf{r}(t_0)]^2 \\rangle}{\\tau}

.. code-block:: python

   from cemd.analysis import msd, diffusion_coefficient, msd_profile

   # Bulk MSD for water oxygen
   msd_data = msd(
       universe,
       atom_type="Ow",
       dt=100.0,  # Timestep in fs
       nblocks=None,  # Auto-calculate
       corrlength=500,
       gaplength=50
   )
   
   # Extract diffusion coefficient
   dc_df = diffusion_coefficient(msd_data, start=1.0)  # Start fit from 1 ps
   print(dc_df)

   # Plot MSD
   from cemd.analysis.diffusion import plot_msd
   plot_msd(msd_data)

Spatial MSD Profile
^^^^^^^^^^^^^^^^^^^

For interfacial systems, compute diffusion as a function of distance from the surface:

.. code-block:: python

   # Compute MSD profile along z-axis
   msd_profile_data, std_profile, _ = msd_profile(
       universe,
       atom_type="Ow",
       dt=100.0,
       axis='z',
       nblocks=None,
       corrlength=100,
       gaplength=50,
       delta=1.0  # Spatial bin size (Å)
   )
   
   # Calculate diffusion coefficient profile
   from cemd.analysis.diffusion import diffusion_coefficient_profile
   dc_profile = diffusion_coefficient_profile(
       msd_profile_data,
       std_profile,
       start=1.0
   )
   
   # Plot the profile
   from cemd.analysis.diffusion import plot_diffusion_profile
   plot_diffusion_profile(dc_profile['iso'])  # Isotropic component

Bond Correlation and Lifetime
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Analyze bond survival probability:

.. code-block:: python

   from cemd.analysis.tcf import bondcorr, lifetime

   # Compute bond correlation function
   tcf_data = bondcorr(
       universe,
       atom_types_a="Si",
       atom_types_b="O",
       distance=2.0,  # Bond cutoff (Å)
       dt=100.0,      # Timestep (fs)
       corrlength=200
   )
   
   # Calculate bond lifetime
   tau, params = lifetime(tcf_data, corrtime=0.1)
   print(f"Bond lifetime: {tau:.3f} ps")

Density Analysis
----------------

1D Density Profiles
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from cemd.analysis import density_profile

   # Density profile of water and calcium along z-axis
   profile = density_profile(
       universe,
       atom_types=["Ow", "Ca"],
       axis="z",
       start=0,
       end=-1,  # All frames
       bin_size=0.1
   )
   
   # Plot the profile
   import matplotlib.pyplot as plt
   profile.plot()
   plt.xlabel("z (Å)")
   plt.ylabel("Density (kg/m³)")
   plt.show()

2D Density Maps
^^^^^^^^^^^^^^^

.. code-block:: python

   from cemd.analysis import density_map

   # 2D density map near an interface
   density_2d = density_map(
       universe,
       atom_types=["Ow", "Ca"],
       interface_coordinate=10.0,  # Interface position along z
       axis="z",
       eps=3.0,  # Distance from interface
       bin_size=0.5
   )
   
   # Plot as heatmap
   import matplotlib.pyplot as plt
   plt.imshow(density_2d.T, origin='lower', aspect='auto')
   plt.colorbar(label="Density (kg/m³)")
   plt.show()

Interface Detection
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from cemd.analysis.density import find_interfaces_coordinates, shift_profile

   # Find solid/liquid interface positions
   solid_types = ["Si", "Ca"]
   liquid_types = ["Ow", "Hw"]
   
   left_liquid, left_solid, right_solid, right_liquid = find_interfaces_coordinates(
       profile,
       solid_types,
       liquid_types
   )
   
   print(f"Left interface: {left_liquid:.2f} - {left_solid:.2f} Å")
   print(f"Right interface: {right_solid:.2f} - {right_liquid:.2f} Å")
   
   # Shift the profile
   shifted = shift_profile(profile, shift=2.0)

Electrostatic Potential
^^^^^^^^^^^^^^^^^^^^^^^

Compute the charge density, electric field, and electrostatic potential
from a set of per-type density profiles and their charges:

.. code-block:: python

   from cemd.analysis import density_profile
   from cemd.analysis.density import electrostatic_potential

   # One profile with a column per atom type
   profile = density_profile(
       universe,
       atom_types=["Ow", "Hw", "Ca"],
       axis="z",
       bin_size=0.1,
   )

   # Charges must be given in the same order as the profile's columns
   charges = [-0.82, 0.41, 2.0]

   charge_density, efield, potential = electrostatic_potential(profile, charges)

   import matplotlib.pyplot as plt
   potential.plot()
   plt.xlabel("z (Å)")
   plt.ylabel("Potential (V)")
   plt.show()