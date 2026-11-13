Kaolinite with ClayFF-MOH
=========================

Kaolinite, Al\ :sub:`2`\ Si\ :sub:`2`\ O\ :sub:`5`\ (OH)\ :sub:`4`, is the main
mineral of the calcined clays used in LC\ :sup:`3` cements. Its hydroxyls carry
most of what changes on heating and grinding (dehydroxylation, disorder of the
octahedral sheet), so their vibrations and orientations matter. ClayFF-MOH
(Cygan et al. 2021, Pouvreau et al. 2017) describes them with a Morse O-H
stretch and metal-O-H angle terms.

This tutorial reads kaolinite from the Crystallography Open Database, types it,
adds the Al-O-H angles with a custom topology rule, substitutes a little Fe
for Al with :meth:`~cemd.core.atomic_system.AtomicSystem.substitute_atoms`,
and assigns ClayFF-MOH.

Prerequisites
-------------

.. code-block:: python

   from cemd import AtomicSystem
   from cemd.analysis import analyze_coordination
   from cemd.topology import NeighborCriterion, TopologyRule

Step 1: Read the structure
--------------------------

COD 9009234 is the kaolinite of Bish (1993), refined with its hydrogen
positions. Replicate the triclinic cell to about 31 × 21 × 15 Å:

.. code-block:: python

   kaol = AtomicSystem.from_cod(9009234)
   kaol = kaol.replicate([6, 4, 2])
   print(kaol.num_atoms, "atoms", kaol.atoms["type"].value_counts().to_dict())

.. code-block:: text

   816 atoms {'O': 432, 'H': 192, 'Al': 96, 'Si': 96}

Step 2: Type with the ClayFF rules
----------------------------------

.. code-block:: python

   kaol.set_types_from_elements()
   kaol.set_topology("clayff")
   print(kaol.atoms["type"].value_counts().to_dict(), kaol.bond_types)

.. code-block:: text

   {'Hh': 192, 'Oh': 192, 'Ob': 144, 'Al': 96, 'Si': 96, 'Osi': 96} ['Hh-Oh']

``Ob`` are the basal Si-O-Si oxygens, ``Osi`` the apical oxygens shared by Si
and the octahedral Al, ``Oh``/``Hh`` the hydroxyls. Only the O-H bonds are
created: ClayFF holds the rest of the structure with non-bonded terms.

All Al are six-fold, as in the raw mineral (27Al NMR: Al(VI) only):

.. code-block:: python

   print(analyze_coordination(kaol, center_types="Al").round(1))

.. code-block:: text

            0    1    2    3    4    5      6    7   8+  mean
   frame
   0      0.0  0.0  0.0  0.0  0.0  0.0  100.0  0.0  0.0   6.0

Step 3: Add the Al-O-H angles
-----------------------------

ClayFF-MOH needs an angle term on every Al-O-H. A
:class:`~cemd.topology.TopologyRule` creates angles between the neighbours of
a central atom: centred on the hydroxyl O, with its two Al and its H as
neighbours, it gives the Al-O-H angles, and also an Al-O-Al angle, which the
force field does not use:

.. code-block:: python

   moh = TopologyRule(
       "type Oh",
       [NeighborCriterion("type Al", 2.1, 2), NeighborCriterion("type Hh", 1.2, 1)],
       angles=True,
   )
   kaol.set_topology(moh)
   print(kaol.angle_types, kaol.num_angles)

   kaol.remove_connection_types(angle_types=["Al-Oh-Al"])
   print(kaol.angle_types, kaol.num_angles)

.. code-block:: text

   ['Al-Oh-Al', 'Al-Oh-Hh'] 576
   ['Al-Oh-Hh'] 384

Each of the 192 hydroxyls bridges two Al in the dioctahedral sheet, hence
384 Al-O-H angles.

Step 4: Isomorphic substitution
-------------------------------

Natural kaolinites hold some Fe\ :sup:`3+` in place of Al.
:meth:`~cemd.core.atomic_system.AtomicSystem.substitute_atoms` picks the sites
at random (``seed`` makes it reproducible) and, by default, never puts two
substituted atoms on the same bridging oxygen, Loewenstein's rule for
tetrahedra:

.. code-block:: python

   fe = kaol.substitute_atoms("Al", "Fe", fraction=0.05, seed=1)
   print(len(fe), "Fe:", fe)
   print(kaol.angle_types)

.. code-block:: text

   5 Fe: [ 69 137 324 477 579]
   ['Al-Oh-Hh', 'Fe-Oh-Hh']

The angles around the new Fe are relabelled ``Fe-Oh-Hh``. ClayFF-MOH has no
Fe-O-H term, so remove them:

.. code-block:: python

   kaol.remove_connection_types(angle_types=["Fe-Oh-Hh"])
   print(kaol.angle_types, kaol.num_angles)

.. code-block:: text

   ['Al-Oh-Hh'] 364

For a heterovalent substitution (Mg for Al in a smectite, Al for Si in a
mica), give the charge of the new atom with ``charge=`` and add the
compensating interlayer cations.

Step 5: Assign ClayFF-MOH
-------------------------

ClayFF uses a single ``ob`` type for every oxygen that is not a hydroxyl:

.. code-block:: python

   kaol.set_ff_from_database({
       "Si": "clayff_moh.st", "Al": "clayff_moh.ao", "Fe": "clayff_moh.feo",
       "Ob": "clayff_moh.ob", "Osi": "clayff_moh.ob",
       "Oh": "clayff_moh.oh", "Hh": "clayff_moh.ho",
   })
   print(f"Total charge: {kaol.total_charge:.3f} e")
   print(kaol.ff_params.bond)
   print(kaol.ff_params.angle)

.. code-block:: text

   Total charge: 0.000 e
   {'Hh-Oh': MorseBondParams(r0=0.945, D=132.2491, alpha=2.1815, ref='https://doi.org/10.1021/acs.jpcc.1c04600', model='clayff_moh')}
   {'Al-Oh-Hh': HarmonicAngleParams(k=15.0, theta0=116.0, ref='https://doi.org/10.1021/acs.jpcc.7b05362', model='clayff_moh')}

The O-H bond is a Morse potential (dioctahedral set) and the Al-O-H angle
takes the bulk value of 116°. Hydroxyls of an edge or a cleaved surface use
110° (see :doc:`ff_database`).

.. code-block:: python

   kaol.write("kaolinite.data")

.. code-block:: text

   Bond Coeffs

   1       132.2491     2.1815     0.9450   # Hh-Oh

   Angle Coeffs

   1        15.0000   116.0000   # Al-Oh-Hh

In the LAMMPS input, declare ``bond_style morse`` and ``angle_style harmonic``.
With water added (harmonic SPC bonds), the data file switches to
``Bond Coeffs # hybrid`` and the input needs ``bond_style hybrid harmonic morse``.

.. image:: /_static/images/tuto_kaolinite.png
   :alt: Kaolinite layers with the Fe substituted in the octahedral sheet
   :align: center
   :width: 700px

Two 1:1 layers seen along *b*: the tetrahedral sheet (Si in yellow) below,
the octahedral sheet (Al in blue, Fe in orange) with its hydroxyls above.
ClayFF draws no Al-O bond, they are held by non-bonded terms.
