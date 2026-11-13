From C-S-H to C-A-S-H
=====================

In blended cements (slag, fly ash, calcined clays) part of the silicon of the
C-S-H is replaced by aluminium. The Al sits in the *bridging* tetrahedra of the
dreierkette chains, as 27Al and 29Si MAS NMR show: four-fold Al(IV) linked to
two Si, which appear as Q\ :sup:`2`\ (1Al) in the 29Si spectrum.

This tutorial builds a C-S-H, substitutes Al for bridging Si with
:meth:`~cemd.build.CSHBuilder.to_cash`, checks the result the way an NMR
spectrum would (chain length, Q\ :sup:`n`, Al coordination), and assigns a
force field.

.. note::

   The C-S-H builder places the interlayer water and calcium with Packmol and
   removes bridging tetrahedra at random: two runs give different structures.
   The numbers below come from one run; the ratios are fixed by the input,
   the Al coordination of the unrelaxed structure varies from run to run.

Prerequisites
-------------

.. code-block:: python

   from cemd.build import CSHBuilder
   from cemd.analysis import analyze_coordination, analyze_silicates

Packmol must be installed (see :doc:`../installation`).

Step 1: Build and type the C-S-H
--------------------------------

Start from the 11 Å tobermorite of Merlino, with Ca/Si = 1.7 and one water
per silicon. The CSH-FF rules then separate intralayer (``Ca``) and interlayer
(``Cw``) calcium, bridging (``Ob``) and non-bridging (``Osi``) oxygens,
silanols, hydroxides and water:

.. code-block:: python

   builder = CSHBuilder(cs_ratio=1.7, ws_ratio=1.0)
   csh = builder.build()

   csh.set_types_from_elements()
   csh.set_topology("cshff")
   print(csh.atoms["type"].value_counts())

.. code-block:: text

   type
   Osi     480
   Hw      228
   Si      180
   Ca      160
   Cw      146
   Oh      130
   Hh      130
   Ob      120
   Ow      114
   Osih      2
   Hsi       2

Check the composition and the silicate network:

.. code-block:: python

   r = analyze_silicates(csh)
   print(f"Ca/Si = {r['Ca/Si']:.2f}, H2O/Si = {r['H2O/Si']:.2f}, MCL = {r['MCL']:.2f}")
   print(f"Qn (%) = {r['Qn_distribution'].round(1)}")

.. code-block:: text

   Ca/Si = 1.70, H2O/Si = 1.00, MCL = 3.00
   Qn (%) = [ 0.  66.7 33.3  0.   0. ]

At Ca/Si = 1.7 two bridging tetrahedra out of three are missing: the chains
are mostly dimers (Q\ :sup:`1`) with some pentamers, hence MCL = 3. The water
count includes the hydroxides as half a water each (H\ :sub:`2`\ O = H / 2).

.. image:: /_static/images/tuto_csh_side.png
   :alt: C-S-H seen along a, with the silicate chains of the two main layers
   :align: center
   :width: 600px

The view is along *a*: the two calcium silicate layers run horizontally, with
the dreierkette chains (Si in yellow) on both faces of each calcium sheet
(green) and the interlayer water, hydroxide and calcium between them.

Step 2: Substitute aluminium
----------------------------

:meth:`~cemd.build.CSHBuilder.to_cash` takes the target Al/Si ratio of the
C-A-S-H. ``compensation="OH"`` removes, for each Al\ :sup:`3+` put in place of
a Si\ :sup:`4+`, the nearest interlayer hydroxide, so Ca/(Si+Al) keeps its
C-S-H value:

.. code-block:: python

   builder = CSHBuilder.from_system(csh)
   cash = builder.to_cash(al_si=0.1, compensation="OH", seed=1)
   print("Al:", cash.get_count("Al"), " Si:", cash.get_count("Si"), " atoms:", cash.num_atoms)

.. code-block:: text

   Al: 16  Si: 164  atoms: 1660

A system that has not been typed yet (straight from ``build()``) is typed with
the CSH-FF rules first, so the result always carries CSH-FF types. The
``seed`` fixes which bridging sites are chosen.

.. code-block:: python

   r = analyze_silicates(cash)
   print(f"Ca/(Si+Al) = {r['Ca/(Si+Al)']:.2f}, Al/Si = {r['Al/Si']:.2f}, MCL = {r['MCL']:.2f}")
   print(f"Qn (%) = {r['Qn_distribution'].round(1)}, Q2(1Al) = {r['Q2(1Al)']:.1f} %")

.. code-block:: text

   Ca/(Si+Al) = 1.70, Al/Si = 0.10, MCL = 3.00
   Qn (%) = [ 0.  73.2 26.8  0.   0. ], Q2(1Al) = 19.5 %

Q\ :sup:`2` drops because the bridging Si became Al, and the Si next to
them are counted as Q\ :sup:`2`\ (1Al). The mean chain length uses
Richardson's formula for C-A-S-H,
MCL = 2 [Q\ :sup:`1` + Q\ :sup:`2`\ (0Al) + 3/2 Q\ :sup:`2`\ (1Al)] / Q\ :sup:`1`,
which counts the bridging Al as part of the chain: replacing a tetrahedron
does not shorten the chain, and MCL stays at 3.

.. image:: /_static/images/tuto_cash_side.png
   :alt: The same structure after substitution, with Al in bridging positions
   :align: center
   :width: 600px

The Al (blue) occupy bridging positions, between two paired tetrahedra.

Step 3: Aluminium coordination
------------------------------

:func:`~cemd.analysis.analyze_coordination` counts the oxygens around each Al,
the number 27Al NMR measures:

.. code-block:: python

   print(analyze_coordination(cash, center_types="Al").round(1))

.. code-block:: text

            0    1    2    3     4     5    6    7   8+  mean
   frame
   0      0.0  0.0  0.0  0.0  43.8  56.2  0.0  0.0  0.0   4.6

In the unrelaxed structure some Al already have a fifth oxygen, a water
or hydroxide of the interlayer within the 2.4 Å cutoff. The fraction changes
from one build to the next and has no meaning before relaxation: compare with
NMR after an energy minimization and MD, by passing the trajectory
(an ``MDAnalysis.Universe``) instead of the system. Every frame then gets a row.

Step 4: Force field
-------------------

CSH-FF has no aluminium. Take the tetrahedral Al of ClayFF and the ClayFF
hydroxide, and CSH-FF for everything else:

.. code-block:: python

   assignments = {
       "Si": "cshff2014.si", "Osi": "cshff2014.o", "Ob": "cshff2014.ob",
       "Osih": "cshff2014.oh", "Hsi": "cshff2014.h",
       "Ca": "cshff2014.ca", "Cw": "cshff2014.cw",
       "Ow": "cshff2014.ow", "Hw": "cshff2014.hw",
       "Oh": "clayff.oh", "Hh": "clayff.ho",
       "Al": "clayff.at",
   }
   cash.set_ff_from_database(assignments)
   print(f"Total charge: {cash.total_charge:.3f} e")

.. code-block:: text

   Total charge: 55.810 e

Both models use the arithmetic (Lorentz-Berthelot) rule, so the cross terms
between Al and the CSH-FF oxygens are mixed from their self terms (see
:doc:`ff_database`).

The base CSH-FF charges were fitted on a tobermorite of Ca/Si = 1; at
Ca/Si = 1.7 the structure carries 55.8 e. Qomi et al. (2014) keep the cations
and spread the excess over the silicate oxygens. Doing the same moves Osi from
-1.14 to about -1.23 e and Ob from -1.00 to about -1.09 e, inside the ranges
they report for their C-S-H models:

.. code-block:: python

   cash.neutralize_charge({"Osi": 1.0, "Ob": 1.0})
   print(f"Total charge: {cash.total_charge:.3f} e")
   cash.write("cash.data")

.. code-block:: text

   Total charge: -0.000 e

.. tip::

   Look at the result before simulating it: ``cash.view()`` opens it in VMD.
   :doc:`visualization` shows how the pictures of this page were rendered.

Going further
-------------

- Relax the structure (energy minimization, then NPT MD), and follow the Al
  coordination along the trajectory with ``analyze_coordination``.
- Cut a pore in the C-A-S-H and fill it with a NaCl solution, as in
  :doc:`tuto_csh_nanopore`.
