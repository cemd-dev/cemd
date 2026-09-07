====
cemd
====

|pypi| |tests| |docs| |python| |licence|

**cemd** builds atomistic models of **cements, clays, carbonates and oxide
glasses**, and writes them out ready for
`LAMMPS <https://www.lammps.org/>`__.

It covers the whole path from a crystal structure to a simulation input:
assembling the solid, cutting surfaces and pores, filling them with a
solution, assigning atom types from force-field rules, resolving the
parameters, and analyzing the result.

Whether the system is a clay interlayer, a calcite surface in brine, an
oxide glass or a cement hydrate, the pieces are the same: a solid, a
solution, a set of atom types and a LAMMPS data file at the end.


What it does
============

**Build.** Every builder is a blueprint you configure once and reuse:

- ``CSHBuilder`` — C-S-H and C-A-S-H from tobermorite, at a target Ca/Si
  and H2O/Si ratio
- ``AFBuilder`` — AFm/AFt phases
- ``SurfaceBuilder`` — surface slabs from Miller indices, ranked by how
  many bonds each termination breaks
- ``SolutionBuilder`` — electrolytes from molarities or explicit counts,
  packed to a target density
- ``GlassBuilder`` — amorphous oxide melts
- ``Splitter`` — cut a system open along an axis and fill the gap with a
  solution

**Type and parameterize.** ``set_topology()`` applies force-field typing
rules (ClayFF, CSHFF and your own), ``set_ff_from_database()`` resolves
the parameters. The bundled database holds ClayFF, CSHFF2014, IFF
(CHARMM and CVFF flavours), GROMOS 54A7/ATB, Raiteri 2015, Guillot &
Sator 2007, and the usual water models (SPC, SPC/E, SPC/Fw, TIP3P,
TIP4P-2005, TIP4P-EW).

**Analyse.** Silicate network statistics (Ca/Si, H2O/Si, Q\ :sup:`n`
distribution, mean chain length), radial distribution functions, density
and diffusion profiles, mean-squared displacement, electrostatic
potential.

**Read and write.** Reads LAMMPS data files (``.data``, ``.lmp``), PDB,
CIF, moltemplate (``.lt``) and SDF; writes LAMMPS data files and PDB. It
also fetches structures directly from
`PubChem <https://pubchem.ncbi.nlm.nih.gov/>`__ and the
`COD <https://www.crystallography.net/cod/>`__, builds them from SMILES,
and converts to and from `MDAnalysis <https://www.mdanalysis.org/>`__ and
`pymatgen <https://pymatgen.org/>`__ objects.

**Inspect.** An optional PySide6/PyVista interface for viewing and editing
structures interactively.


Two examples
============

A mineral surface in contact with an electrolyte -- the archetypal
interface problem:

.. code-block:: python

   from cemd import AtomicSystem
   from cemd.build import SolutionBuilder, SurfaceBuilder

   calcite = AtomicSystem.from_file("calcite.cif")
   slabs, *_ = SurfaceBuilder(calcite).build((1, 0, 4), min_slab_size=12.0)

   brine = SolutionBuilder(density=1.0, molarities={"Na": 0.5, "Cl": 0.5})
   system = slabs[0].add_liquid_layer(brine, thickness=20.0)

   system.set_topology("clayff")
   system.write("calcite_brine.data")

And a cement hydrate, where cemd does what nothing else does -- build the
C-S-H at a target Ca/Si, then measure where it can be cut before opening
a pore, since under CSHFF the silicate framework carries no explicit bond
to guide you:

.. code-block:: python

   from cemd.build import CSHBuilder, SolutionBuilder, Splitter

   csh = CSHBuilder(cs_ratio=1.5, ws_ratio=1.0).build(model="tob11a_merlino.cif")
   csh.set_types_from_elements()
   csh.set_topology("cshff")

   pore = SolutionBuilder(density=1.0, counts={"Na": 8, "HO": 8})
   system = Splitter(csh, coordinate=17.0, axis="z", gap_size=25.0) \
       .add_solution(pore, padding=2.0) \
       .split()

   system.write("csh_nanopore.data")


Installation
============

Python 3.11 or later. Install into an environment of its own, so cemd and
its dependencies stay out of the way of your other work:

.. code-block:: bash

   python -m venv cemd_env
   source cemd_env/bin/activate          # cemd_env\Scripts\activate on Windows

   pip install cemd

For the graphical interface, which then opens with ``cemd-gui``:

.. code-block:: bash

   pip install "cemd[gui]"
   cemd-gui

To work on cemd itself, clone it and install editable, which points the
environment at your working copy rather than copying the code:

.. code-block:: bash

   git clone https://github.com/cemd-dev/cemd.git
   cd cemd
   pip install -e ".[dev,gui,docs]"

The `installation guide <https://cemd-dev.github.io/cemd/installation.html>`__
covers conda environments and the prerequisites in more detail.

Everything comes with it, including
`Packmol <https://github.com/m3g/packmol>`__ -- every builder that packs
molecules into a volume shells out to that binary, and pip puts it in the
environment alongside the Python dependencies (MDAnalysis, pymatgen,
RDKit, NumPy, SciPy, pandas, matplotlib, Dask). Should no wheel match your
platform, install Packmol yourself and put it on your ``$PATH``.


Documentation
=============

**https://cemd-dev.github.io/cemd/**

The user guide covers building, analysis and the force-field database, and
four worked tutorials go from a calcite surface to a C-S-H nanopore filled
with an alkaline pore solution. Every output shown in them was captured
from a real run.

To build the documentation locally:

.. code-block:: bash

   sphinx-build -b html docs docs/_build/html


Tests
=====

.. code-block:: bash

   pytest

The suite runs against the real bundled data and force fields rather than
mocks. Tests that shell out to Packmol are skipped when it is not
installed.


Citing
======

See ``CITATION.cff``.


License
=======

GPL-3.0-only. See ``LICENCE``.


.. |pypi| image:: https://img.shields.io/pypi/v/cemd.svg
   :target: https://pypi.org/project/cemd/
   :alt: PyPI

.. |tests| image:: https://github.com/cemd-dev/cemd/actions/workflows/tests.yml/badge.svg
   :target: https://github.com/cemd-dev/cemd/actions/workflows/tests.yml
   :alt: Tests

.. |docs| image:: https://github.com/cemd-dev/cemd/actions/workflows/docs.yml/badge.svg
   :target: https://cemd-dev.github.io/cemd/
   :alt: Documentation

.. |python| image:: https://img.shields.io/badge/python-3.11%2B-blue.svg
   :target: https://www.python.org/downloads/
   :alt: Python 3.11+

.. |licence| image:: https://img.shields.io/badge/licence-GPL--3.0-blue.svg
   :target: https://github.com/cemd-dev/cemd/blob/main/LICENCE
   :alt: Licence: GPL-3.0-only
