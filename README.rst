====
cemd
====

**cemd** builds atomistic models of **cements, clays, carbonates and oxide
glasses**, and writes them out ready for
`LAMMPS <https://www.lammps.org/>`__.

It covers the whole path from a crystal structure to a simulation input:
assembling the solid, cutting surfaces and pores, filling them with a
solution, assigning atom types from force-field rules, resolving the
parameters, and analysing the result.

.. code-block:: python

   from cemd.build import CSHBuilder, SolutionBuilder, Splitter

   # A C-S-H matrix at Ca/Si = 1.5
   csh = CSHBuilder(cs_ratio=1.5, ws_ratio=1.0).build(model="tob11a_merlino.cif")
   csh.set_types_from_elements()
   csh.set_topology("cshff")

   # Open a 25 A nanopore along an interlayer and fill it with NaOH
   pore = SolutionBuilder(density=1.0, counts={"Na": 8, "HO": 8})
   system = Splitter(csh, coordinate=17.0, axis="z", gap_size=25.0) \
       .add_solution(pore, padding=2.0) \
       .split()

   system.set_ff_from_database()
   system.write("csh_nanopore.data")


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

**Type and parameterise.** ``set_topology()`` applies force-field typing
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


Installation
============

.. code-block:: bash

   pip install .               # library
   pip install ".[gui]"        # with the graphical interface

Python 3.11 or later. The graphical interface then starts with:

.. code-block:: bash

   cemd-gui

To work on cemd itself, install it editable, with the test and
documentation extras:

.. code-block:: bash

   pip install -e ".[dev,gui,docs]"

`Packmol <https://github.com/m3g/packmol>`__ must be installed and on your
``$PATH``: every builder that packs molecules into a volume shells out to
it. The Python dependencies (MDAnalysis, pymatgen, RDKit, NumPy, SciPy,
pandas, matplotlib, Dask) are installed automatically.


Documentation
=============

Build it locally with:

.. code-block:: bash

   sphinx-build -b html docs docs/_build/html

The user guide covers building, analysis and the force-field database, and
four worked tutorials go from a calcite surface to a C-S-H nanopore filled
with an alkaline pore solution.


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
