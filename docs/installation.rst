Installation
============

Prerequisites
-------------

What you need in place before installing **cemd**:

* **Python 3.11 or later** — required.
* **Conda** — recommended for creating the environment
  (`Miniconda or Anaconda <https://www.anaconda.com/docs/main>`__), though a
  plain ``venv`` works too.
* **VMD** — *optional*, and only for
  :meth:`~cemd.core.atomic_system.AtomicSystem.view`, which opens the system in
  VMD (`download <https://www.ks.uiuc.edu/Research/vmd/>`__). It must be on your
  ``$PATH``. The graphical interface renders with PyVista and needs no VMD.

`Packmol <https://github.com/m3g/packmol>`__ does not need installing: it is a
dependency of **cemd**, and pip puts its binary in your environment alongside
the Python packages. Install it yourself, and put it on your ``$PATH``, only if
no wheel matches your platform.


User Installation
-----------------

It is strongly recommended to install **cemd** in a dedicated virtual environment.

Create a virtual environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**cemd** pulls in pymatgen, MDAnalysis and RDKit, each of which pins versions of
its own. Installing it beside unrelated work asks pip to satisfy every project's
constraints at once, which is where incompatibilities appear. A dedicated
environment avoids the question entirely.

**Option 1: Using conda (recommended)**

.. code-block:: bash

   conda create -n cemd python=3.11
   conda activate cemd

**Option 2: Using venv (Python built-in)**

.. code-block:: bash

   python -m venv cemd
   source cemd/bin/activate          # Linux/macOS
   # cemd\Scripts\activate           # Windows

Install the package
^^^^^^^^^^^^^^^^^^^

Once your environment is active, install **cemd** from PyPI:

.. code-block:: bash

   pip install cemd

For GUI support:

.. code-block:: bash

   pip install "cemd[gui]"

That installs a ``cemd-gui`` command, which opens the interface:

.. code-block:: bash

   cemd-gui

Developer Installation
----------------------

If you plan to contribute to **cemd** development, clone the repository and install in editable mode:

.. code-block:: bash

   git clone https://github.com/cemd-dev/cemd.git
   cd cemd

Then create the development environment. This installs **cemd** editable,
with the test, GUI and documentation extras, so there is nothing to run
afterwards:

.. code-block:: bash

   conda env create -f environment.yml
   conda activate cemd

Editable means the environment points at your working copy: an edit to the
source is visible on the next import, with no reinstall.

Without conda, the same thing in a virtual environment:

.. code-block:: bash

   python -m venv .venv
   source .venv/bin/activate
   pip install -e ".[dev,gui,docs]"

.. note::

   Do not also add the repository to ``PYTHONPATH``. It would make cemd
   importable from every environment on the machine, including the throwaway
   ones used to check what a user actually receives -- and those checks then
   silently test the working copy instead of the installed package.

Verifying the Installation
--------------------------

Run the following to confirm that **cemd** is correctly installed:

.. code-block:: python

   import cemd
   from cemd import AtomicSystem
   print(f"CEMD {cemd.__version__} is ready.")