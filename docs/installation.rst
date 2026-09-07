Installation
============

Prerequisites
-------------

The following software must be installed before setting up **cemd**:

* **Python 3.11** — required Python version.
* **Conda** — package and environment manager (`Miniconda or Anaconda <https://www.anaconda.com/docs/main>`__).
* **Packmol** — required for automated system construction (`download <https://m3g.github.io/packmol/>`__). Must be accessible in your ``$PATH``.
* **VMD** — required for system visualization (`download <https://www.ks.uiuc.edu/Research/vmd/>`__). Must be accessible in your ``$PATH``.


User Installation
-----------------

It is strongly recommended to install **cemd** in a dedicated virtual environment.

Create a virtual environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Choose the method that suits you best:

**Option 1: Using venv (Python built-in)**

.. code-block:: bash

   python -m venv cemd_env
   source cemd_env/bin/activate      # Linux/macOS
   # cemd_env\Scripts\activate       # Windows

**Option 2: Using conda**

.. code-block:: bash

   conda create -n cemd_env python=3.11
   conda activate cemd_env

Install the package
^^^^^^^^^^^^^^^^^^^

Once your environment is active, install **cemd** from the repository:

.. code-block:: bash

   pip install "git+https://github.com/cemd-dev/cemd.git"

For GUI support:

.. code-block:: bash

   pip install "cemd[gui] @ git+https://github.com/cemd-dev/cemd.git"

.. note::

   **cemd** is not on PyPI yet, so it is installed from GitHub. Once it is
   published, ``pip install cemd`` and ``pip install "cemd[gui]"`` will
   work instead.

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