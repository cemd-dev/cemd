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

Create a dedicated conda environment:

.. code-block:: bash

   conda env create -f environment.yml
   conda activate cemd

Then install the package in editable mode:

.. code-block:: bash

   pip install -e .

For GUI development, use the dedicated environment file:

.. code-block:: bash

   conda env create -f environment_gui.yml
   conda activate cemd_ui
   pip install -e ".[gui]"

Verifying the Installation
--------------------------

Run the following to confirm that **cemd** is correctly installed:

.. code-block:: python

   import cemd
   from cemd import AtomicSystem
   print(f"CEMD {cemd.__version__} is ready.")