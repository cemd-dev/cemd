CEMD Documentation
==================

Atomistic models of **cements, clays, carbonates and oxide glasses**, built
for `LAMMPS <https://www.lammps.org/>`__.

CEMD covers the whole path from a crystal structure to a simulation input:
assembling the solid, cutting surfaces and pores, filling them with a
solution, assigning atom types from force-field rules, resolving the
parameters, and analyzing the result.

.. grid:: 2
   :gutter: 3

   .. grid-item-card:: :material-regular:`construction;2em` Installation
      :link: installation
      :link-type: doc
      :text-align: center

      Step-by-step instructions to set up CEMD and its dependencies
      on your system.

      +++
      .. button-ref:: installation
         :expand:
         :color: primary
         :click-parent:

         To the installation guide

   .. grid-item-card:: :material-regular:`rocket_launch;2em` Getting Started
      :link: getting_started
      :link-type: doc
      :text-align: center

      New to CEMD? Start here for a quick introduction to the main
      concepts and your first atomic system.

      +++
      .. button-ref:: getting_started
         :expand:
         :color: primary
         :click-parent:

         To the getting started guide

   .. grid-item-card:: :material-regular:`auto_stories;2em` User Guide
      :link: user_guide/index
      :link-type: doc
      :text-align: center

      In-depth explanations of CEMD concepts — building systems,
      applying force fields, and running analyses.

      +++
      .. button-ref:: user_guide/index
         :expand:
         :color: primary
         :click-parent:

         To the user guide

   .. grid-item-card:: :material-regular:`biotech;2em` API Reference
      :link: api/index
      :link-type: doc
      :text-align: center

      Detailed description of all classes, methods and functions
      in CEMD. Assumes familiarity with the key concepts.

      +++
      .. button-ref:: api/index
         :expand:
         :color: primary
         :click-parent:

         To the API reference

.. toctree::
   :hidden:
   :caption: Contents:

   installation
   getting_started
   user_guide/index
   api/index