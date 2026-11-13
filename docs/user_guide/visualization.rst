Visualizing and Rendering
=========================

cemd shows a system in `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_ with a
ready-made scene: element colours, CPK spheres sized by element, dynamic bonds
for the covalent units, orthographic projection, white background and the
periodic box. The pictures of the tutorials were rendered from that scene.

Opening a system
----------------

.. code-block:: python

   from cemd import AtomicSystem

   system = AtomicSystem.from_file("cash.data")
   system.view()

   # With a trajectory read by MDAnalysis (DCD, XTC, LAMMPS dump...)
   system.view(trajectory="production.dcd")

VMD must be on the ``PATH``. ``view()`` writes the system to a temporary
LAMMPS data file, loads it with TopoTools and applies the representations; the
Python call returns when VMD is closed.

Commands available in the VMD console
-------------------------------------

The scene defines two Tcl procedures, typed in the VMD console (*Extensions >
Tk Console*):

``render_tachyon [file] [resx] [resy]``
   Renders the current view with Tachyon, with ambient occlusion and 12
   antialiasing samples, to a PNG. Defaults: ``output.png``, 2048 × 2048 px.

   .. code-block:: tcl

      rotate x by -90
      rotate y by -90
      render_tachyon cash.png 1600 1600

``move_by x y z``
   Shifts every atom of every frame by a vector (Å) and wraps them back in the
   box: the way to bring a pore, an interface or a molecule to the middle of the
   picture.

   .. code-block:: tcl

      move_by 0 0 10

Orienting the view
------------------

The view starts along *c*, *a* to the right. Rotations are about the screen
axes and add up:

- ``rotate x by -90``: *c* vertical, looking along *b*, the usual view of a
  layered mineral or a slab;
- then ``rotate y by -90``: looking along *a*, the view of the C-S-H chains in
  :doc:`tuto_cash`;
- ``scale by 1.2`` zooms in, ``display resetview`` goes back to the start
  (it also cancels the rotations).

Cropping the picture
--------------------

Tachyon renders the whole window, with a wide white margin. Trim it with
ImageMagick:

.. code-block:: bash

   magick cash.png -fuzz 5% -trim +repage cash_cropped.png    # ImageMagick 7
   convert cash.png -fuzz 5% -trim +repage cash_cropped.png   # ImageMagick 6

``-fuzz 5%`` treats near-white pixels as background, so the antialiased edges
are trimmed too. For the web, ``pngquant`` reduces the file two to three times
without visible change:

.. code-block:: bash

   pngquant --quality=80-95 --output cash_small.png cash_cropped.png
