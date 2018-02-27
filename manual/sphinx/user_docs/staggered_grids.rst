.. _sec-staggergrids:

Staggered grids
===============

Until now all quantities have been cell-centred i.e. both velocities and
conserved quantities were defined at the same locations. This is because
these methods are simple and this was the scheme used in the original
BOUT. This class of methods can however be susceptible to grid-grid
oscillations, and so most shock-capturing schemes involve densities and
velocities (for example) which are not defined at the same location:
their grids are staggered.

By default BOUT++ runs with all quantities at cell centre. To enable
staggered grids, set

::

    StaggerGrids = true

in the top section of the ``BOUT.inp`` file. The **test-staggered**
example illustrates how to use staggered grids in BOUT++.

There are four possible locations in a grid cell where a quantity can be
defined in BOUT++: centre, lower X, lower Y, and lower Z. These are
illustrated in :numref:`staggergrids-location`.

.. _staggergrids-location:
.. figure:: ../figs/stagLocations.*
   :alt: Staggered grid cell locations

   The four possible cell locations for defining quantities

To specify the location of a variable, use the method
``setLocation()`` with one of the locations ``CELL_CENTRE``,
``CELL_XLOW``, ``CELL_YLOW`` , or ``CELL_ZLOW`` .

The key lines in the **test-staggered** example which specify the
locations of the evolving variables are

::

    Field3D n, v;

    int physics_init(bool restart) {
      v.setLocation(CELL_YLOW); // Staggered relative to n
      SOLVE_FOR2(n, v);
      ...

which makes the velocity ``v`` staggered to the lower side of the cell
in Y, whilst the density :math:`n` remains cell centred.

.. note:: If BOUT++ was configued ``--with-checks``,
          ``Field3D::setLocation`` will throw an exception if you
          don't have staggered grids turned on and try to set the
          location to something other than ``CELL_CENTRE``. If you
          want to be able to run your model with and without staggered
          grids, you should do something like::

            if (v.getMesh()->StaggerGrids) {
              v.setLocation(CELL_YLOW);
            }

          Compiling BOUT++ with checks turned off will instead cause
          ``Field3D::setLocation`` to silently set the location to
          ``CELL_CENTRE`` if staggered grids are off, regardless of
          what you pass it.


Arithmetic operations between staggered quantities are handled by
interpolating them to the same location according to the algorithm in
:numref:`fig-stagArith`.

.. _fig-stagArith:
.. figure:: ../figs/stagArith.*
   :alt: How the cell location of an arithmetic operation is decided

   How the cell location of an arithmetic operation (``+, -, *, /,
   \pow``) is decided


If performing an operation between variables defined at two different
locations, the order of the variables matter: the result will be defined
at the locations of the **left** variable. For example, ``n*v`` would be
``CELL_CENTRE`` because this is the location of ``n`` , whilst ``v*n``
would be ``CELL_YLOW`` . Relying on this behaviour could lead to
trouble, to make your code clearer it’s probably best to use the
interpolation routines. Include the header file

::

    #include <interpolation.hxx>

then use the ``interp_to(field, location)`` function. Using this,
``interp_to(n, CELL_YLOW)*v`` would be ``CELL_YLOW`` as ``n`` would be
interpolated.

Differential operators by default return fields which are defined at the
same location as their inputs, so here ``Grad_par(v)`` would be
``CELL_YLOW`` . If this is not what is wanted, give the location of the
result as an additional argument: ``Grad_par(v, CELL_CENTRE)`` uses
staggered differencing to produce a result which is defined at the cell
centres. As with the arithmetic operators, if you ask for the result to
be staggered in a different direction from the input then the
differencing will be to cell centre and then be interpolated. For
example ``Grad_par(v, CELL_XLOW)`` would first perform staggered
differencing from ``CELL_YLOW`` to get a result at ``CELL_CENTRE`` , and
then interpolate the result to ``CELL_XLOW`` .

Advection operators which take two arguments return a result which is
defined at the location of the field being advected. For example
``Vpar_Grad_par(v, f)`` calculates :math:`v \nabla_{||} f` and returns a
result at the same location as ``f``. If ``v`` and ``f`` are defined at
the same locations then centred differencing is used, if one is centred
and the other staggered then staggered differencing is used, and if both
are staggered to different locations then the behaviour is less well
defined (don’t do it). As with other differential operators, the
required location of the result can be given as an optional argument.
