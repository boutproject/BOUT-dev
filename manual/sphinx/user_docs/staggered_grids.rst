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
staggered grids, set::

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
`Field3D::setLocation` with one of the `CELL_LOC` locations
`CELL_CENTRE`, `CELL_XLOW`, `CELL_YLOW`, or `CELL_ZLOW`.

The key lines in the **staggered_grid** example which specify the
locations of the evolving variables are::

    Field3D n, v;

    int init(bool restart) {
      v.setLocation(CELL_YLOW); // Staggered relative to n
      SOLVE_FOR(n, v);
      ...

which makes the velocity ``v`` staggered to the lower side of the cell
in Y, whilst the density :math:`n` remains cell centred.

.. note:: If BOUT++ was configued ``--with-checks``,
          `Field3D::setLocation` will throw an exception if you don't
          have staggered grids turned on and try to set the location
          to something other than `CELL_CENTRE`. If you want to be
          able to run your model with and without staggered grids, you
          should do something like::

            if (v.getMesh()->StaggerGrids) {
              v.setLocation(CELL_YLOW);
            }

          Compiling BOUT++ with checks turned off will instead cause
          `Field3D::setLocation` to silently set the location to
          `CELL_CENTRE` if staggered grids are off, regardless of what
          you pass it.


Arithmetic operations can only be performed between variables with the same
location. When performing a calculation at one location, to include a variable
from a different location, use the interpolation routines. Include the header
file

::

    #include <interpolation.hxx>

then use the ``interp_to(field, location, region)`` function. For example,
given a `CELL_CENTRE` field ``n`` and a `CELL_YLOW` field ``v``, to calculate
``n*v`` at `CELL_YLOW`, call ``interp_to(n, CELL_YLOW)*v`` whose result will be
`CELL_YLOW` as ``n`` is interpolated.

.. note:: The region argument is optional but useful (see :ref:`sec-iterating`
          for more on regions). The default `RGN_ALL` reproduces the historical
          behaviour of BOUT++, which communicates before returning the result
          from ``interp_to``. Communication is necessary because the result of
          interpolation in the guard cells depends on data from another process
          (except, currently, in the case of interpolation in the z-direction
          which can be done without communication because all the z-points are
          on the same process).

          Using RGN_NOBNDRY no communication is performed
          (so interp_to is faster, potentially significantly faster when using
          many processes) and all the guard cells are invalid. Whichever region
          is used, the boundary guard cells are invalid since no boundary
          condition is applied in interp_to. If the guard cells are needed
          (e.g. to calculate a derivative) a boundary condition must be applied
          explicitly to the result.

          RGN_NOX and RGN_NOY currently have identical behaviour to RGN_ALL
          because at present BOUT++ has no functions for single-direction
          communication which could in principle be used in these cases (if the
          combination of region and direction of interpolation allows it). x-
          or y-interpolation can never be calculated in guard cells without
          communication because the corner guard cells are never valid.

Differential operators by default return fields which are defined at
the same location as their inputs, so here ``Grad_par(v)`` would be
`CELL_YLOW` . If this is not what is wanted, give the location of the
result as an additional argument: ``Grad_par(v, CELL_CENTRE)`` uses
staggered differencing to produce a result which is defined at the
cell centres. It is an error to ask for the result to be staggered in
a different direction from the input as the best that could be done
would be to calculate output at ``CELL_CENTRE`` and then interpolate
this to the requested location, but the interpolation would in general
require boundary conditions to be applied first.

Advection operators which take two arguments return a result which is
defined at the location of the field being advected. For example
``Vpar_Grad_par(v, f)`` calculates :math:`v \nabla_{||} f` and returns a
result at the same location as ``f``. If ``v`` and ``f`` are defined at
the same locations then centred differencing is used, if one is centred
and the other staggered then staggered differencing is used; it is an
error for both to be staggered to different locations. As with other
differential operators, the required location of the result can be
given as an optional argument, but at least for now it is an error for
this to be different from the location of the field being advected
(``f`` here).
