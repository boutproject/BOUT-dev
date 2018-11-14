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

in the top section of the ``BOUT.inp`` file and enable the locations you will
use with a call to mesh->addCoordinates(location) (see below).  The
**test-staggered** example illustrates how to use staggered grids in BOUT++.

There are four possible locations in a grid cell where a quantity can be
defined in BOUT++: centre, lower X, lower Y, and lower Z. These are
illustrated in :numref:`staggergrids-location`.

.. _staggergrids-location:
.. figure:: ../figs/stagLocations.*
   :alt: Staggered grid cell locations

   The four possible cell locations for defining quantities

The possible locations are specified with the `CELL_LOC` type, which has the
possible values `CELL_CENTRE`, `CELL_XLOW`, `CELL_YLOW`, or `CELL_ZLOW`.
`CELL_CENTRE` is enabled by default, but a call to
mesh->addCoordinates(location) is required to enable the others.

The key lines in the **staggered_grid** example which specify the
locations of the evolving variables are::

    Field3D n, v;

    int init(bool restart) {

      mesh->addCoordinates(CELL_YLOW);

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

.. note:: For advanced users:
          If you change members of the Coordinates object manually, you should
          change the CELL_CENTRE Coordinates and only call addCoordinates() for
          other locations after you call Coordinates::geometry() on the
          CELL_CENTRE Coordinates. Then the Coordinates at staggered locations
          will be interpolated from the correct, final CELL_CENTRE version.

          The example uses the global Mesh object 'mesh'. If you are using any
          other Mesh objects, you need to initialize the Coordinates objects in
          their coords_map members by calling the addCoordinates(CELL_LOC
          location) method for each location you will use.

          An exception will be thrown if you call Coordinates::geometry() from
          any Coordinates object at a staggered location, since these are
          expected to be consistent with (and calculated from) the CELL_CENTRE
          Coordinates. If you need to change them, you can set the option
          mesh:allow_geometry_without_recalculate_staggered=true to disable
          this check; you must then ensure that all the Coordinates objects in
          Mesh::coords_map are consistent with each other. Setting this option
          also allows changes to be made and geometry() to be called on the
          CELL_CENTRE Coordinates after other locations have been added to
          coords_map; in this case you will need to update the other locations
          explicitly by calling Mesh::addCoordinates(location, true) - the
          optional second argument causes addCoordinates to overwrite any
          existing Coordinates object at location and replace it with one
          calculated from the current CELL_CENTRE Coordinates.


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

.. note:: The region argument is optional but useful (see :ref:sec_iterating
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
cell centres. As with the arithmetic operators, if you ask for the
result to be staggered in a different direction from the input then
the differencing will be to cell centre and then be interpolated. For
example ``Grad_par(v, CELL_XLOW)`` would first perform staggered
differencing from `CELL_YLOW` to get a result at `CELL_CENTRE` , and
then interpolate the result to `CELL_XLOW` .

Advection operators which take two arguments return a result which is
defined at the location of the field being advected. For example
``Vpar_Grad_par(v, f)`` calculates :math:`v \nabla_{||} f` and returns a
result at the same location as ``f``. If ``v`` and ``f`` are defined at
the same locations then centred differencing is used, if one is centred
and the other staggered then staggered differencing is used, and if both
are staggered to different locations then the behaviour is less well
defined (don’t do it). As with other differential operators, the
required location of the result can be given as an optional argument.
