.. _sec-gridgen:

Generating input grids
======================

The simulation mesh describes the number and topology of grid points,
the spacing between them, and the coordinate system. For many problems,
a simple mesh can be created using options.

.. code-block:: bash

    [mesh]
    nx = 260  # X grid size
    ny = 256  # Y grid size

    dx = 0.1  # X mesh spacing
    dy = 0.1  # Y mesh spacing

The above options will create a :math:`260\times 256` mesh in X and Y
(MZ option sets Z resolution), with mesh spacing of :math:`0.1` in both
directions. By default the coordinate system is Cartesian (metric tensor
is the identity matrix), but this can be changed by specifying the
metric tensor components.

Integer quantities such as ``nx`` must be numbers (like “260”), not
expressions (like “256 + 2\*MXG”). Real (floating-point) values can be
expressions, allowing quite complicated analytic inputs. For example in
the example ``test-griddata``:

.. code-block:: bash

    # Screw pinch

    rwidth = 0.4

    Rxy = 0.1 + rwidth*x  # Radius from axis     [m]
    L   = 10              # Length of the device [m]

    dy = L/ny
    hthe = 1.0

    Zxy = L * y / (2*pi)

    Bpxy = 1.0      # Axial field [T]
    Btxy = 0.1*Rxy  # Azimuthal field [T]
    Bxy = sqrt(Btxy^2 + Bpxy^2)

    dr = rwidth / nx
    dx = dr * Bpxy * Rxy

These expressions use the same mechanism as used for variable
initialisation (:ref:`sec-expressions`): ``x`` is a variable from
:math:`0` to :math:`1` in the domain which is uniform in index space;
``y`` and ``z`` go from :math:`0` to :math:`2\pi`. As with variable
initialisation, common trigonometric and mathematical functions can be
used. In the above example, some variables depend on each other, for
example ``dy`` depends on ``L`` and ``ny``. The order in which these
variables are defined doesn’t matter, so ``L`` could be defined below
``dy``, but circular dependencies are not allowed. If the variables are
defined in the same section (as ``dy`` and ``L``) then no section prefix
is required. To refer to a variable in a different section, prefix the
variable with the section name e.g. “``section:variable``”.

More complex meshes can be created by supplying an input grid file to
describe the grid points, geometry, and starting profiles. Currently
BOUT++ supports either NetCDF, HDF5 format binary files. During startup,
BOUT++ looks in the grid file for the following variables. If any are
not found, a warning will be printed and the default values used.

-  X and Y grid sizes (integers) ``nx`` and ``ny`` **REQUIRED**

-  Differencing quantities in 2D arrays ``dx[nx][ny]`` and
   ``dy[nx][ny]``. If these are not found they will be set to 1.

-  Diagonal terms of the metric tensor :math:`g^{ij}` ``g11[nx][ny]``,
   ``g22[nx][ny]``, and ``g33[nx][ny]``. If not found, these will be set
   to 1.

-  Off-diagonal metric tensor :math:`g^{ij}` elements ``g12[nx][ny]``,
   ``g13[nx][ny]``, and ``g23[nx][ny]``. If not found, these will be set
   to 0.

-  Z shift for interpolation between field-aligned coordinates and
   shifted coordinates (see ``manual/coordinates.pdf``). Perpendicular
   differential operators are calculated in shifted coordinates when
   ``ShiftXderivs`` in ``mesh/mesh.hxx`` is enabled. ``ShiftXderivs``
   can be set in the root section of ``BOUT.inp`` as
   ``ShiftXderivs = true``. The shifts must be provided in the gridfile
   in a field ``zshift[nx][ny]``. If not found, ``zshift`` is set to
   zero.

The remaining quantities determine the topology of the grid. These are
based on tokamak single/double-null configurations, but can be adapted
to many other situations.

-  Separatrix locations ``ixseps1``, and ``ixseps2`` If neither is
   given, both are set to nx (i.e. all points in closed “core” region).
   If only ``ixseps1`` is found, ``ixseps2`` is set to nx, and if only
   ixseps2 is found, ixseps1 is set to -1.

-  Branch-cut locations ``jyseps1_1``, ``jyseps1_2``, ``jyseps2_1``, and
   ``jyseps2_2``

-  Twist-shift matching condition ``ShiftAngle[nx]`` for field aligned
   coordinates. This is applied in the “core” region between indices
   ``jyseps2_2``, and ``jyseps1_1 + 1``, if either ``TwistShift = True``
   enabled in the options file or in general the ``TwistShift`` flag in
   ``mesh/impls/bout/boutmesh.hxx`` is enabled by other means. BOUT++
   automatically reads the twist shifts in the gridfile if the shifts
   are stored in a field in a field ShiftAngle[nx]. If not given, this
   is set to zero.

The only quantities which are required are the sizes of the grid. If
these are the only quantities specified, then the coordinates revert to
Cartesian.

This section describes how to generate inputs for tokamak equilibria. If
you’re not interested in tokamaks then you can skip to the next section.

The directory ``tokamak_grids`` contains code to generate input grid
files for tokamaks. These can be used by the ``2fluid`` and
``highbeta_reduced`` modules, and are (mostly) compatible with inputs to
the BOUT-06 code.


BOUT++ Topology
---------------

Basic
~~~~~

In order to handle tokamak geometry BOUT++ contains an internal topology
which is determined by the branch-cut locations (``jyseps1_1``,
``jyseps1_2``, ``jyseps2_1``, and ``jyseps2_2``) and separatrix
locations (``ixseps1`` and ``ixseps2``).

The separatrix locations, ``ixseps1`` and ``ixseps2``, give the indices
in the ``x`` domain where the first and second separatrices are located.

If ``ixseps1 == ixseps2`` then there is a single separatrix representing
the boundary between the core region and the SOL region and the grid is
a connected double null configuration. If ``ixseps1 > ixseps2`` then
there are two separatrices and the inner separatrix is ``ixseps2`` so
the tokamak is an upper double null. If ``ixseps1 < ixseps2`` then there
are two separatrices and the inner separatrix is ``ixseps1`` so the
tokamak is a lower double null.

In other words: Let us for illustrative purposes say that
``ixseps1 > ixseps2`` (see :numref:`fig-topology-cross-section`). Let
us say that we have a field ``f(x,y,z)`` with a global ``x``-index which
includes ghost points. ``f(x<=xseps1,y,z)``) will then be periodic in
the ``y``-direction, ``f(xspes1<x<=xseps2,y,z)``) will have boundary
condition in the ``y``-direction set by the lowermost ``ydown`` and
``yup``. If ``f(xspes2<x,y,z)``) the boundary condition in the
``y``-direction will be set by the uppermost ``ydown`` and ``yup``. As
for now, there is no difference between the two sets of upper and lower
``ydown`` and ``yup`` boundary conditions (unless manually specified,
see :ref:`sec-custom-BC`).

These values are set either in the grid file or in ``BOUT.inp``.
:numref:`fig-topology-cross-section` shows schematically how ``ixseps`` is
used.

The branch cut locations, ``jyseps1_1``, ``jyseps1_2``, ``jyseps2_1``,
and ``jyseps2_2``, split the ``y`` domain into logical regions defining
the SOL, the PFR (private flux region) and the core of the tokamak. This
is illustrated also in :numref:`fig-topology-cross-section`. If
``jyseps1_2 == jyseps2_1`` then the grid is a single null configuration,
otherwise the grid is a double null configuration.

.. _fig-topology-cross-section:
.. figure:: ../figs/topology_cross_section.*
   :alt: Cross-section of the tokamak topology used in BOUT++

   Deconstruction of a poloidal tokamak cross-section into logical
   domains using the parameters ``ixseps1``, ``ixseps2``,
   ``jyseps1_1``, ``jyseps1_2``, ``jyseps2_1``, and ``jyseps2_2``

Advanced
~~~~~~~~

The internal domain in BOUT++ is deconstructed into a series of
logically rectangular sub-domains with boundaries determined by the
``ixseps`` and ``jyseps`` parameters. The boundaries coincide with
processor boundaries so the number of grid points within each sub-domain
must be an integer multiple of ``ny/nypes`` where ``ny`` is the number
of grid points in ``y`` and ``nypes`` is the number of processors used
to split the y domain. Processor communication across the domain
boundaries is then handled internally. :numref:`fig-topology-schematic`
shows schematically how the different regions of a double null tokamak
with ``ixseps1 = ixseps2`` are connected together via communications.

.. note::
   To ensure that each subdomain follows logically, the
   ``jyseps`` indices must adhere to the following conditions:

    - ``jyseps1_1 > -1``
    - ``jyseps2_1 >= jyseps1_1 + 1``
    - ``jyseps1_2 >= jyseps2_1``
    - ``jyseps2_2 <= ny - 1``

   To ensure that communications work branch cuts must align with
   processor boundaries.

.. _fig-topology-schematic:
.. figure:: ../figs/topology_schematic.*

   Schematic illustration of domain decomposition and communication in
   BOUT++ with ``ixseps1 = ixseps2``

Implementations
~~~~~~~~~~~~~~~

In BOUT++ each processor has a logically rectangular domain, so any
branch cuts needed for X-point geometry (see
:numref:`fig-topology-schematic`) must be at processor boundaries.

In the standard “bout” mesh (``src/mesh/impls/bout/``), the
communication is controlled by the variables

::

    int UDATA_INDEST, UDATA_OUTDEST, UDATA_XSPLIT;
    int DDATA_INDEST, DDATA_OUTDEST, DDATA_XSPLIT;
    int IDATA_DEST, ODATA_DEST;

These control the behavior of the communications as shown in
:numref:`fig-boutmesh-comms`.

.. _fig-boutmesh-comms:
.. figure:: ../figs/boutmesh-comms.*
   :alt: Communication of guard cells in BOUT++

   Communication of guard cells in BOUT++. Boundaries in X have only
   one neighbour each, but boundaries in Y can be split into two,
   allowing branch cuts

In the Y direction, each boundary region (**U**\ p and **D**\ own in Y)
can be split into two, with ``0 <= x < UDATA_XSPLIT`` going to the
processor index ``UDATA_INDEST``, and ``UDATA_INDEST <= x < LocalNx`` going
to ``UDATA_OUTDEST``. Similarly for the Down boundary. Since there are
no branch-cuts in the X direction, there is just one destination for the
**I**\ nner and **O**\ uter boundaries. In all cases a negative
processor number means that there’s a domain boundary so no
communication is needed.

The communication control variables are set in the ``topology()``
function, in ``src/mesh/impls/bout/boutmesh.cxx`` starting around line
2056. First the function ``default_connections()`` sets the topology to
be a rectangle

To change the topology, the function ``set_connection`` checks that the
requested branch cut is on a processor boundary, and changes the
communications consistently so that communications are two-way and there
are no “dangling” communications.

3D variables
------------

BOUT++ was originally designed for tokamak simulations where the input
equilibrium varies only in X-Y, and Z is used as the axisymmetric
toroidal angle direction. In those cases, it is often convenient to have
input grids which are only 2D, and allow the Z dimension to be specified
independently, such as in the options file. The problem then is how to
store 3D variables in the grid file?

Two representations are now supported for 3D variables:

#. A Fourier representation. If the size of the toroidal domain is not
   specified in the grid file (``nz`` is not defined), then 3D fields
   are stored as Fourier components. In the Z dimension the coefficients
   must be stored as

   .. math::

      [n = 0, n = 1 (\textrm{real}), n = 1 (\textrm{imag}), n = 2
      (\textrm{real}), n = 2 (\textrm{imag}), \ldots ]

   where :math:`n` is the toroidal mode number. The size of the array
   must therefore be odd in the Z dimension, to contain a constant
   (:math:`n=0`) component followed by real/imaginary pairs for the
   non-axisymmetric components.

   If you are using IDL to create a grid file, there is a routine in
   ``tools/idllib/bout3dvar.pro`` for converting between BOUT++’s real
   and Fourier representation.

#. Real space, as values on grid points. If ``nz`` is set in the grid
   file, then 3D variables in the grid file must have size
   ``nx``\ :math:`\times`\ ``ny``\ :math:`\times`\ ``nz``. These are
   then read in directly into ``Field3D`` variables as required.

From EFIT files
---------------

An IDL code called “Hypnotoad” has been developed to create BOUT++ input
files from R-Z equilibria. This can read EFIT ’g’ files, find flux
surfaces, and calculate metric coefficients. The code is in
``tools/tokamak_grids/gridgen``, and has its own manual under the
``doc`` subdirectory.

From ELITE and GATO files
-------------------------

Currently conversions exist for ELITE ``.eqin`` and GATO ``dskgato``
equilibrium files. Conversion of these into BOUT++ input grids is in two
stages: In the first, both these input files are converted into a common
NetCDF format which describes the Grad-Shafranov equilibrium. These
intermediate files are then converted to BOUT++ grids using an
interactive IDL script.

Generating equilibria
---------------------

The directory ``tokamak_grids/shifted_circle`` contains IDL code to
generate shifted circle (large aspect ratio) Grad-Shafranov equilibria.

.. figure:: ../figs/grid_gen.*
    :alt: IDL routines and file formats used in taking output from
          different codes and converting into input to BOUT++.

    IDL routines and file formats used in taking output from different
    codes and converting into input to BOUT++.
