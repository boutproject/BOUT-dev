Mesh
====

The mesh is used in pretty much all parts of the code, and deals with
things like the geometry of the mesh (metric tensors etc.), and how
the mesh is divided between processors (communications). The `Mesh`
class defines an interface, and there is currently a single
implementation:

- `BoutMesh` (``src/mesh/boutmesh.cxx``) which is backwards compatible
  with the BOUT and BOUT-06 codes. This is a logically rectangular
  mesh so the number of radial points (x) can’t change in the
  poloidal direction (y).

Grid data sources
-----------------

All data sources inherit from `GridDataSource`. They must supply a
method to test if a variable exists, `GridDataSource::hasVar`::

    bool hasVar(const string &name);

and then use the `get <GridDataSource::get>` methods to get
integers or reals::

    bool get(Mesh *m, <type> &variable, const string &name);

Loading a mesh
--------------

The `Mesh constructor <Mesh::Mesh>` takes `GridDataSource` and
`Options` objects. You can also call `Mesh::create` with just one of
these objects, which will call out to the `MeshFactory` singleton to
create a mesh "automatically". This is the way that it is done in
:doc:`bout++.cxx <../_breathe_autogen/file/bout_09_09_8cxx>`. Once you
have instantiated a `Mesh` object, you can then call `Mesh::load` to
read in all the appropriate variables from the `GridDataSource`::

    mesh = Mesh::create();  ///< Create the mesh
    mesh->load();           ///< Load from sources. Required for Field initialisation

For post-processing of the results, it’s useful to have mesh
quantities in the dump files along with the results. To do this,
there’s the function `Mesh::outputVars` (see also `Datafile` and
`Options`)::

    // Create an output file from an Options object
    dump = Datafile(options->getSection("output"));

    // Possibly add some other variables to the output file
    ...

    // Save mesh configuration into output file
    mesh->outputVars(dump);

which is called during BOUT++ initialisation.

Implementation: BoutMesh
~~~~~~~~~~~~~~~~~~~~~~~~

`BoutMesh` class uses the BOUT indices (which trace back to UEDGE)::

    int ixseps1, ixseps2, jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2;

`ixseps1 <BoutMesh::ixseps1>` and `ixseps2 <BoutMesh::ixseps2>` give
the X location of the separatrices, and are equal in the case of
single-null configurations. The indexing is such that all points ``0
<= x < ixseps1`` are inside the separatrix, whilst ``ixseps1 <= x <
LocalNx`` are outside. See :ref:`sec-bout-topology` for more details.

Index ranges
------------

The `Mesh` class includes several public members which describe the size
of the mesh, and are used all over BOUT++ to loop over variables::

    /// Size of the mesh on this processor including guard/boundary cells
    int LocalNx, LocalNy, LocalNz;
    /// Local ranges of data (inclusive), excluding guard cells
    int xstart, xend, ystart, yend;

Getting data
------------

The `Mesh::load` code above needs to read data for the mesh, and
physics codes usually need to read their initial profiles during
initialisation.  To do this, Mesh provides an overloaded function
`Mesh::get`::

    int get(var, const char *name); // Request data from mesh file

where ``var`` can be just about any BOUT++ datatype (`Field2D`,
`Vector3D` etc.).

Implementation: BoutMesh
~~~~~~~~~~~~~~~~~~~~~~~~

For integers and BoutReals, the implementation is fairly trivial. Uses
the Mesh protected functions to find a data source and read data from
it::

    GridDataSource* s = findSource(name);  // Find a source of data
    s->open(name);                          // Open the source
    bool success = s->fetch(&ival, name);   // Get the data
    s->close();                             // Close the source

To read 2D and 3D fields, the branch-cuts need to be taken into account.

Communications
--------------

The most common type of communication is to just exchange all guard
cells with neighboring processors. Mesh provides the following commands
for doing this::

    template <typename... Ts>
    int communicate(Ts&... ts);      // Communicate one or more fields
    int communicate(FieldGroup);     // Communicate a group of fields
    comm_handle send(FieldGroup);    // Send data
    int wait(comm_handle);           // Receive data

`Mesh::communicate` can be used to communicate any number of variables
together, and makes the code quite clear. For example in
``examples/DriftInstability/2fluid.cxx`` around line 360::

    // Need to communicate jpar
    mesh->communicate(jpar);

Since this uses the `FieldData` interface like Datafile, this can be
used to communicate all BOUT++ field data types. You can also create a
`FieldGroup` object to group fields together, then communicate them
all together::

    FieldGroup comgrp;  // Group of variables for communication
    Field3D P;
    Vector3D V;

    comgrp.add(P); // Add the variables
    comgrp.add(V); // Usually done in PhysicsModel::init

    mesh->communicate(comgrp); // Communicate in PhysicsModel::rhs

Internally, this is how the templated `Mesh::communicate` works.

If you want to overlap communications with calculations then use the
`Mesh::send` and `Mesh::wait` functions instead of
`Mesh::communicate`::

    comm_handle ch = mesh->send(comgrp); // Start the communications
    // Calculations which don't need variables in comgrp
    wait(ch); // Wait for all communications to finish

Implementation: BoutMesh
~~~~~~~~~~~~~~~~~~~~~~~~

In `BoutMesh`, the communication is controlled by the variables::

    int UDATA_INDEST, UDATA_OUTDEST, UDATA_XSPLIT;
    int DDATA_INDEST, DDATA_OUTDEST, DDATA_XSPLIT;
    int IDATA_DEST, ODATA_DEST;
    int lower_inner_corner_dest, upper_inner_corner_dest, lower_outer_corner_dest,
        upper_outer_corner_dest; // destinations for the corner cells
    int lower_inner_corner_orig, upper_inner_corner_orig, lower_outer_corner_orig,
        upper_outer_corner_orig; // origins for the corner guard cells
    // y-limits of buffers communicated in x-direction. Include y-boundary cells but not
    // y-guard cells. Need different variables for sending and receiving because y-boundary
    // might be present on sending proc but not receiving proc or vice versa
    int IDATA_buff_lowerY_send, IDATA_buff_upperY_send, ODATA_buff_lowerY_send,
        ODATA_buff_upperY_send;
    int IDATA_buff_lowerY_recv, IDATA_buff_upperY_recv, ODATA_buff_lowerY_recv,
        ODATA_buff_upperY_recv;
    // x-limits of buffers communicated in y-direction. Include x-boundary cells but not
    // x-guard cells.
    int YDATA_buff_innerX, YDATA_buff_outerX;

In the Y direction, each boundary region (**U**\ p and **D**\ own in Y)
can be split into two, with ``x < UDATA_XSPLIT`` going to the
processor index ``UDATA_INDEST``, and ``UDATA_INDEST <= x`` going
to ``UDATA_OUTDEST``. Similarly for the Down boundary. Since there are
no branch-cuts in the X direction, there is just one destination for the
**I**\ nner and **O**\ uter boundaries. In all cases a negative
processor number means that there’s a domain boundary.

The corners (cells that are both x-guards and y-guards) are handled specially.  Away from
the boundaries the MXG*MYG cells at the corner of the processor's grid are sent to
``*_corner_dest``, and the corner cells are received from ``*_corner_orig``; the two may
not be the same at X-points, if the separatrix location it at a processor boundary. The
sending and receiving locations are set to be consistent with the behaviour if the
separatrix is in the interior of a processor's grid, which means communication as if all
guard cells were first communicated in y, and then communicated in x (this is not actually
done so that all the communications can be done at the same time, reducing latency).

Where there is a physical boundary (where boundary conditions are applied), the boundary
cells should be communicated to fill the corner cells. Since the communication pattern is
different from the one for the corner cells in the interior, the simplest implementation
is:
* add x-boundary cells to the y-communications
* add y-boundary cells to the x-communications
* set the ``*_corner_dest`` and ``*_corner_orig`` corresponding to the boundary to -1 so
  that the 'corner communications' are not used for boundaries.

If the option ``mesh::include_corner_cells`` is set to ``false`` (default is ``true``),
then the previous behaviour (up to BOUT++ v4) is restored:
* In the y-direction
     * [0,UDATA_XSPLIT) is communicated to UDATA_INDEST at the upper boundary
     * [UDATA_XSPLIT, LocalNy) is communicated to UDATA_OUTDEST at the upper boundary
     * [0,DDATA_XSPLIT) is communicated to DDATA_INDEST at the lower boundary
     * [DDATA_XSPLIT, LocalNy) is communicated to DDATA_OUTDEST at the lower boundary
* in the x-direction [MYG, MYG + MYSUB) is communicated at both boundaries

X communications
----------------

For parallel Laplacian inversions, communication is needed in the X
direction only, and involves quantities which are not in Fields::

    bool firstX();  // True if at the inner X boundary
    bool lastX();   // True if at the outer X boundary
    int NXPE, PE_XIND; // Number of processors in X, and X processor index
    int sendXOut(BoutReal *buffer, int size, int tag);
    sendXIn(BoutReal *buffer, int size, int tag);
    comm_handle irecvXOut(BoutReal *buffer, int size, int tag);
    comm_handle irecvXIn(BoutReal *buffer, int size, int tag);

The variables `Mesh::NXPE` and `Mesh::PE_XIND` shouldn’t really be
there, but are currently needed because the SPT algorithm in
`LaplaceSPT` needs to know when it’s going to be next and so keep
track of which processor number is currently working. This logic to
pass a problem along a chain in X should really be moved into Mesh.

Y-Z surface communications
--------------------------

Some operations (like parallel inversions in
``bout++/src/invert/invert_parderiv.cxx``) need to be performed on Y-Z
surfaces, i.e. slices at constant X. This needs to be able to handle
open and closed surfaces, and that closed surfaces may need a shift in
the Z direction to match one end onto the other (a twist-shift
condition).

The simplest operation is to average a quantity over Y with
`averageY`.

To test if a particular surface is closed, there is the function
`periodicY`.

The most general way to access data on surfaces is to use the
`SurfaceIter` iterator, which can be created using
`SurfaceIter::SurfaceIter`::

    SurfaceIter* surface(mesh);

This then allows looping over the surfaces in the usual way::

    for(surf->first(); !surf->isDone(); surf->next()) {
      ...
    }

To test if the surface is closed, there’s the test `SurfaceIter::closed`::

    bool surf->closed(BoutReal &ts)

which returns true if the surface is closed, along with the twist-shift
angle.

Initial profiles
----------------

The initial profiles code needs to construct a solution which is
smooth everywhere, with a form of perturbation specified in the input
file for each direction. In order to do this, it needs a continuous
function to use as an index. This is supplied by the functions
`Mesh::GlobalX` and `Mesh::GlobalY`::

    BoutReal GlobalX(int jx); // Continuous X index between 0 and 1
    BoutReal GlobalY(int jy); // Continuous Y index (0 -> 1)

which take a local x or y index and return a globally continuous x or y
index.

Differencing
------------

The mesh spacing is given by the public members `Mesh::dx`, `Mesh::dy`
and `Mesh::dx`::

    // These used for differential operators
    Field2D dx, dy;
    Field2D d2x, d2y;    // 2nd-order correction for non-uniform meshes
    BoutReal zlength, dz;    // Derived from options (in radians)

Metrics
-------

While `Mesh` handles the numerical details of the mesh, the "physical"
details are handled by `Coordinates`. The contravariant and covariant
metric tensor components are public members of `Coordinates`::

    // Contravariant metric tensor (g^{ij})
    Field2D g11, g22, g33, g12, g13, g23; // These are read in grid.cxx

    // Covariant metric tensor
    Field2D g_11, g_22, g_33, g_12, g_13, g_23;

    int calcCovariant();     // Invert contravatiant metric to get covariant
    int calcContravariant(); // Invert covariant metric to get contravariant

If only one of these sets is modified by an external code, then
`Coordinates::calcCovariant` and `Coordinates::calcContravariant` can
be used to calculate the other (uses Gauss-Jordan currently).

From the metric tensor components, `Coordinates` calculates several
other useful quantities::

    int jacobian(); // Calculate J and Bxy
    Field2D J; // Jacobian
    Field2D Bxy; // Magnitude of B = nabla z times nabla x

    /// Calculate differential geometry quantities from the metric tensor
    int geometry();

    // Christoffel symbol of the second kind (connection coefficients)
    Field2D G1_11, G1_22, G1_33, G1_12, G1_13;
    Field2D G2_11, G2_22, G2_33, G2_12, G2_23;
    Field2D G3_11, G3_22, G3_33, G3_13, G3_23;

    Field2D G1, G2, G3;

These quantities are public and accessible everywhere, but this is
because they are needed in a lot of the code. They shouldn’t change
after initialisation, unless the physics model starts doing fancy things
with deforming meshes.

Miscellaneous
-------------

There are some public members of `Mesh` which are there for some
specific task and don’t really go anywhere else (yet).

To perform radial derivatives in tokamak geometry, interpolation is
needed in the Z direction. This is done by shifting in Z by a phase
factor, performing the derivatives, then shifting back. The following
public variables are currently used for this::

    bool ShiftXderivs; // Use shifted X derivatives
    int  ShiftOrder;   // Order of shifted X derivative interpolation
    Field2D zShift;    // Z shift for each point (radians)

    Field2D ShiftTorsion; // d <pitch angle> / dx. Needed for vector differentials (Curl)
    Field2D IntShiftTorsion; // Integrated shear (I in BOUT notation)
    bool IncIntShear; // Include integrated shear (if shifting X)

    int  TwistOrder;   // Order of twist-shift interpolation

This determines what order method to use for the interpolation at the
twist-shift location, with ``0`` meaning FFT during communication. Since
this must be 0 at the moment it’s fairly redundant and should be
removed.

A (currently experimental) feature is::

    bool StaggerGrids;    ///< Enable staggered grids (Centre, Lower). Otherwise all vars are cell centred (default).
