Mesh
====

The mesh is used in pretty much all parts of the code, and deals with
things like the geometry of the mesh (metric tensors etc.), and how
the mesh is divided between processors (communications). The
:cpp:class:`Mesh` class defines an interface, and there is currently a
single implementation:

- :cpp:class:`BoutMesh` (``src/mesh/boutmesh.cxx``) which is backwards
   compatible with the BOUT and BOUT-06 codes. This is a logically
   rectangular mesh so the number of radial points (x) can’t change in
   the poloidal direction (y).

Grid data sources
-----------------

All data sources inherit from :cpp:class:`GridDataSource`, defined in
:doc:`grid.hxx<../_breathe_autogen/file/griddata_8hxx>` at
line 43. They must supply a method to test if a variable exists::

    bool GridDataSource::hasVar(const char *name);

a method to get the size of the variable

::

    vector<int> GridDataSource::getSize(const char *name);

To fetch data, first the (x,y,z) origin must be set::

    bool GridDataSource::setOrigin(int x = 0, int y = 0, int z = 0);

and then use methods to fetch integers or reals::

    bool GridDataSource::fetch(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0);
    bool GridDataSource::fetch(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0);

In addition, GridDataSource implementations can have methods which
should be called before and after variables are accessed::

    void GridDataSource::open(const char *name = NULL);
    void GridDataSource::close();

Loading a mesh
--------------

To load in a mesh from a file or other source, there are the commands::

    int addSource(GridDataSource);   // Add a data source
    int load();                      // Load from added data sources
    int load(GridDataSource);        // Load from specified data source

all of which return an error code (0 if successful). ``addSource`` is
used to add a set of input data sources which inherit from
:cpp:class:`GridDataSource`. ``load()`` loads the mesh from these
sources, querying each data source in turn for the required variables
(in the order in which they were added). ``load(GridDataSource)``
loads the mesh from only the supplied data source.

In :doc:`bout++.cxx<../_breathe_autogen/file/bout_09_09_8cxx>`, this
is used to initialise the mesh::

    mesh->addSource(new GridFile(data_format(grid_name), grid_name));
    if(mesh->load()) {
      output << "Failed to read grid. Aborting\n";
      return 1;
    }

which creates a :cpp:class:`GridFile` object based on the data format
of the grid file name, then adds that as a source of data for Mesh.

For post-processing of the results, it’s useful to have mesh quantities
in the dump files along with the results. To do this, there’s the
function

::

    void outputVars(Datafile &file); // Add mesh vars to file

which is called during BOUT++ initialisation.

Implementation: BoutMesh
~~~~~~~~~~~~~~~~~~~~~~~~

BoutMesh class uses the BOUT indices (which trace back to UEDGE)::

    int ixseps1, ixseps2, jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2;

``ixseps1`` and ``ixseps2`` give the X location of the separatrices, and
are equal in the case of single-null configurations. The indexing is
such that all points ``0 <= x < ixseps1`` are inside the separatrix,
whilst ``ixseps1 <= x < LocalNx`` are outside.

Index ranges
------------

The Mesh class includes several public members which describe the size
of the mesh, and are used all over BOUT++ to loop over variables::

    /// Size of the mesh on this processor including guard/boundary cells
    int LocalNx, LocalNy, LocalNz;
    /// Local ranges of data (inclusive), excluding guard cells
    int xstart, xend, ystart, yend;

Getting data
------------

The ``load()`` code above needs to read data for the mesh, and physics
codes usually need to read their initial profiles during initialisation.
To do this, Mesh provides an overloaded function ``get``::

    int get(var, const char *name); // Request data from mesh file

where ``var`` can be just about any BOUT++ datatype
(:cpp:class:`Field2D`, :cpp:class:`Vector3D` etc.).

Implementation: BoutMesh
~~~~~~~~~~~~~~~~~~~~~~~~

For integers and BoutReals, the implementation is fairly trivial. Uses
the Mesh protected functions to find a data source and read data from
it.

::

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

    int communicate(FieldData, ...); // Communicate one or more fields
    int communicate(FieldGroup);     // Communicate a group of fields
    int communicate(FieldData);      // Returns error code
    comm_handle send(FieldGroup);    // Send data
    int wait(comm_handle);           // Receive data

``communicate(FieldData)`` can (currently) be used to communicate up to
4 variables together, and makes the code quite clear. For example in
``examples/DriftInstability/2fluid.cxx`` around line 360::

    // Need to communicate jpar
    mesh->communicate(jpar);

Since this uses the :cpp:class:`FieldData` interface like Datafile,
this can be used to communicate all BOUT++ field data types. The limit
of 4 is because the C-style ``varargs`` system doesn’t work with “non
POD” variables, i.e. classes. To communicate a larger number of
variables, create a :cpp:class:`FieldGroup` object to group fields
together, then communicate them all together::

    FieldGroup comgrp;  // Group of variables for communication
    Field3D P;
    Vector3D V;

    comgrp.add(P); // Add the variables
    comgrp.add(V); // Usually done in physics_init

    mesh->communicate(comgrp); // Communicate in physics_run

If you want to overlap communications with calculations then use the
``send`` and ``wait`` functions instead of ``communicate``.

::

    comm_handle ch = mesh->send(comgrp); // Start the communications
    // Calculations which don't need variables in comgrp
    wait(ch); // Wait for all communications to finish

Implementation: BoutMesh
~~~~~~~~~~~~~~~~~~~~~~~~

In BoutMesh, the communication is controlled by the variables

::

    int UDATA_INDEST, UDATA_OUTDEST, UDATA_XSPLIT;
    int DDATA_INDEST, DDATA_OUTDEST, DDATA_XSPLIT;
    int IDATA_DEST, ODATA_DEST;

In the Y direction, each boundary region (**U**\ p and **D**\ own in Y)
can be split into two, with ``0 <= x < UDATA_XSPLIT`` going to the
processor index ``UDATA_INDEST``, and ``UDATA_INDEST <= x < LocalNx`` going
to ``UDATA_OUTDEST``. Similarly for the Down boundary. Since there are
no branch-cuts in the X direction, there is just one destination for the
**I**\ nner and **O**\ uter boundaries. In all cases a negative
processor number means that there’s a domain boundary.

X communications
----------------

For parallel Laplacian inversions, communication is needed in the X
direction only, and involves quantities which are not in Fields.

::

    bool firstX();  // True if at the inner X boundary
    bool lastX();   // True if at the outer X boundary
    int NXPE, PE_XIND; // Number of processors in X, and X processor index
    int sendXOut(BoutReal *buffer, int size, int tag);
    sendXIn(BoutReal *buffer, int size, int tag);
    comm_handle irecvXOut(BoutReal *buffer, int size, int tag);
    comm_handle irecvXIn(BoutReal *buffer, int size, int tag);

The variables ``NXPE`` and ``PE_XIND`` shouldn’t really be there, but
are currently needed because the SPT algorithm in :doc:`invert_laplace.cxx<../_breathe_autogen/file/invert__laplace_8cxx>`
needs to know when it’s going to be next and so keep track of which
processor number is currently working. This logic to pass a problem
along a chain in X should really be moved into Mesh.

Y-Z surface communications
--------------------------

Some operations (like parallel inversions in
``bout++/src/invert/invert_parderiv.cxx``) need to be performed on Y-Z
surfaces, i.e. slices at constant X. This needs to be able to handle
open and closed surfaces, and that closed surfaces may need a shift in
the Z direction to match one end onto the other (a twist-shift
condition).

The simplest operation is to average a quantity over Y::

    const Field2D averageY(const Field2D &f); // Average in Y

Currently this is only implemented for 2D fields. More generally a set
of FieldData objects could be used.

To test if a particular surface is closed, there is the function

::

    bool surfaceClosed(int jx, BoutReal &ts); // Test if a surface is closed, and if so get the twist-shift angle

The most general way to access data on surfaces is to use an iterator,
which can be created using::

    SurfaceIter* iterateSurfaces();

This then allows looping over the surfaces in the usual way

::

    for(surf->first(); !surf->isDone(); surf->next()) {
      ...
    }

**NB**: This iterator splits the surfaces between processors, so each
individual processor will iterate over a different set of surfaces. This
is to allow automatic load balancing when gathering and scattering data
from an entire surface onto one processor using::

    surf->gather(FieldData, BoutReal *recvbuffer);
    surf->scatter(BoutReal *sendbuffer, Field result);

The buffer is assumed to be large enough to hold all the data. To get
the number of points in Y for this surface, use

::

    int ysize = surf->ysize();

To test if the surface is closed, there’s the test

::

    bool surf->closed(BoutReal &ts)

which returns true if the surface is closed, along with the twist-shift
angle.

Initial profiles
----------------

The initial profiles code needs to construct a solution which is smooth
everywhere, with a form of perturbation specified in the input file for
each direction. In order to do this, it needs a continuous function to
use as an index. This is supplied by the functions::

    BoutReal GlobalX(int jx); // Continuous X index between 0 and 1
    BoutReal GlobalY(int jy); // Continuous Y index (0 -> 1)

which take a local x or y index and return a globally continuous x or y
index.

Differencing
------------

The mesh spacing is given by the public members

::

    // These used for differential operators
    Field2D dx, dy;
    Field2D d2x, d2y;    // 2nd-order correction for non-uniform meshes
    BoutReal zlength, dz;    // Derived from options (in radians)

Metrics
-------

The contravariant and covariant metric tensor components are public
members of :cpp:class:`Mesh`::

    // Contravariant metric tensor (g^{ij})
    Field2D g11, g22, g33, g12, g13, g23; // These are read in grid.cxx

    // Covariant metric tensor
    Field2D g_11, g_22, g_33, g_12, g_13, g_23;

    int calcCovariant();     // Invert contravatiant metric to get covariant
    int calcContravariant(); // Invert covariant metric to get contravariant

If only one of these sets is modified by an external code, then
``calc_covariant`` and ``calc_contravariant`` can be used to calculate
the other (uses Gauss-Jordan currently).

From the metric tensor components, Mesh calculates several other useful
quantities::

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

There are some public members of Mesh which are there for some specific
task and don’t really go anywhere else (yet).

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

::

    int  TwistOrder;   // Order of twist-shift interpolation

This determines what order method to use for the interpolation at the
twist-shift location, with ``0`` meaning FFT during communication. Since
this must be 0 at the moment it’s fairly redundant and should be
removed.

A (currently experimental) feature is::

    bool StaggerGrids;    ///< Enable staggered grids (Centre, Lower). Otherwise all vars are cell centred (default).
