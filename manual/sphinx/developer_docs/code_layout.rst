Code layout
===========

BOUT++ is organised into classes and groups of functions which operate
on them: It’s not purely object-oriented, but takes advantage of many of
C++’s object-oriented features.

:numref:`fig-layout1` shows the most important parts of BOUT++ and how
they fit together.

.. _fig-layout1:
.. figure:: ../figs/layout1.*
   :alt: Overview of BOUT++ control flow

   Overview of BOUT++ control flow during initialisation (red), and
   running (blue)

The initialisation process is shown in red: basic information is first
read from the grid file (e.g. size of the grid, topology etc.), then
the user-supplied initialisation code is called. This code can read
other variables from the grid, and makes at least one call to
`PhysicsModel::bout_solve` to specify a variable to be evolved. The
main thing `bout_solve <PhysicsModel::bout_solve>` does is to add
these variables to the solver.

The process of running a timestep is shown in blue in
:numref:`fig-layout1`: The main loop calls the solver, which in turn
calls PVODE. To evolve the system PVODE makes calls to the RHS
function inside solver. This moves data between PVODE and BOUT++, and
calls the user-supplied `PhysicsModel::rhs` code to calculate
time-derivatives. Much of the work calculating time-derivatives
involves differential operators.

Calculation of the `RHS function <PhysicsModel::rhs>`, and handling of
data in BOUT++ involves many different
components. :numref:`fig-layout2` shows (most) of the classes and
functions involved, and the relationships between them. Some thought
was put into how this should be organised, but it has also changed
over time, so some parts could be cleaner.

.. _fig-layout2:
.. figure:: ../figs/layout2.*
   :alt: Relationships used in calculating the RHS function

   Relationship between important classes and functions used in
   calculating the RHS function

Directories
-----------

The source code for the core of BOUT++ is divided into include files
(which can be used in physics models) in ``bout++/include``, and source
code and low-level includes in ``bout++/src``. Many parts of the code
are defined by their interface, and can have multiple different
implementations. An example is the time-integration solvers: many
different implementations are available, some of which use external
libraries, but all have the same interface and can be used
interchangeably. This is reflected in the directory structure inside
``bout++/src``. A common pattern is to store individual implementations
of an interface in a subdirectory called ``impls``.

::

    include/foo.hxx
    src/.../foo.cxx
    src/.../foo_factory.hxx
    src/.../foo_factory.cxx
    src/.../impls/one/one.hxx
    src/.../impls/one/one.cxx

where ``foo.hxx`` defines the interface, ``foo.cxx`` implements common
functions used in several implementations. ``foo_factory`` creates new
implementations, and is the only file which includes all the
implementations. Individual implementations are stored in their own
subdirectories of ``impls``. Components which follow this pattern
include ``fileio`` formats, ``invert/laplace`` and ``invert/parderiv``
inversion codes, ``mesh``, and ``solver``.

The current source code files are:

- :doc:`bout++.cxx<../_breathe_autogen/file/bout_09_09_8cxx>`: Main
  file which initialises, runs and finalises BOUT++. Currently
  contains a `main()` function, though this is being removed shortly.

- field

   - :doc:`field2d.cxx<../_breathe_autogen/file/field2d_8cxx>`
     implements the `Field2D` class. This is a scalar field which
     varies only in :math:`x` and :math:`y` and is used for things
     like metric tensor components and initial profiles. It supplies
     lots of overloaded operators and functions on these objects.

   - :doc:`field3d.cxx<../_breathe_autogen/file/field3d_8cxx>`
     implements the `Field3D` class, which varies in :math:`x`,
     :math:`y` and :math:`z`. Since these handle a lot more memory
     than Field2D objects, the memory management is more complicated
     and includes reference counting. See section
     :ref:`sec-memorymanage` for more details.

   - :doc:`field_data.cxx<../_breathe_autogen/file/field__data_8cxx>`
     Implements some functions in the `FieldData` class. This is a
     mainly pure virtual interface class which is inherited by
     `Field2D` and `Field3D`.

   - :doc:`fieldperp.cxx<../_breathe_autogen/file/fieldperp_8cxx>`
     implements a `FieldPerp` class to store slices perpendicular to
     the magnetic field i.e. they are a function of :math:`x` and
     :math:`z` only. This is mainly used for Laplacian inversion
     routines, and needs to be integrated with the other fields
     better.

   - :doc:`initialprofiles.cxx<../_breathe_autogen/file/initialprofiles_8cxx>`
     routines to set the initial values of fields when a simulation
     first starts. Reads settings from the option file based on the name
     of the variable.

   - :doc:`vecops.cxx<../_breathe_autogen/file/vecops_8cxx>` a
     collection of function to operate on vectors.  Contains things
     like ``Grad``, ``Div`` and ``Curl``, and uses a combination of
     field differential operators (in
     :doc:`difops.cxx<../_breathe_autogen/file/difops_8cxx>`) and
     metric tensor components (in `Mesh`).

   - :doc:`vector2d.cxx<../_breathe_autogen/file/vector2d_8cxx>`
     implements the `Vector2D` class, which uses a `Field2D` object
     for each of its 3 components. Overloads operators to supply
     things like dot and cross products.

   - :doc:`vector3d.cxx<../_breathe_autogen/file/vector3d_8cxx>`
     implements `Vector3D` by using a `Field3D`
     object for each component.

   - :doc:`where.cxx<../_breathe_autogen/file/where_8cxx>` supplies
     functions for choosing between values based on selection
     criteria.

- fileio

   - :doc:`datafile.cxx<../_breathe_autogen/file/datafile_8cxx>`
     supplies an abstract `DataFile` interface for data
     input and output. Handles the conversion of data in fields and
     vectors into blocks of data which are then sent to a specific
     file format.

   - :doc:`formatfactory.cxx<../_breathe_autogen/file/formatfactory_8cxx>`

   - :doc:`formatfactory.hxx<../_breathe_autogen/file/formatfactory_8hxx>`

   - impls

      - :doc:`emptyformat.hxx<../_breathe_autogen/file/emptyformat_8hxx>`

      - hdf5

         - :doc:`h5_format.cxx<../_breathe_autogen/file/h5__format_8cxx>` implements an
           interface to the HDF5 library

         - :doc:`h5_format.hxx<../_breathe_autogen/file/h5__format_8hxx>`

      - netcdf

         - :doc:`nc_format.cxx<../_breathe_autogen/file/nc__format_8cxx>` implements an
           interface to the NetCDF-4 library

         - :doc:`nc_format.hxx<../_breathe_autogen/file/nc__format_8hxx>`

      - netcdf4

         - :doc:`ncxx<../_breathe_autogen/file/ncxx4_8cxx>`
           implements an interface to the NetCDF-4 library using the
           C++ API

         - :doc:`ncxx<../_breathe_autogen/file/ncxx4_8hxx>`

      - pnetcdf

         - :doc:`pnetcdf.cxx<../_breathe_autogen/file/pnetcdf_8cxx>`
           Parallel NetCDF interface

         - :doc:`pnetcdf.hxx<../_breathe_autogen/file/pnetcdf_8hxx>`

- invert

   - :doc:`fft_fftw.cxx<../_breathe_autogen/file/fft__fftw_8cxx>`
     implements the :doc:`fft.hxx<../_breathe_autogen/file/fft_8hxx>`
     interface by calling the Fastest Fourier Transform in the West
     (FFTW) library.

- invert / laplace

   - :doc:`invert_laplace.cxx<../_breathe_autogen/file/invert__laplace_8cxx>` uses Fourier
      decomposition in :math:`z` combined with tri- and band-diagonal
      solvers in :math:`x` to solve Laplacian problems.

   - impls

      - serial\_tri

         - :doc:`serial_tri.hxx<../_breathe_autogen/file/serial__tri_8hxx>`

         - :doc:`serial_tri.cxx<../_breathe_autogen/file/serial__tri_8cxx>`

      - serial\_band

         - :doc:`serial_band.hxx<../_breathe_autogen/file/serial__band_8hxx>`

         - :doc:`serial_band.cxx<../_breathe_autogen/file/serial__band_8cxx>`

      - spt

         - :doc:`spt.hxx<../_breathe_autogen/file/spt_8hxx>`

         - :doc:`spt.cxx<../_breathe_autogen/file/spt_8cxx>`

      - pdd

         - :doc:`pdd.hxx<../_breathe_autogen/file/pdd_8hxx>`

         - :doc:`pdd.cxx<../_breathe_autogen/file/pdd_8cxx>`

- invert / parderiv

   -
     :doc:`invert_parderiv.cxx<../_breathe_autogen/file/invert__parderiv_8cxx>`
     inverts a problem involving only parallel :math:`y`
     derivatives. Intended for use in some preconditioners.

   - impls

      - cyclic

         - :doc:`cyclic.cxx<../_breathe_autogen/file/cyclic_8cxx>`

         - :doc:`cyclic.hxx<../_breathe_autogen/file/cyclic_8hxx>`

- :doc:`lapack_routines.cxx<../_breathe_autogen/file/lapack__routines_8cxx>` supplies an
   interface to the LAPACK linear solvers, which are used by the
   ``invert_laplace`` routines.

- mesh

   - :doc:`boundary_factory.cxx<../_breathe_autogen/file/boundary__factory_8cxx>` creates boundary
     condition operators which can then be applied to
     fields. Described in section :ref:`sec-BoundaryFactory`.

   - :doc:`boundary_region.cxx<../_breathe_autogen/file/boundary__region_8cxx>` implements a way
     to describe and iterate over boundary regions. Created by the
     mesh, and then used by boundary conditions. See
     section :ref:`sec-BoundaryRegion` for more details.

   - :doc:`boundary_standard.cxx<../_breathe_autogen/file/boundary__standard_8cxx>` implements some
     standard boundary operations and modifiers such as ``Neumann``
     and ``Dirichlet``.

   - :doc:`difops.cxx<../_breathe_autogen/file/difops_8cxx>` is a
     collection of differential operators on scalar fields. It uses
     the differential methods in :doc:`derivs.cxx<../_breathe_autogen/file/derivs_8cxx>` and the metric tensor
     components in `Mesh` to compute operators.

   - :doc:`interpolation.cxx<../_breathe_autogen/file/interpolation_8cxx>` contains functions
     for interpolating fields

   - :doc:`mesh.cxx<../_breathe_autogen/file/mesh_8cxx>` is the base
     class for the `Mesh` object. Contains routines useful
     for all `Mesh` implementations.

   - impls

      - bout

         - :doc:`boutmesh.cxx<../_breathe_autogen/file/boutmesh_8cxx>`
           implements a mesh interface which is compatible with BOUT
           grid files.

         - :doc:`boutmesh.hxx<../_breathe_autogen/file/boutmesh_8hxx>`

- physics

   - :doc:`gyro_average.cxx<../_breathe_autogen/file/gyro__average_8cxx>`
      gyro-averaging operators

   - :doc:`smoothing.cxx<../_breathe_autogen/file/smoothing_8cxx>`
     provides smoothing routines on scalar fields

   - :doc:`sourcex.cxx<../_breathe_autogen/file/sourcex_8cxx>` contains
     some useful routines for creating sources and sinks in physics
     equations.

- solver

   - :doc:`solver.cxx<../_breathe_autogen/file/solver_8cxx>` is the
     interface for all solvers

   - impls

      - cvode

         - :doc:`cvode.cxx<../_breathe_autogen/file/cvode_8cxx>` is the
           implementation of `Solver` which interfaces with
           the SUNDIALS CVODE library.

         - :doc:`cvode.hxx<../_breathe_autogen/file/cvode_8hxx>`

      - ida

         - :doc:`ida.cxx<../_breathe_autogen/file/ida_8cxx>` is the
           implementation which interfaces with the SUNDIALS IDA
           library

         - :doc:`ida.hxx<../_breathe_autogen/file/ida_8hxx>`

      - petsc

         - :doc:`petsc.cxx<../_breathe_autogen/file/petsc_8cxx>` is the
           interface to the PETSc time integration routines

         - :doc:`petsc.hxx<../_breathe_autogen/file/petsc_8hxx>`

      - pvode

         - :doc:`pvode.cxx<../_breathe_autogen/file/pvode_8cxx>`
           interfaces with the 1998 (pre-SUNDIALS) version of PVODE
           (which became CVODE).

         - :doc:`pvode.hxx<../_breathe_autogen/file/pvode_8hxx>`

- sys

   - :doc:`boutcomm.cxx<../_breathe_autogen/file/boutcomm_8cxx>`

   - :doc:`boutexception.cxx<../_breathe_autogen/file/boutexception_8cxx>`
     is an exception class which are used for error handling

   - :doc:`derivs.cxx<../_breathe_autogen/file/derivs_8cxx>` contains
     basic derivative methods such as upwinding, central difference
     and WENO methods. These are then used by
     :doc:`difops.cxx<../_breathe_autogen/file/difops_8cxx>`. Details are
     given in section :ref:`sec-derivatives`.

   - :doc:`msg_stack.cxx<../_breathe_autogen/file/msg__stack_8cxx>` is
     part of the error handling system. It maintains a stack of
     messages which can be pushed onto the stack at the start of a
     function, then removed (popped) at the end. If an error occurs or
     a segmentation fault is caught then this stack is printed out and
     can help to find errors.

   - :doc:`options.cxx<../_breathe_autogen/file/options_8cxx>` provides
     an interface to the BOUT.inp option file and the command-line
     options.

   - :doc:`optionsreader.cxx<../_breathe_autogen/file/optionsreader_8cxx>`

   - :doc:`output.cxx<../_breathe_autogen/file/output_8cxx>`

   - :doc:`range.cxx<../_breathe_autogen/file/range_8cxx>` Provides the
     RangeIterator class, used to iterate over a set of
     ranges. Described in section :ref:`sec-rangeiterator`

   - :doc:`timer.cxx<../_breathe_autogen/file/timer_8cxx>` a class for
     timing parts of the code like communications and file
     I/O. Described in section :ref:`sec-timerclass`

   - :doc:`utils.cxx<../_breathe_autogen/file/utils_8cxx>` contains
     miscellaneous small useful routines such as allocating and
     freeing arrays.

   - options

      - :doc:`optionparser.hxx<../_breathe_autogen/file/optionparser_8hxx>`

      - :doc:`options_ini.cxx<../_breathe_autogen/file/options__ini_8cxx>`

      - :doc:`options_ini.hxx<../_breathe_autogen/file/options__ini_8hxx>`

