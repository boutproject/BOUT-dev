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
(which can be used in physics models) in ``bout++/include/bout``, and
source code and low-level includes in ``bout++/src``. Many parts of
the code are defined by their interface, and can have multiple
different implementations. An example is the time-integration solvers:
many different implementations are available, some of which use
external libraries, but all have the same interface and can be used
interchangeably. This is reflected in the directory structure inside
``bout++/src``. A common pattern is to store individual
implementations of an interface in a subdirectory called ``impls``.

::

    include/bout/foo.hxx
    src/.../foo.cxx
    src/.../impls/one/one.hxx
    src/.../impls/one/one.cxx

where ``foo.hxx`` defines the interface, ``foo.cxx`` implements common
functions used in several implementations. Individual implementations
are stored in their own subdirectories of ``impls``. Components which
follow this pattern include ``invert/laplace`` and ``invert/parderiv``
inversions, ``mesh``, and ``solver``.

The layout of the ``src/`` directory is as follows:

- ``src/bout++.cxx``: Main file which initialises, runs and finalises
  BOUT++. The two most important public functions are `BoutInitialise`
  and `BoutFinalise` for starting/stopping the library.

- ``src/field``

  - Implementations of "fields" (the scalars `Field2D`, `Field3D`,
    `FieldPerp`, and vectors `Vector2D`, `Vector3D`), as well as basic
    operations (arithmetic, :ref:`sec-algebraic-ops`, and vector
    calculus), and :ref:`initialisation <sec-variable-init>`.

- ``src/invert``

  - Implementations of different inverse operators, including Fourier
    transforms (via an interface to `FFTW <http://www.fftw.org>`_, see
    `bout::fft::rfft` and `bout::fft::irfft`).

- ``src/invert/laplace``

  - Implementations of some inverse generalised Laplacian operator
    variations, where some directions may be constant. See
    :ref:`sec-laplacian` for more details.

- ``src/invert/parderiv``

  - Inversion of parallel derivatives, intended for use in
    preconditioners. See `InvertPar`

- ``src/invert/pardiv``

  - Inversion of parallel divergence, intended for use in
    :ref:`sec-nonlocal-heatflux`.

- ``src/mesh``

  - Implementations of low-level numerical routines. This includes
    things like the :ref:`inter-process communications
    <sec-communications>`, :ref:`boundary conditions <sec-bndryopts>`,
    :ref:`derivative operators <sec-diffops>`, interpolation, the
    coordinate system (including :ref:`sec-parallel-transforms`), and
    the :ref:`sec-mesh` itself.

- ``src/physics``

  - This contains some specialised physics operators and routines,
    such as gyro-averaging and :ref:`sec-nonlocal-heatflux`.

- ``src/solver``

  - Implementations of :ref:`time integration solvers
    <sec-time-integration>`

- ``src/sys``

  - General purpose utilities used throughout the library, such as
    `BoutException`, wrappers for C libraries like ``PETSc`` and
    ``HYPRE``, screen and file input and output.

