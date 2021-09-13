.. _sec-gpusupport:

GPU support
===========

This section describes work in progress to develop GPU support in BOUT++ models.
It includes both configuration and compilation on GPU systems, but also ways to
write physics models which are designed to give higher performance. These methods
may also be beneficial for CPU architectures.

CMake configuration
-------------------

To compile BOUT++ components into GPU kernels a few different pieces need to be configured to work together:
RAJA, Umpire, and a CUDA compiler.


.. _tab-gpusupport-cmake:
.. table:: Useful CMake configuration settings

   +----------------------+-----------------------------------------+------------------------+
   | Name                 | Description                             | Default                |
   +======================+=========================================+========================+
   | BOUT_ENABLE_RAJA     | Include RAJA header, use RAJA loops     | Off                    |
   +----------------------+-----------------------------------------+------------------------+
   | BOUT_ENABLE_UMPIRE   | Use Umpire to allocate Array memory     | Off                    |
   +----------------------+-----------------------------------------+------------------------+
   | BOUT_ENABLE_CUDA     | Compile with nvcc compiler              | Off                    |
   +----------------------+-----------------------------------------+------------------------+
   | CUDA_ARCH            | Set CUDA architecture to compile for    | compute_70,code=sm_70  |
   +----------------------+-----------------------------------------+------------------------+
   | BOUT_ENABLE_WARNINGS | nvcc has incompatible warning flags     | On (turn Off for CUDA) |
   +----------------------+-----------------------------------------+------------------------+


Single index operators
----------------------

BOUT++ models are typically implemented using operations which take in
fields, perform an operation, and return a new field. These are
convenient, but the consequence is that an expression like
``Grad_par(phi) * B0`` contains two loops over the domain, one for the
gradient operator ``Grad_par``, and another for the
multiplication. Complex models can contain dozens of these loops. When
using OpenMP or GPU threads this results in many small kernels being
launched, and typically poor efficiency.

The "single index operators" in BOUT++ offer a way to manually combine
loops over the domain into a smaller number of loops. It is perhaps
less convenient than a template expression system might be, but
considerably easier to debug and maintain.

Single index operators have the same name as field operations, but the interface
has two key differences:

1. The functions take an index as an additional final argument
2. Rather than fields (e.g Field2D, Field3D types), these operate on
   field accessors (Field2DAccessor, FieldAccessor types). These offer
   faster, more direct, but less safe access to the underlying data
   arrays.

For example a simple field operation::

  Field3D phi;
  ...
  Field3D result = DDX(phi);

might be written as::

  Field3D phi;
  ...
  Field3D result;

  // Create accessors for fast (unsafe) data access
  auto phi_acc = FieldAccessor<>(phi);
  auto result_acc = FieldAccessor<>(result);

  BOUT_FOR(i, result.region("RGN_NOBNDRY")) {
    result_acc[i] = DDX(phi_acc, i);
  }

For a simple example like this the code is harder to read, and there
is not much benefit because there is one loop over the domain in both
cases. The benefit becomes more apparent when multiple operations are
combined.

The operators are implemented in a header file, so that they can be
inlined. These are in ``include/bout/single_index_ops.hxx``. Table
:numref:`tab-gpusupport-singleindexfunctions` lists the functions
which have been implemented.

.. _tab-gpusupport-singleindexfunctions:
.. table:: Single index operator functions

   +------------------------------ +-------------------------------------------+
   | Function                      | Description                               |
   +===============================+===========================================+
   | DDX(F3D, i)                   | Derivative in X with ``ddx:first=C2``     |
   +-------------------------------+-------------------------------------------+
   | DDY(F3D, i)                   | Derivative in Y with ``ddy:first=C2``     |
   +-------------------------------+-------------------------------------------+
   | DDZ(F3D, i)                   | Derivative in Z with ``ddz:first=C2``     |
   +-------------------------------+-------------------------------------------+
   | bracket(F2D, F3D, i)          | bracket(F2D, F3D, BRACKET_ARAKAWA)        |
   +-------------------------------+-------------------------------------------+
   | bracket(F3D, F3D, i)          | bracket(F3D, F2D, BRACKET_ARAKAWA)        |
   +-------------------------------+-------------------------------------------+
   | Delp2(F3D, i)                 | Delp2 with useFFT=false, C2 central diff. |
   +-------------------------------+-------------------------------------------+
   | Div_par_Grad_par(F3D, i)      | 2nd order central difference              |
   +-------------------------------+-------------------------------------------+
   | b0xGrad_dot_Grad(F3D, F2D, i) | C2 central diff. for all derivatives      |
   +-------------------------------+-------------------------------------------+
   | b0xGrad_dot_Grad(F2D, F3D, i) | C2 central diff. for all derivatives      |
   +-------------------------------+-------------------------------------------+
   | D2DY2(F3D, i)                 | C2 2nd-order derivative ``ddy:second=C2`` |
   +-------------------------------+-------------------------------------------+
   | Grad_par(F3D, i)              | C2 derivative, ``ddy:first=C2``           |
   +-------------------------------+-------------------------------------------+

Note that for efficiency the method used in the single index operators
cannot be changed at runtime: The numerical method is fixed at compile
time. The ``DDX`` single index operator, for example, always uses 2nd
order central difference.

Unit tests which ensure that the single index operators produce the
same result as the original field operations are in
``tests/unit/include/bout/test_single_index_ops.cxx``. Note that there
are limitations to these unit tests, in particular the geometry may
not be fully exercised. Additional errors are possible when combining
these methods, or porting code from field operations to single index
operations. It is therefore essential to have integrated tests and
benchmarks for each model implementation: Unit tests are necessary
but not sufficient for correctness.

CoordinatesAccessor
-------------------

The differential operators used in physics models typically need
access to grid spacing (eg. dx), non-uniform grid corrections
(e.g. d1_dx), and multiple coordinate system fields (metric tensor
components). When a ``FieldAccessor`` is created from a field, it uses the
field's coordinate system to create a ``CoordinateAccessor``, which
provides fast access to this extra data.

The coordinate system data is usually stored in separate arrays, so
that many different pointers would need to be passed to the CUDA
kernels to use this data directly. This was found to cause compilation
errors with ``nvcc`` along the lines of "Formal parameter space
overflowed".

The ``CoordinatesAccessor`` reduces the number of parameters (data
pointers) by packing all ``Coordinates`` data (grid spacing, metric
tensor components) into a single array. The ordering of this data in
the array has not been optimised, but is currently striped: Data for
the same grid cell is close to each other in memory. Some guidance on
memory layout can be found `on the NVidia website
<https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html#coalesced-access-to-global-memory>`_ and might be used to improve performance in future. It is
likely that the results might be architecture dependent.

To minimise the number of times this data needs to be copied from
individual fields into the single array, and then copied from CPU to
GPU, ``CoordinatesAccessor``s are cached. A map (``coords_store``
defined in ``coordinates_accessor.cxx``) associates
``Array<BoutReal>`` objects (containing the array of data) to
``Coordinates`` pointers. If a ``CoordinatesAccessor`` is constructed
with a ``Coordinates`` pointer which is in the cache, then the
previously created ``Array`` data is used.
Some care is therefore needed if the ``Coordinates`` data is modified,
to ensure that a new ``CoordinatesAccessor`` data array is created by
clearing the old data from the cache.

The easiest way to clear the cache is to call the static function
``CoordinatesAccessor::clear()``, which will delete all arrays from
the cache. To remove a single ``Coordinates`` key from the cache, pass
the pointer to ``CoordinatesAccessor::clear(coordinates_ptr)``.  In
both cases the number of keys removed from the cache will be returned.

Examples
--------

The ``blob2d-outerloop`` example is the simplest one which uses single index operators
and (optionally) RAJA. It should solve the same set of equations, with the same inputs,
as `blob2d`.

``hasegawa-wakatani-3d`` is a 3D turbulence model, typically solved in a slab geometry.

``elm-pb-outerloop`` is a much more complicated model, which should solve the same
equations, and have the same inputs, as ``elm-pb``. Note that there are some differences:

* The numerical methods used in ``elm-pb`` can be selected at
  run-time, and typically include WENO schemes e.g W3. In
  ``elm-pb-outerloop`` the methods are fixed to C2 in all cases.
* The equations solved by ``elm-pb`` can be changed by modifying input settings.
  To achieve higher performance, ``elm-pb-outerloop`` does this at compile time.
  There are checks to ensure that the code has been compiled with flags consistent
  with the input settings. See the README file for more details.


Notes:

* When RAJA is used in a physics model, all members of the PhysicsModel
  should be public. If this is not done, then a compiler error like
  "error: The enclosing parent function ("rhs") for an extended __device__ lambda
  cannot have private or protected access within its class" may be encountered.

* Class member variables cannot in general be used inside a RAJA loop: The ``this``
  pointer is captured by value in the lambda function, not the value of the member variable. 
  When the member variable is used on the GPU, the ``this`` pointer is generally not valid
  (e.g. on NERSC Cori GPUs). Some architectures have Address Translation Services (ATS)
  which enable host pointers to be resolved on the GPU.


