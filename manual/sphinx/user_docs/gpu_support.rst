.. _sec-gpusupport:

GPU support
===========

This section describes work in progress to develop GPU support in BOUT++ models.

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
