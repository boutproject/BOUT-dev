.. _sec-newv5:

=============================
 New Features in BOUT++ v5.0
=============================

BOUT++ v5.0 is a new major release, adding tons of new features, improving
existing code, as well as removing some old and deprecated things. There are
some breaking changes which will require modifications to your physics model,
but for the vast majority we have provided tooling to automate this as much as
possible.


3D Metrics
==========

Up until now, BOUT++ has been limited to varying the metric components only in
the XY plane. This release now introduces 3D metrics as a compile-time option,
allow simulations of devices such as stellarators.

To enable 3D metrics, build BOUT++ like:

.. code-block:: console

   ./configure --enable-metric-3D

or, with CMake:

.. code-block:: console

   cmake . -B build -DBOUT_ENABLE_METRIC_3D=ON

Changes
-------

Types
~~~~~

Adding 3D metrics to BOUT++ has been a substantial effort, requiring many
changes to a significant amount of the source code. The main change is that the
metric components, ``g11``, ``g22``, and so on, as well as the grid spacing,
``dx``, ``dy``, ``dz``, have changed from `Field2D` to
`Coordinates::FieldMetric`: a type alias for either `Field3D` or `Field2D`
depending on if BOUT++ was built with or without 3D metrics respectively.

.. note::
   `Coordinates::dz` has also changed to be `Field2D` even without 3D
   metrics. This is a breaking change, in that it may be necessary to change
   user code in order to keep working. If you don't use 3D metrics, wrapping the
   use of ``dz``, and similarly `Coordinates::zlength()`, in a call to
   `getUniform` will return a ``BoutReal``.


The use of `Coordinates::FieldMetric` has been followed through the rest of the
code base. If a metric component enters an expression that previously contained
only `Field2D` and `BoutReal` types, the result is now a
`Coordinates::FieldMetric`. This means that functions that previously both took
and returned a `Field2D` now return a `Coordinates::FieldMetric` (we could have
chosen to make the return type ``auto`` instead and rely on the compiler to
deduce the correct type, but we have chosen to make the dependence on the metric
dimensionality more explicit).

Because almost any operation on a vector involves the metric, the individual
components of `Vector2D` are now also of type
`Coordinates::FieldMetric`. Realistically, the use of `Vector2D` in a model
making use of 3D metrics is probably ill-advised.

Indexing
~~~~~~~~

3D metrics also requires changes in how fields are indexed. In `BOUT_FOR` loops,
generally no changes are required, as they already do The Right Thing. In other
cases, simply changing, for example, ``dx(x, y)`` to ``dx(x, y, z)`` is
sufficient: in the 2D metric case, the third index is accepted and discarded.

Many methods and operators have been upgraded to deal with 3D metrics. For
example, the `LaplaceXZpetsc` implementation has been modified to deal with
non-zero ``g_{xz}`` terms.

FIXME WHEN COORDINATES REFACTORED
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to simplify a lot of code, the call to `Coordinates::geometry`, which
calculates the connection coefficients `Coordinates::G1_11` and so on, has been
moved out of the `Coordinates` constructor. This is because computing the
coefficients involves derivatives which requires `Coordinates` and causes all
sorts of headaches and shims. As most users do not call the constructor
themselves anyway, this change should not be much of an issue.


Incompatibilities
-----------------

Many features of BOUT++ have been written assuming an axisymmetric coordinate
system. Once 3D metrics are enabled, this is no longer (necessarily) true which
breaks several features. For instance, many of the Laplacian inversion solvers
use intrinsically 2D methods, and so are not available when using 3D
metrics. Most of these features are runtime options, and therefore will throw an
exception if you try to use them. To get a list of available Laplacian solvers,
for example, you can pass the ``--list-laplacians`` flag to a compiled BOUT++
executable, which will print all the Laplacian solvers, noting which are
unavailable and why.

Several boundary conditions are also incompatible with 3D metrics, unfortunately
at the time of writing there is no easy way to list those that are. Several of
these, such as ``zerolaplace`` have no alternative implementations, so this may
mean it is not possible to run a given model with 3D metrics.

There are a few tests that don't work with 3D metrics, mostly because they rely
on one of the above incompatible methods or operators.

There is a preprocessor macro, ``BOUT_USE_METRIC_3D``, and a ``constexpr bool``,
`bout::build::use_metric_3d`, which can be used to guard code that doesn't
compile or work with 3D metrics, or perhaps needs to be handled differently.

Caution should be exercised with FFT-based methods. Technically, FFTs do work
with 3D metrics, but will not give the correct answer with non-constant ``dz``.
