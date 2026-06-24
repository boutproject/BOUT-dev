.. _sec-field-expressions:

Field Expressions
=================

BOUT++ field algebra now supports *lazy expressions* for many common
operations. Instead of creating a temporary field for every ``+``,
``-``, ``*``, ``/``, ``sqrt`` or ``abs``, BOUT++ can keep the expression
symbolic and evaluate it only when a concrete field or scalar result is
needed.

This keeps ordinary model code readable while reducing temporary
allocations and extra loops over the mesh. It is especially helpful for
accelerator backends, where launching fewer kernels matters.

What stays lazy
---------------

The following operations can form lazy expressions over `Field2D`,
`Field3D`, `Field3DParallel`, and `FieldPerp` where the combination makes
sense:

- Arithmetic operators: ``+``, ``-``, ``*``, ``/``
- Unary algebraic operators such as ``sqrt``, ``abs``, ``exp``, ``log``,
  ``sin``, ``cos``, ``tan``, ``sinh``, ``cosh``, ``tanh``, ``floor``,
  and ``SQ``
- Simple conditionals with ``if_else`` and ``if_else_zero``
- Reductions such as ``min``, ``max``, and ``mean``

For example::

  Field3D n, T;
  Field3D result;

  result = sqrt(SQ(n) + SQ(T));

The right-hand side can stay lazy until the assignment to ``result``.

When evaluation happens
-----------------------

An expression is evaluated when BOUT++ needs actual storage or a scalar
answer. Common triggers are:

- assigning to a field
- constructing a field from an expression
- assigning a field expression into an `Options` object
- calling scalar reductions such as ``min``, ``max``, or ``mean``

Examples::

  Field3D result = n + T;
  options["rhs"] = n + T;
  BoutReal max_value = max(abs(n + T), true);

Region-limited expressions
--------------------------

Many algebraic operators take a ``region`` argument, usually defaulting
to ``RGN_ALL``. A lazy expression keeps track of that region.

Only values inside the requested region are guaranteed to be valid after
materialization. This is useful for skipping guard-cell work when the
result will only be used in a smaller region::

  Field3D interior = abs(n, "RGN_NOBNDRY");

As with other region-limited field operations in BOUT++, code that later
uses guard cells should communicate or otherwise fill those cells before
relying on them.

Metadata propagation
--------------------

When an expression is materialized into a field, BOUT++ propagates the
field metadata carried by the expression:

- mesh pointer
- cell location
- field directions
- for `FieldPerp`, the y-index

This means expressions are intended to behave like ordinary field
operations in user code. Compatibility checks still apply: combining
fields on different meshes or incompatible staggered locations is an
error.

Mixed field types
-----------------

Several mixed-type combinations are supported directly:

- `Field2D` with `Field3D`: the 2D quantity is broadcast in ``z``
- `FieldPerp` with matching perpendicular data: the operation uses the
  `FieldPerp` y-index
- expressions involving metric components may return
  `Coordinates::FieldMetric`, which is `Field2D` or `Field3D` depending
  on how BOUT++ was built

In practice, this means code such as::

  Coordinates::FieldMetric grad = coords->J / coords->g_22;
  Field3D rhs = density * temperature + background_2d;

can use the same algebraic style even when metric dimensionality or
field rank differs.

Conditionals
------------

``if_else`` selects between two algebraic branches without forcing the
branches to be precomputed::

  Field3D rhs = if_else(use_source, source * density, sink * density);

``if_else_zero(condition, expr)`` is a shorthand for selecting either an
expression or zero::

  Field3D rhs = if_else_zero(include_drive, drive * profile);

This is particularly convenient when optional source terms are enabled
or disabled by compile-time or run-time logic.

Reductions on expressions
-------------------------

Reductions can operate directly on expressions instead of requiring an
intermediate field::

  BoutReal rms = sqrt(mean(SQ(n - n0), true, "RGN_NOBNDRY"));
  BoutReal max_error = max(abs(lhs - rhs), true);

This is often clearer than explicitly constructing a temporary field,
and it avoids extra storage.

Relation to GPU execution
-------------------------

Lazy field expressions are the high-level path to reducing temporary
work. They are a good default when ordinary field algebra expresses the
operation clearly.

For more control, especially when you want to fuse derivative operators
into a single explicit loop, see :ref:`sec-gpusupport`.

See also
--------

- :doc:`algebraic_operators`
- :doc:`gpu_support`
- :doc:`differential_operators`
