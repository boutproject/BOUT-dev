.. _sec-algebraic-ops:

Algebraic operators
===================

BOUT++ provides a wide variety of algebraic operators acting on fields.

Most of these operators can participate in the lazy field-expression
system described in :ref:`sec-field-expressions`. In practice this means
you can usually write ordinary algebraic code and let BOUT++ delay
evaluation until assignment or reduction.

For a completely up-to-date list, see the ``Non-member functions`` part
of :doc:`field2d.hxx<../_breathe_autogen/file/field2d_8hxx>`,
:doc:`field3d.hxx<../_breathe_autogen/file/field3d_8hxx>`, and
:doc:`fieldperp.hxx<../_breathe_autogen/file/fieldperp_8hxx>`.

Common operators
----------------

.. _tab-algebraic-ops:
.. table:: Algebraic operators

   +------------------------------------------+------------------------------------------------------+
   | Name                                     | Description                                          |
   +==========================================+======================================================+
   | ``min(f, allpe=true, region)``           | Minimum (optionally over all processes)              |
   +------------------------------------------+------------------------------------------------------+
   | ``max(f, allpe=true, region)``           | Maximum (optionally over all processes)              |
   +------------------------------------------+------------------------------------------------------+
   | ``mean(f, allpe=true, region)``          | Mean (optionally over all processes)                 |
   +------------------------------------------+------------------------------------------------------+
   | ``pow(lhs, rhs, region)``                | :math:`\mathtt{lhs}^\mathtt{rhs}`                    |
   +------------------------------------------+------------------------------------------------------+
   | ``SQ(f, region)``                        | Square of ``f``                                      |
   +------------------------------------------+------------------------------------------------------+
   | ``sqrt(f, region)``                      | :math:`\sqrt{(f)}`                                   |
   +------------------------------------------+------------------------------------------------------+
   | ``abs(f, region)``                       | :math:`|f|`                                          |
   +------------------------------------------+------------------------------------------------------+
   | ``exp(f, region)``                       | :math:`e^f`                                          |
   +------------------------------------------+------------------------------------------------------+
   | ``log(f, region)``                       | :math:`\log(f)`                                      |
   +------------------------------------------+------------------------------------------------------+
   | ``sin(f, region)``                       | :math:`\sin(f)`                                      |
   +------------------------------------------+------------------------------------------------------+
   | ``cos(f, region)``                       | :math:`\cos(f)`                                      |
   +------------------------------------------+------------------------------------------------------+
   | ``tan(f, region)``                       | :math:`\tan(f)`                                      |
   +------------------------------------------+------------------------------------------------------+
   | ``sinh(f, region)``                      | :math:`\sinh(f)`                                     |
   +------------------------------------------+------------------------------------------------------+
   | ``cosh(f, region)``                      | :math:`\cosh(f)`                                     |
   +------------------------------------------+------------------------------------------------------+
   | ``tanh(f, region)``                      | :math:`\tanh(f)`                                     |
   +------------------------------------------+------------------------------------------------------+
   | ``floor(f, region)``                     | Returns a field with the floor of `f` at each point  |
   +------------------------------------------+------------------------------------------------------+
   | ``filter(f, n, region)``                 | Calculate the amplitude of the Fourier mode in the   |
   |                                          | z-direction with mode number `n`                     |
   +------------------------------------------+------------------------------------------------------+
   | ``lowpass(f, nmax, region)``             | Remove Fourier modes (in the z-direction) with mode  |
   |                                          | number higher than `zmax`                            |
   +------------------------------------------+------------------------------------------------------+
   | ``lowpass(f, nmax, nmin, region)``       | Remove Fourier modes (in the z-direction) with mode  |
   |                                          | number higher than `zmax` or lower than `zmin`       |
   +------------------------------------------+------------------------------------------------------+
   | ``shiftZ(f, angle, region)``             | Rotate `f` by `angle` in the z-direction.            |
   |                                          | :math:`\mathtt{angle}/2\pi` is the fraction of the   |
   |                                          | domain multiplied by :math:`2\pi` so angle is in     |
   |                                          | radians if the total size of the domain is           |
   |                                          | :math:`2\pi`                                         |
   +------------------------------------------+------------------------------------------------------+
   | ``DC(f, region)``                        | The average in the z-direction of `f`                |
   |                                          | (DC stands for direct current, i.e. the constant part|
   |                                          | of `f` as opposed to the AC, alternating current, or |
   |                                          | fluctuating part)                                    |
   +------------------------------------------+------------------------------------------------------+
   | ``if_else(cond, lhs, rhs)``              | Select between two algebraic branches                |
   +------------------------------------------+------------------------------------------------------+
   | ``if_else_zero(cond, expr)``             | Select either ``expr`` or zero                       |
   +------------------------------------------+------------------------------------------------------+

These operators can usually be combined directly in expressions::

  Field3D rhs = sqrt(SQ(n) + SQ(T));
  Field3D masked = if_else(use_drive, source * profile, sink * profile);
  BoutReal max_error = max(abs(lhs - rhs), true);

Reductions such as ``min``, ``max``, and ``mean`` can operate directly
on an expression, so an intermediate field is often unnecessary.

Region arguments
----------------

These operators take a ``region`` argument. Common values are [#]_ (see
:ref:`sec-iterating`):

- ``RGN_ALL``, which is the whole mesh
- ``RGN_NOBNDRY``, which skips all boundaries
- ``RGN_NOX``, which skips the x boundaries
- ``RGN_NOY``, which skips the y boundaries

The default is usually ``RGN_ALL``. Restricting the region can improve
performance when guard-cell values will not be used.

When a region-limited expression is materialized into a field, only the
selected region is guaranteed to contain valid values. This is the same
performance-oriented convention used by other field operators.

Further reading
---------------

- :ref:`sec-field-expressions`
- :ref:`sec-gpusupport`

.. [#] More regions may be added in future, for example to act on only
       subsets of the physical domain.
