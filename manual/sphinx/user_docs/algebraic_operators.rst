.. _sec-algebraic-ops:

Algebraic operators
=========================

BOUT++ provides a wide variety of algebraic operators acting on fields.

The algebraic operators are listed in :numref:`tab-algebraic-ops`.
For a completely up-to-date list, see the ``Non-member functions``
part of :doc:`field2d.hxx<../_breathe_autogen/file/field2d_8hxx>`,
:doc:`field3d.hxx<../_breathe_autogen/file/field3d_8hxx>`,
:doc:`fieldperp.hxx<../_breathe_autogen/file/fieldperp_8hxx>`.

.. _tab-algebraic-ops:
.. table:: Algebraic operators

   +------------------------------------------+------------------------------------------------------+ 
   |  Name                                    | Description                                          |
   +==========================================+======================================================+
   | ``min(f, allpe=true, region)``           | Minimum (optionally over all processes)              | 
   +------------------------------------------+------------------------------------------------------+
   | ``max(f, allpe=true, region)``           | Maximum (optionally over all processes)              |
   +------------------------------------------+------------------------------------------------------+
   | ``pow(lhs, rhs, region)``                | :math:`\mathtt{lhs}^\mathtt{rhs}`                    |
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

These operators take a ``region`` argument, whose values can be [#]_ (see
:ref:`sec-iterating`)

-  `RGN_ALL`, which is the whole mesh;

-  `RGN_NOBNDRY`, which skips all boundaries;

-  `RGN_NOX`, which skips the x boundaries

-  `RGN_NOY`, which skips the y boundaries

The default value for the region argument is `RGN_ALL` which should work in all
cases.  However, the region argument can be used for optimization, to skip
calculations in guard cells if it is known that those results will not be
needed (for example, if no derivatives of the result will be calculated). Since
these operators can be relatively expensive compared to addition, subtraction,
multiplication this can be a useful performance improvement.

.. [#] More regions may be added in future, for example to act on only subsets of the
       physical domain.
