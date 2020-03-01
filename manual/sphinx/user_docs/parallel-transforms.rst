.. _sec-parallel-transforms:

Parallel Transforms
===================

In most BOUT++ simulations the Y coordinate is parallel to the
magnetic field. In particular if the magnetic field :math:`\mathbf{B}`
can be expressed as

.. math::
   \mathbf{B} = \nabla z \times \nabla x

then the Clebsch operators can be used. See section :ref:`sec-diffops`
for more details.

The structure of the magnetic field can be simple, as in a slab
geometry, but in many cases it is quite complicated. In a tokamak, for
example, the magnetic shear causes deformation of grid cells and
numerical issues. One way to overcome this is to transform between
local coordinate systems, interpolating in the toroidal (Z) direction
when calculating gradients along the magnetic field. This is called
the *shifted metric* method. In more general geometries such as
stellarators, the magnetic field can have a 3D structure and
stochastic regions. In this case the interpolation becomes 2D (in X
and Z), and is known as the Flux Coordinate Independent (FCI) method.

To handle these different cases in the same code, the BOUT++ mesh
implements different `ParallelTransform` classes. Each `Field3D` class
contains a pointer to the values up and down in the Y direction,
called yup and ydown.  These values are calculated during
communication::

   Field3D f(0.0);  // f allocated, set to zero
   f.yup();    // error: f.yup not allocated

   mesh->communicate(f);
   f.yup(); // ok

   f.ydown()(0,1,0); // ok


In the case of slab geometry, yup and ydown point to the original
field (f).  For this reason the value of f along the magnetic field
from f(x,y,z) is given by f.ydown(x,y-1,z) and f.yup(x,y+1,z). To take
a second derivative along Y using the Field3D iterators (section
:ref:`sec-iterating`)::

   Field3D result;
   result.allocate(); // Need to allocate before indexing

   for(const auto &i : result.region(RGN_NOBNDRY)) {
     result[i] = f.yup()[i.yp()] - f.ydown()[i.ym()];
   }

Note the use of yp() and ym() to increase and decrease the Y index.

Field-aligned grid
------------------

The default `ParallelTransform` is the identity transform, which sets
yup() and ydown() to point to the same field. In the input options the
setting is

.. code-block:: bash

   [mesh]
   paralleltransform = identity


This then uses the `ParallelTransformIdentity` class to calculate the
yup and ydown fields.

This is mostly useful for slab geometries, where for a straight magnetic field
the grid is either periodic in the y-direction or ends on a y-boundary. By
setting the global option ``TwistShift = true`` and providing a ``ShiftAngle``
in the gridfile or ``[mesh]`` options a branch cut can be introduced between
the beginning and end of the y-domain.

`ParallelTransformIdentity` can also be used in non-slab geometries. Then
``TwistShift = true`` should be set so that a twist-shift boundary condition is
applied on closed field lines, as field-line following coordinates are not
periodic in poloidal angle. Note that it is not recommended to use
`ParallelTransformIdentity` with toroidal geometries, as magnetic shear will
make the radial derivatives inaccurate away from the outboard midplane (which
is normall chosen as the zero point for the integrated shear).

Shifted metric
--------------

The shifted metric method is selected using:

.. code-block:: bash

   [mesh]
   paralleltransform = shifted

so that mesh uses the `ShiftedMetric` class to calculate parallel
transforms.  During initialisation, this class reads a quantity zShift
from the input or grid file. If zShift is not found then qinty is read
instead. If qinty is not found then the angle is zero, and this method
becomes the same as the identity transform.  For each X and Z index,
the zShift variable should contain the toroidal angle of a magnetic
field line at :math:`z=0` starting at :math:`\phi=0` at a reference
location :math:`\theta_0`:

.. math::

   \mathtt{zShift} = \int_{\theta_0}^\theta \frac{B_\phi h_\theta}{B_\theta R} d\theta

Note that here :math:`\theta_0` does not need to be constant in X
(radius), since it is only the relative shifts between Y locations
which matters.

FCI method
----------

To use the FCI method for parallel transforms, set

.. code-block:: bash

   [mesh]
   paralleltransform = fci

which causes the `FCITransform` class to be used for parallel
transforms.  This reads four variables (3D fields) from the input
grid: `forward_xt_prime`, `forward_zt_prime`, `backward_xt_prime`, and
`backward_zt_prime`. These give the cell indices, not in general
integers, in the forward (yup) and backward (ydown) directions. These
are arranged so that forward_xt_prime(x,y,z) is the x index at
y+1. Hence f.yup()(x,y+1,z) is calculated using
forward_xt_prime(x,y,z) and forward_zt_prime(x,y,z), whilst
f.ydown()(x,y-1,z) is calculated using backward_xt_prime(x,y,z) and
backward_zt_prime(x,y,z).

Tools for calculating these mappings include Zoidberg, a Python tool
which carries out field-line tracing and generates FCI inputs.
