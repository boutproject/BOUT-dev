.. _sec-derivatives-of-fft:

Derivatives of the Fourier transform
====================================

By using the definition of the Fourier transformed, we have

.. math::

   F(x,y,\xi) = {\int_{-\infty}^{\infty} {f(x,y,z)\exp(-2\pi iz\xi)} \; \text{d} {z}}

this gives

.. math::
   :label: f_derivative

   &{\int_{-\infty}^{\infty} {(\partial_zf[x,y,z])\exp(-2\pi iz\xi)} \; \text{d} {z}}\\
   =& {\int_{-\infty}^{\infty} {\partial_z(f[x,y,z]\exp[-2\pi iz\xi])} \; \text{d} {z}}
   - {\int_{-\infty}^{\infty} {f(x,y,z)\partial_z\exp(-2\pi iz\xi)} \; \text{d} {z}}\\
   =& (f[x,y,z]\exp[-2\pi iz\xi])\bigg|_{-\infty}^{\infty} - (-2\pi
   i\xi){\int_{-\infty}^{\infty} {f(x,y,z)\exp(-2\pi iz\xi)} \; \text{d} {z}}\\
   =& 2\pi i\xi F(x,y,\xi)

where we have used that :math:`f(x,y,\pm\infty)=0` in order to have a
well defined Fourier transform. This means that

.. math::

   \partial_z^n F(x,y,\xi) = (2\pi i \xi)^n F(x,y,\xi)

In our case, we are dealing with periodic boundary conditions. Strictly
speaking, the Fourier transform does not exist in such cases, but it is
possible to define a Fourier transform in the limit which in the end
lead to the Fourier series [1]_ By discretising the spatial domain, it
is no longer possible to represent the infinite amount of Fourier modes,
but only :math:`N+1` number of modes, where :math:`N` is the number of
points (this includes the modes with negative frequencies, and the
zeroth offset mode). For the discrete Fourier transform, we have

.. math::
   :label: DFT

   F(x,y)_{k} = \frac{1}{N}\sum_{Z=0}^{N-1}f(x,y)_{Z}\exp(\frac{-2\pi i k Z}{N})

where :math:`k` is the mode number, :math:`N` is the number of points
in :math:`z`. If we call the sampling points of :math:`z` for
:math:`z_Z`, where :math:`Z = 0, 1 \ldots N-1`, we have that
:math:`z_Z = Z \text{d}z`. As our domain goes from :math:`[0, 2\pi[`,
we have that (since we have one less line segment than point)
:math:`\text{d}z (N-1) = L_z = 2\pi - \text{d}z`, which gives
:math:`\text{d}z = \frac{2\pi}{N}`.  Inserting this is equation
(:eq:`DFT`) yields

.. math::

   F(x,y)_{k} = \frac{1}{N}\sum_{Z=0}^{N-1}f(x,y)_{Z}\exp( - i k
   Z\text{d}z) = \frac{1}{N}\sum_{Z=0}^{N-1}f(x,y)_{Z}\exp( - i k z_Z)

The discrete version of equation (:eq:`f_derivative`) thus gives

.. math::

   \partial_z^n F(x,y)_k = (i k)^n F(x,y)_k

.. [1] For more detail see Bracewell, R. N. - The Fourier Transform
       and Its Applications 3rd Edition chapter 10
