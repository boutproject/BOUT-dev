.. default-role:: math

=======================================
Geometry and Differential Operator
=======================================

:Author: X. Q. Xu

Geometry
========

In a axisymmetric toroidal system, the magnetic field can be expressed
as

.. math:: {\bf B}=I(\psi)\nabla\zeta+\nabla\zeta\times\nabla\psi,

where `\psi` is the poloidal flux, `\theta` is the
poloidal angle-like coordinate, and `\zeta` is the toroidal
angle. Here, `I(\psi)=RB_t`. The two important geometrical
parameters are: the curvature, `\bf \kappa`, and the local
pitch, `\nu(\psi,\theta)`,

.. math:: \nu(\psi,\theta)= {I(\psi){\bf \cal J}/R^2}.

The local pitch `\nu(\psi,\theta)` is related to the MHD safety
q by `\hat q(\psi)={2\pi}^{-1}\oint\nu(\psi,\theta) d\theta` in
the closed flux surface region, and `\hat
q(\psi)={2\pi}^{-1}\int_{inboard}^{outboard}\nu(\psi,\theta) d\theta`
in the scrape-off-layer. Here `{\bf \cal
J}=(\nabla\psi\times\nabla\theta\cdot\nabla\zeta)^{-1}` is the
coordinate Jacobian, `R` is the major radius, and `Z` is
the vertical position.

Geometry and Differential Operators
===================================

In a axisymmetric toroidal system, the magnetic field can be expressed
as `{\bf B}=I(\psi)\nabla\zeta+\nabla\zeta\times\nabla\psi`, where
`\psi` is the poloidal flux, `\theta` is the poloidal
angle-like coordinate, and `\zeta` is the toroidal angle. Here,
`I(\psi)=RB_t`. The two important geometrical parameters are: the
curvature, `\bf \kappa`, and the local pitch,
`\nu(\psi,\theta)`, and
`\nu(\psi,\theta)= {I(\psi){\bf \cal J}/R^2}`. The local pitch
`\nu(\psi,\theta)` is related to the MHD safety q by
`\hat q(\psi)={2\pi}^{-1}\oint\nu(\psi,\theta) d\theta` in the
closed flux surface region, and
`\hat q(\psi)={2\pi}^{-1}\int_{inboard}^{outboard}\nu(\psi,\theta) d\theta`
in the scrape-off-layer. Here
`{\bf \cal J}=(\nabla\psi\times\nabla\theta\cdot\nabla\zeta)^{-1}`
is the coordinate Jacobian, `R` is the major radius, and `Z`
is the vertical position.

Differential Operators
----------------------

For such an axisymmetric equilibrium the metric coefficients are only
functions of `\psi` and `\theta`. Three spatial differential
operators appear in the equations given as:
`{\bf v_E}\cdot\nabla_\perp`, `\nabla_\|` and
`\nabla_\perp^2`.

.. math::

   \begin{aligned}
   \nabla_\|&=&{\bf b_0}\cdot\nabla={1\over {\cal J}B}{\partial\over\partial\theta}+{I\over BR^2}{\partial\over\partial\zeta}={B_p\over hB}{\partial\over\partial\theta}+{B_t\over RB}{\partial\over\partial\zeta}, \\
   {\cal J}\nabla^2&=&
   {\partial\over\partial\psi}\left({\cal J}J_{11}{\partial\over\partial\psi}\right)
   +{\partial\over\partial\psi}\left({\cal J}J_{12}{\partial\over\partial\theta}\right) \nonumber\\
   &+&{\partial\over\partial\theta}\left({\cal J}J_{22}{\partial\over\partial\theta}\right)
   +{\partial\over\partial\theta}\left({\cal J}J_{12}{\partial\over\partial\psi}\right)  \nonumber\\
   &+&{1\over R^2}{\partial^2\over\partial\zeta^2}. \\
   \nabla_\|^2&=&{\bf b}_0\cdot\nabla({\bf b}_0\cdot\nabla)={1\over {\cal J}B}{\partial\over\partial\theta}\left({1\over {\cal J}B}{\partial\over\partial\theta}\right)
   +{1\over {\cal J}B}{\partial\over\partial\theta}\left({B_t\over RB}{\partial\over\partial\zeta}\right) \\
   &+&{B_t\over {\cal J}RB^2}{\partial^2\over\partial\theta\partial\zeta}
   +\left({B_t\over {\cal J}RB}\right)^2{\partial^2\over\partial\zeta^2}, \\
   \nabla_\perp^2\Phi&=&-\nabla\cdot[{\bf b}\times({\bf b}\times\nabla\Phi)]=\nabla^2\Phi-(\nabla\cdot{\bf b})({\bf b}\cdot\nabla\Phi)-\nabla_\|^2\Phi\end{aligned}

where the coordinate Jacobian and metric coefficients are defined as
following:

.. math::

   \begin{aligned}
   {\cal J}&=&\nabla\psi\times\nabla\theta\cdot\nabla\zeta={h\over B_p}, \\
   h&=&\sqrt{Z_\theta^2+R_\theta^2}, \\
   J_{11}&=&|\nabla\psi|^2={R^2\over {\cal J}^2}(Z_\theta^2+R_\theta^2), \\
   J_{12}&=&J_{21}=\nabla\psi\cdot\nabla\theta=-{R^2\over {\cal J}^2}(Z_\theta Z_\psi+R_\psi R_\theta), \\
   J_{13}&=&J_{31}=0, \\
   J_{22}&=&|\nabla\theta|^2={R^2\over {\cal J}^2}(Z_\psi^2+R_\psi^2), \\
   J_{23}&=&J_{32}=0, \\
   J_{33}&=&|\nabla\zeta|^2={1\over R^2}.\end{aligned}

Concentric circular cross section inside the separatrix without the SOL
-----------------------------------------------------------------------

For concentric circular cross section inside the separatrix without the
SOL, the differential operators are reduced to:

.. math::

   R &= R_0+r\cos\theta, \\
   Z &= r\sin\theta, \\
   B_t &= {B_{t0}R_0\over R}, \\
   B_p &= {1\over R}{\partial\psi\over\partial r}, \\
   R_\psi &= {\cos\theta\over RB_p}, \\
   R_\theta &= -r\sin\theta, \\
   Z_\psi &= {\sin\theta\over RB_p}, \\
   Z_\theta &= r\cos\theta, \\
   {\cal J} &= {r\over B_p}, \\
   h &= r, \\
   J_{11} &= |\nabla\psi|^2=r^2B_p^2, \\
   J_{12} = J_{21} &= \nabla\psi\cdot\nabla\theta=0,\\
   J_{13} = J_{31} &= 0, \\
   J_{22} &= |\nabla\theta|^2={1\over r^2}, \\
   J_{23} = J_{32} &= 0, \\
   J_{33} &= |\nabla\zeta|^2={1\over R^2},\\
   \nabla^2 &\simeq {1\over r}{\partial\over\partial r}\left(r{\partial\over\partial r}\right)+{1\over r^2}{\partial^2\over\partial \theta^2}+{1\over R^2}{\partial^2\over\partial \zeta^2}


Field-aligned coordinates with `\theta` as the coordinate along the field line
----------------------------------------------------------------------------------------

A suitable coordinate mapping between field-aligned ballooning
coordinates (`x`, `y`, `z`) and the usual flux
coordinates (`\psi`, `\theta`, `\zeta`) is

.. math::

   \begin{aligned}
   x&=&\psi-\psi_s, \nonumber \\
   y&=&\theta, \nonumber \\
   z&=&\zeta-\int_{\theta_0}^\theta \nu(x,y)dy.\end{aligned}

as shown in Fig. 1. The covering area given by the square ABCD in the
usual flux coordinates is the same as the parallelogram ABEF in the
field-aligned coordinates. The magnetic separatrix is denoted by
`\psi=\psi_s`. In this choice of coordinates, `x` is a
flux surface label, `y`, the poloidal angle, is also the
coordinate along the field line, and `z` is a field line label
within the flux surface.

The coordinate Jacobian and metric coefficients are defined as
following:

.. math::

   \begin{aligned}
   {\cal J}&=&\nabla\psi\times\nabla\theta\cdot\nabla\zeta={h\over B_p}, \\
   h&=&\sqrt{Z_\theta^2+R_\theta^2}, \\
   {\cal J}_{11}&=&|\nabla x|^2={R^2\over {\cal J}^2}(Z_\theta^2+R_\theta^2), \\
   {\cal J}_{12}&=&{\cal J}_{21}=\nabla x\cdot\nabla y=-{R^2\over {\cal J}^2}(Z_\theta Z_\psi+R_\psi R_\theta), \\
   {\cal J}_{22}&=&|\nabla y|^2={R^2\over {\cal J}^2}(Z_\psi^2+R_\psi^2), \\
   {\cal J}_{13}&=&{\cal J}_{31}=\nabla x\cdot\nabla z=-\nu\nabla x\cdot\nabla y-|\nabla x|^2\left(\int_{y_0}^y {\partial \nu(x,y)\over\partial\psi}dy\right)=-|\nabla x|^2I_s, \\
   {\cal J}_{23}&=&{\cal J}_{32}=\nabla y\cdot\nabla z=-\nu|\nabla y|^2-\nu\nabla x\cdot\nabla y\left(\int_{y_0}^y {\partial \nu(x,y)\over\partial\psi}dy\right), \\
   {\cal J}_{33}&=&|\nabla z|^2=\left |\nabla\zeta-\nu\nabla \theta-\nabla\psi\left(\int_{y_0}^y {\partial \nu(x,y)\over\partial\psi}dy\right)\right |^2, \\
   I_s &=&  {{\cal J}_{12}\over|\nabla\psi|^2}\nu(x,y)+\left(\int_{y_0}^y {\partial \nu(x,y)\over\partial\psi}dy\right).\end{aligned}

Here `h` is the local minor radius, `I_s` is the
integrated local shear, and `y_0` is an arbitrary integration
parameter, which, depending on the choice of Jacobian, determines the
location where `I_s=0`. The disadvantage of this choice of
coordinates is that the Jacobian diverges near the X-point as
`B_p\rightarrow 0` and its effect spreads over the entire flux
surafces near the separatrix as the results of coordinate transform
`z`. Therefore a better set of coordinates is needed for X-point
divertor geometry. The derivatives are obtained from the chain rule as
follows:

.. math::

   \begin{aligned}
   {d\over d\psi}&=&{\partial\over \partial x} - \left(\int_{y_0}^y {\partial \nu(x,y)\over\partial\psi}dy\right){\partial\over \partial z},   \\
   {d\over d\theta}&=&{\partial\over \partial y} - \nu(x,y){\partial\over \partial z},   \\
   {d\over d\zeta}&=&{\partial\over \partial z}.\end{aligned}

In the field-aligned ballooning coordinates, the parallel differential
operator is simple, involving only one coordinate `y`

.. math::

   \begin{aligned}
   \partial_\|^0 &=&  {\bf b}_0\cdot\nabla_\|=\left({B_p\over hB}\right){\partial\over\partial y}.\end{aligned}

which requires a few grid points. The total axisymmetric drift
operator becomes

The perturbed `{\bf E}\times {\bf B}` drift operator becomes

.. math::

   \begin{aligned}
   {\delta\bf v_E}\cdot\nabla_\perp&=&
   {c\over BB_\|^*}\left\{
   -{I\over J}{\partial\langle\delta\phi\rangle\over\partial\theta}
   +{B_p^2}
   {\partial\langle\delta\phi\rangle\over\partial z}
   \right\}{\partial\over\partial\psi} \nonumber\\
   &+&{c\over BB_\|^*}\left\{{I\over{\cal J}}
   {\partial\langle\delta\phi\rangle\over\partial\psi}
   +{{\cal J}_{12}\over R^2}
   {\partial\langle\delta\phi\rangle\over\partial z}
   \right\}{\partial\over\partial\theta} \nonumber\\
   &-&{c\over BB_\|^*}\left\{B_p^2
   {\partial\langle\delta\phi\rangle\over\partial\psi}
   +{{\cal J}_{12}\over R^2}
   {\partial\langle\delta\phi\rangle\over\partial\theta}
   \right\}{\partial\over\partial z},\end{aligned}

when the conventional turbulence ordering (`k_\|\ll k_\perp`) is
used, the perturbed `{\bf E}\times {\bf B}` drift operator can
be further reduced to a simple form

.. math::

   \begin{aligned}
   {\delta\bf v_E}\cdot\nabla_\perp&=&
   {cB\over B_\|^*}\left(
   {\partial\langle\delta\phi\rangle\over\partial z}{\partial\over\partial x}
   -{\partial\langle\delta\phi\rangle\over\partial x}{\partial\over\partial z}\right)\end{aligned}

where `\partial/\partial\theta\simeq -\nu\partial/\partial z` is
used. In the perturbed `{\bf E}\times {\bf B}` drift operator
the poloidal and radial derivatives are written in the usual flux
`(\psi,\theta,\zeta)` coordinates in order to have various
options for valid discretizations. The general Laplacian operator for
potential is

.. math::

   \begin{aligned}
   {\cal J}\nabla^2\Phi&=&{\partial\over\partial x}\left({\cal J}{\cal J}_{11}{\partial\Phi\over\partial x}
   +{\cal J}{\cal J}_{12}{\partial\Phi\over\partial y}
   +{\cal J}{\cal J}_{13}{\partial\Phi\over\partial z}\right) \nonumber\\
   &+&{\partial\over\partial y}\left({\cal J}{\cal J}_{21}{\partial\Phi\over\partial x}
   +{\cal J}{\cal J}_{22}{\partial\Phi\over\partial y}
   +{\cal J}{\cal J}_{23}{\partial\Phi\over\partial z}\right) \nonumber\\
   &+&{\partial\over\partial z}\left({\cal J}{\cal J}_{31}{\partial\Phi\over\partial x}
   +{\cal J}{\cal J}_{32}{\partial\Phi\over\partial y}
   +{\cal J}{\cal J}_{33}{\partial\Phi\over\partial z}\right).\end{aligned}

 The general perpendicular Laplacian operator for potential is

.. math::

   \begin{aligned}
   {\cal J}\nabla_\perp^2\Phi&=&{\partial\over\partial x}\left({\cal J}{\cal J}_{11}{\partial\Phi\over\partial x}
   +{\cal J}{\cal J}_{12}{\partial\Phi\over\partial y}
   +{\cal J}{\cal J}_{13}{\partial\Phi\over\partial z}\right) \nonumber\\
   &+&{\partial\over\partial y}\left({\cal J}{\cal J}_{21}{\partial\Phi\over\partial x}
   +{\cal J}{\cal J}_{22}{\partial\Phi\over\partial y}
   +{\cal J}{\cal J}_{23}{\partial\Phi\over\partial z}\right) \nonumber\\
   &+&{\partial\over\partial z}\left({\cal J}{\cal J}_{31}{\partial\Phi\over\partial x}
   +{\cal J}{\cal J}_{32}{\partial\Phi\over\partial y}
   +{\cal J}{\cal J}_{33}{\partial\Phi\over\partial z}\right) \nonumber\\
   &-&\left({B_p\over hB}\right){\partial\over\partial y}
   \left[\left({B_p\over hB}\right){\partial\Phi\over\partial y}\right] \nonumber\\
   &-&\left({B_p\over hB}\right)^2{\partial\ln B\over\partial y}{\partial\Phi\over\partial y}.\end{aligned}

The general perpendicular Laplacian operator for axisymmetric
potential `\Phi_0(x,y)` is

.. math::

   \begin{aligned}
   {\cal J}\nabla_\perp^2\Phi_0&=&{\partial\over\partial x}\left({\cal J}{\cal J}_{11}{\partial\Phi_0\over\partial x}
   +{\cal J}{\cal J}_{12}{\partial\Phi_0\over\partial y}\right) \nonumber\\
   &+&{\partial\over\partial y}\left({\cal J}{\cal J}_{21}{\partial\Phi_0\over\partial x}
   +{\cal J}{\cal J}_{22}{\partial\Phi_0\over\partial y}\right) \nonumber\\
   &-&\left({B_p\over hB}\right){\partial\over\partial y}
   \left[\left({B_p\over hB}\right){\partial\Phi_0\over\partial y}\right]  \nonumber\\
   &-&\left({B_p\over hB}\right)^2{\partial\ln B\over\partial y}{\partial\Phi\over\partial y}.\end{aligned}

For the perturbed potential `\delta\phi`, we can drop the
`\partial/\partial y` terms in Eq. (69) due to the elongated
nature of the turbulence (`k_\|/k_\perp\ll1`). The general
perpendicular Laplacian operator for perturbed potential
`\delta\phi` reduces to

.. math::

   \begin{aligned}
   {\cal J}\nabla_\perp^2\delta\phi&=&
   {\partial\over\partial x}\left({\cal J}{\cal J}_{11}{\partial\delta\phi\over\partial x}
   +{\cal J}{\cal J}_{13}{\partial\delta\phi\over\partial z}\right) \nonumber\\
   &+&{\partial\over\partial z}\left({\cal J}{\cal J}_{31}{\partial\delta\phi\over\partial x}
   +{\cal J}{\cal J}_{33}{\partial\delta\phi\over\partial z}\right).\end{aligned}

If the non-split potential `\Phi` is a preferred option, the
gyrokinetic Poisson equation (18) and the general perpendicular
Laplacian operator Eq. (69) have to be used. Then the assumption
`k_\|/k_\perp\ll1` is not used to simplify the perpendicular
Laplacian operator.
