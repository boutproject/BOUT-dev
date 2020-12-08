.. default-role:: math

.. _sec-field-aligned-coordinates:

=========================
Field-aligned coordinates
=========================

:Author: B.Dudson§, M.V.Umansky, L.C.Wang, X.Q.Xu, L.L.LoDestro
         §Department of Physics, University of York, UK
         Lawrence Livermore National Laboratory, USA
         IFTS, China

Introduction
============

This manual covers the field-aligned coordinate system used in many
BOUT++ tokamak models, and useful derivations and expressions.

Orthogonal toroidal coordinates
===============================

Starting with an orthogonal, right-handed toroidal coordinate system
`\left(r, \theta, \phi\right)`. `\theta` is the poloidal angle (from
`0` to `2\pi`) in the clockwise direction in the right R-Z
plane. `\phi` is the toroidal angle (also `0` to `2\pi`) going
anti-clockwise from the top of the tokamak.

We define the poloidal magnetic field `B_{pol}` as the component of
the magnetic field in the `\theta` direction, and the toroidal field `B_{tor}`
as the component of the magnetic field in the `\phi` direction.

We now introduce the poloidal flux `\psi` as the new radial
coordinate.  If the poloidal magnetic field `B_{pol}` is positive
then `\psi` increases with radius; if `B_{pol}` is negative then
`\psi` decreases with radius. To keep the coordinate system
right-handed, we define a new toroidal coordinate `\zeta` which is
defined as `\zeta = \sigma_{B\theta}\phi`, where the sign of the
poloidal magnetic field is `\sigma_{B\theta} \equiv {B_{\text{pol}}}/
\left|{B_{\text{pol}}}\right|`. If `B_{pol} > 0` then `\zeta` is
anti-clockwise looking down from above the tokamak, and if `B_{pol} <
0` then `\zeta` is clockwise. This coordinate system `\left(\psi,
\theta, \zeta\right)` is orthogonal and right-handed.

The magnitudes of the basis vectors are

.. math::

   \begin{aligned}
   \left|{\boldsymbol{e}}_\psi\right| = \frac{1}{R\left|{B_{\text{pol}}}\right|} \qquad \left|\boldsymbol{e}_\theta\right| = {h_\theta}
   \qquad \left|\boldsymbol{e}_\zeta\right| = R
   \end{aligned}

where `{h_\theta}` is the poloidal arc length per radian.
The non-zero covariant metric coefficients are

.. math::

   \begin{aligned}
   g_{\psi\psi} = \frac{1}{\left(R\left|{B_{\text{pol}}}\right|\right)^2} \qquad g_{\theta\theta} =
   h_\theta^2 \qquad g_{\zeta\zeta} = R^2\end{aligned}

and the magnitudes of the reciprocal vectors are therefore

.. math::

   \begin{aligned}
   \left|\nabla\psi\right| = R\left|{B_{\text{pol}}}\right| \qquad \left|\nabla\theta\right| = \frac{1}{h_\theta}
   \qquad \left|\nabla\zeta\right| = \frac{1}{R}\end{aligned}

The cross products are:

.. math::

   \boldsymbol{e}_\psi\times\boldsymbol{e}_\theta = J \nabla\zeta \qquad 
   \boldsymbol{e}_\psi\times\boldsymbol{e}_\zeta = -J \nabla\theta \qquad
   \boldsymbol{e}_\theta\times\boldsymbol{e}_\zeta = J \nabla\psi

where `J = h_\theta / \left|{B_{\text{pol}}}\right|` is the Jacobian, which is
always positive. Similarly,

.. math::

   \begin{aligned}
   \nabla\psi \times \nabla\theta = \frac{1}{J} \boldsymbol{e}_\zeta \qquad
   \nabla\psi \times \nabla\zeta = - \frac{1}{J} \boldsymbol{e}_\theta \qquad
   \nabla\theta \times \nabla\zeta = \frac{1}{J} \boldsymbol{e}_\psi
   \end{aligned}

The magnetic field `{\boldsymbol{B}}` can be expressed as

.. math::

   \begin{aligned}
    {\boldsymbol{B}}=& B_\theta \nabla \theta + B_\zeta \nabla \zeta \\
    =& B_\theta \boldsymbol{e}^\theta + B_\zeta \boldsymbol{e}^\zeta \\
    =& B_{\text{pol}} \frac{\boldsymbol{e}_\theta}{h_\theta} + B_{\text{tor}} \frac{\boldsymbol{e}_\zeta}{R} \\
    =& {B_{\text{pol}}}\hat{{\boldsymbol{e}}}_\theta + {B_{\text{tor}}}\hat{{\boldsymbol{e}}}_\zeta\end{aligned}


Field-aligned coordinates
=========================

In order to efficiently simulate (predominantly) field-aligned
structures, the standard coordinate system used by BOUT++ models is a
Clebsch system where grid-points are aligned to the magnetic field
along the `y` coordinate.

To align to the magnetic field we define a local field line pitch `\nu`:

.. math::

   \begin{aligned}
   \nu\left(\psi, \theta\right) = \frac{{\boldsymbol{B}}\cdot\nabla\phi}{{\boldsymbol{B}}\cdot\nabla\theta} =
   \frac{{B_{\text{tor}}}{h_\theta}}{{B_{\text{pol}}}R}
   \end{aligned}
The sign of the poloidal field `{B_{\text{pol}}}` and toroidal field 
`{B_{\text{tor}}}` can be either + or -.

The field-aligned coordinates `\left(x,y,z\right)` are defined by:

.. math::
   :label: eq:coordtransform

   \begin{aligned}
   x = {\sigma_{B\theta}}\left(\psi - \psi_0\right) \qquad y = \theta \qquad z = \sigma_{B\theta}
   \left(\phi - \int_{\theta_0}^{\theta}\nu\left(\psi,\theta\right)d\theta\right)
   \end{aligned}

The coordinate system is chosen so that `x` increases radially
outwards, from plasma to the wall. The `y` coordinate increases in the
same direction as `\theta` i.e. clockwise in the right-hand poloidal
plane. The `z` coordinate increases in the same direction as `\zeta`
i.e.  anti-clockwise looking from the top if `B_{pol}>0` and clockwise
if `B_{pol} < 0`.

This coordinate system is right-handed if `B_{pol}>0`, and left-handed if `B_{pol}<0`. 

The reciprocal basis vectors are therefore

.. math::

   \begin{aligned}
   \nabla x = {\sigma_{B\theta}}\nabla \psi \qquad
   \nabla y = \nabla \theta \qquad
   \nabla z = \nabla\zeta - \sigma_{B\theta}\left[\int_{\theta_0}^\theta{\frac{\partial \nu\left(\psi,\theta\right)}{\partial \psi}} d\theta\right] \nabla\psi
   - \sigma_{B\theta}\nu\left(\psi, \theta\right)\nabla\theta
   \end{aligned}
  
The term in square brackets is the integrated local shear:

.. math::

   \begin{aligned}
   I = \int_{y_0}^y\frac{\partial\nu\left(x, y\right)}{\partial\psi}dy\end{aligned}

  
The basis vectors are:

.. math::
   
   \begin{aligned}
   \boldsymbol{e}_x =& J\left(\nabla y \times \nabla z\right) = {\sigma_{B\theta}} {\boldsymbol{e}}_\psi + I{\boldsymbol{e}}_\zeta \\
   \boldsymbol{e}_y =& J\left(\nabla z \times \nabla x\right) = {\boldsymbol{e}}_\theta + \nu{\boldsymbol{e}}_\phi \\
   \boldsymbol{e}_z =& J\left(\nabla x \times \nabla y\right) = {\boldsymbol{e}}_\zeta
   \end{aligned}
 
where `{\boldsymbol{e}}_\phi =
{\sigma_{B\theta}}{\boldsymbol{e}}_\zeta` is always anticlockwise when
seen from above the tokamak looking down. The direction of
`{\boldsymbol{e}}_\zeta` depends on the sign of the poloidal field
`\sigma_{B\theta}`.

Magnetic field
--------------

Magnetic field is given in Clebsch form by:

.. math::

   \begin{aligned}
   {\boldsymbol{B}}= \nabla z\times \nabla x = \frac{1}{J}{\boldsymbol{e}}_y\end{aligned}

The contravariant components of this are then

.. math::

   \begin{aligned}
   B^y = \frac{{B_{\text{pol}}}}{{h_\theta}} \qquad B^x = B^z = 0\end{aligned}

i.e. `{\boldsymbol{B}}` can be written as

.. math::

   \begin{aligned}
   {\boldsymbol{B}}= \frac{{B_{\text{pol}}}}{{h_\theta}}{\boldsymbol{e}}_y\end{aligned}

and the covariant components calculated using `g_{ij}` as

.. math::

   \begin{aligned}
   B_x = {\sigma_{B\theta}}{B_{\text{tor}}}I R \qquad B_y = \frac{B^2 {h_\theta}}{{B_{\text{pol}}}} \qquad B_z = {\sigma_{B\theta}}{B_{\text{tor}}}R\end{aligned}

The unit vector in the direction of equilibrium `{\boldsymbol{B}}` is
therefore

.. math::

   \begin{aligned}
   {\boldsymbol{b}} = \frac{1}{JB}{\boldsymbol{e}}_y = \frac{1}{JB}\left[g_{xy}\nabla x + g_{yy}\nabla y
   + g_{yz}\nabla z\right]\end{aligned}

Jacobian and metric tensors
---------------------------

The Jacobian of this coordinate system is

.. math::

   \begin{aligned}
   J^{-1} \equiv \left(\nabla x\times\nabla y\right)\cdot\nabla z = {B_{\text{pol}}}/ {h_\theta}\end{aligned}

which can be either positive or negative, depending on the sign of
`{B_{\text{pol}}}`. The contravariant metric tensor is
given by:

.. math::

   \begin{aligned}
   g^{ij} \equiv {\boldsymbol{e}}^i \cdot{\boldsymbol{e}}^j \equiv \nabla u^i \cdot \nabla u^j = \left(%
   \begin{array}{ccc}
   \left(R{B_{\text{pol}}}\right)^2 & 0 & -I\left(R{B_{\text{pol}}}\right)^2 \\
   0 & 1 / {h_\theta}^2 & -{\sigma_{B\theta}}\nu / {h_\theta}^2 \\
   -I\left(R{B_{\text{pol}}}\right)^2 & -{\sigma_{B\theta}}\nu / {h_\theta}^2 & I^2\left(R{B_{\text{pol}}}\right)^2 + B^2 /
   \left(R{B_{\text{pol}}}\right)^2
   \end{array}
   %
    \right)\end{aligned}

and the covariant metric tensor:

.. math::

   \begin{aligned}
   g_{ij} \equiv {\boldsymbol{e}}_i \cdot{\boldsymbol{e}}_j = \left(%
   \begin{array}{ccc}
   I^2 R^2 + 1 / {\left({R{B_{\text{pol}}}}\right)^2}& {\sigma_{B\theta}}{B_{\text{tor}}}{h_\theta}I R / {B_{\text{pol}}}& I R^2 \\
   {\sigma_{B\theta}}{B_{\text{tor}}}{h_\theta}I R / {B_{\text{pol}}}& B^2{h_\theta}^2 / {B_{\text{pol}}}^2 & {\sigma_{B\theta}}{B_{\text{tor}}}{h_\theta}R / {B_{\text{pol}}}\\
   I R^2 & {\sigma_{B\theta}}{B_{\text{tor}}}{h_\theta}R / {B_{\text{pol}}}& R^2
   \end{array}
   %
    \right)\end{aligned}

or equivalently:

.. math::

   \begin{aligned}
   g_{ij} = \left(%
   \begin{array}{ccc}
   I^2 R^2 + 1 / {\left({R{B_{\text{pol}}}}\right)^2}& {\sigma_{B\theta}} I \nu R^2 & I R^2 \\
   {\sigma_{B\theta}} I \nu R^2 & J^2B^2 & {\sigma_{B\theta}} \nu R^2 \\
   I R^2 & {\sigma_{B\theta}}\nu R^2 & R^2
   \end{array}
   %
   \right)\end{aligned}

Right-handed field-aligned coordinates
======================================

If the poloidal magnetic field is negative, i.e. anti-clockwise in the right-hand R-Z plane, then the above
coordinate system is left-handed and the Jacobian `J` is negative.
To obtain a consistently right-handed coordinate system, we have to reverse the direction of the `y` coordinate
when `B_{pol} < 0`:

This `\left(x,y,z\right)` coordinate system is defined by:

.. math::
   :label: eq:coordtransform2

   \begin{aligned}
   x = {\sigma_{B\theta}}\left(\psi - \psi_0\right) \qquad y = {\sigma_{B\theta}}\theta \qquad z = \sigma_{B\theta}
   \left(\phi - \int_{\theta_0}^{\theta}\nu\left(\psi,\theta\right)d\theta\right)
   \end{aligned}

The radial coordinate `x` always points outwards. The `y` coordinate
increases in the direction of the poloidal magnetic field: clockwise
in the right-hand poloidal plane if `B_{pol} > 0`, and anti-clockwise
otherwise.  The `z` coordinate increases in the same direction as
`\zeta` i.e.  anti-clockwise looking from the top if `B_{pol}>0` and
clockwise if `B_{pol} < 0`.

This is still a Clebsch coordinate system:

.. math::

   \begin{aligned}
   {\boldsymbol{B}}= \nabla z\times \nabla x = \frac{1}{J}{\boldsymbol{e}}_y
   \end{aligned}

but the Jacobian is now always positive:

.. math::

   \begin{aligned}
   J = h_\theta / \left|B_{\text{pol}}\right|
   \end{aligned}


The reciprocal basis vectors are

.. math::
   
   \begin{aligned}
   \nabla x =& {\sigma_{B\theta}} \nabla \psi \\
   \nabla y =& {\sigma_{B\theta}} \nabla \theta \\
   \nabla z =& \nabla \zeta - {\sigma_{B\theta}} I \nabla \psi - {\sigma_{B\theta}}\nu\nabla\theta
   \end{aligned}

and basis vectors

.. math::
   
   \begin{aligned}
   \boldsymbol{e}_x =& J\left(\nabla y \times \nabla z\right) = {\sigma_{B\theta}} {\boldsymbol{e}}_\psi + I{\boldsymbol{e}}_\zeta \\
   \boldsymbol{e}_y =& J\left(\nabla z \times \nabla x\right) = {\sigma_{B\theta}} {\boldsymbol{e}}_\theta + \nu{\boldsymbol{e}}_\zeta \\
   \boldsymbol{e}_z =& J\left(\nabla x \times \nabla y\right) = {\boldsymbol{e}}_\zeta
   \end{aligned}


The contravariant metric tensor is:

.. math::

   \begin{aligned}
   g^{ij} \equiv {\boldsymbol{e}}^i \cdot{\boldsymbol{e}}^j \equiv \nabla u^i \cdot \nabla u^j = \left(%
   \begin{array}{ccc}
   \left(R{B_{\text{pol}}}\right)^2 & 0 & -I\left(R{B_{\text{pol}}}\right)^2 \\
   0 & 1 / {h_\theta}^2 & -\nu / {h_\theta}^2 \\
   -I\left(R{B_{\text{pol}}}\right)^2 & -\nu / {h_\theta}^2 & I^2\left(R{B_{\text{pol}}}\right)^2 + B^2 /
   \left(R{B_{\text{pol}}}\right)^2
   \end{array}
   %
   \right)\end{aligned}

and the covariant metric tensor:

.. math::

   \begin{aligned}
   g_{ij} = \left(%
   \begin{array}{ccc}
   I^2 R^2 + 1 / {\left({R{B_{\text{pol}}}}\right)^2}& I \nu R^2 & I R^2 \\
   I \nu R^2 & J^2B^2 & \nu R^2 \\
   I R^2 & \nu R^2 & R^2
   \end{array}
   %
   \right)\end{aligned}

The differences from the previous coordinate system are that `g_{xy}`,
`g_{yz}`, `g^{yz}` and `J` are multiplied by `{\sigma_{B\theta}}`. If
`B_{pol} < 0` so the poloidal magnetic field is anticlockwise in the
right-hand R-Z plane, then the `y` direction changes.
   
Differential operators in field-aligned coordinates
===================================================

The derivative of a scalar field `f` along the *unperturbed*
magnetic field `{\boldsymbol{b}}_0` is given by

.. math::

   \begin{aligned}
   \partial^0_{||}f \equiv {\boldsymbol{b}}_0 \cdot\nabla f =
   \frac{1}{JB}{\frac{\partial f}{\partial y}} = \frac{{B_{\text{pol}}}}{B{h_\theta}}{\frac{\partial f}{\partial y}}\end{aligned}

Note that J could be positive or negative. The parallel divergence is given by

.. math::

   \begin{aligned}
   \nabla^0_{||}f = B_0\partial^0_{||}\left(\frac{f}{B_0}\right)\end{aligned}

Using equation :eq:`eq:general_laplacian`,
the Laplacian operator is given by

.. math::

   \begin{aligned}
   \nabla^2 = &\frac{\partial^2}{\partial x^2}\left|\nabla x\right|^2 +
       \frac{\partial^2}{\partial y^2}\left|\nabla y\right|^2 +
       \frac{\partial^2}{\partial z^2}\left|\nabla z\right|^2 \nonumber \\
       &-2\frac{\partial^2}{\partial x\partial z}I\left(R{B_{\text{pol}}}\right)^2 -
       2\frac{\partial^2}{\partial y\partial z}\frac{\nu}{h_\theta^2}\\
       &+\frac{\partial}{\partial x}\nabla^2x + \frac{\partial}{\partial
   y}\nabla^2y + \frac{\partial}{\partial z}\nabla^2z \nonumber\end{aligned}

Using equation :eq:`eq:laplace_expand` for
`\nabla^2x = G^x` etc, the values are

.. math::

   \begin{aligned}
   \nabla^2x = \frac{{B_{\text{pol}}}}{h_\theta}\frac{\partial}{\partial x}\left(h_\theta
   R^2{B_{\text{pol}}}\right) \qquad \nabla^2y = \frac{{B_{\text{pol}}}}{h_\theta}\frac{\partial}{\partial
   y}\left(\frac{1}{{B_{\text{pol}}}h_\theta}\right)\end{aligned}

.. math::

   \begin{aligned}
   \nabla^2z = -\frac{{B_{\text{pol}}}}{h_\theta}\left[\frac{\partial}{\partial x}\left(IR^2{B_{\text{pol}}}
   h_\theta\right) + \frac{\partial}{\partial y}\left(\frac{\nu}{{B_{\text{pol}}}h_\theta}\right)\right]\end{aligned}

Neglecting some parallel derivative terms, the perpendicular Laplacian
can be written:

.. math::

   \begin{aligned}
   \nabla_\perp^2= {\left({R{B_{\text{pol}}}}\right)^2}\left[{\frac{\partial^2 }{\partial {x}^2}} - 2I\frac{\partial^2}{\partial z\partial x} +
   \left(I^2 + \frac{B^2}{\left({R{B_{\text{pol}}}}\right)^4}\right){\frac{\partial^2 }{\partial {z}^2}}\right] + \nabla^2 x {\frac{\partial }{\partial x}} +
   \nabla^2 z{\frac{\partial }{\partial z}}\end{aligned}

The second derivative along the equilibrium field

.. math::

   \begin{aligned}
   \partial^2_{||}\phi = \partial^0_{||}\left(\partial^0_{||}\phi\right) =
   \frac{1}{JB}{\frac{\partial }{\partial y}}\left(\frac{1}{JB}\right){\frac{\partial \phi}{\partial y}}
   + \frac{1}{g_{yy}}\frac{\partial^2\phi}{\partial y^2}\end{aligned}

A common expression (the Poisson bracket in reduced MHD) is (from
equation :eq:`eq:brackets`)):

.. math::

   \begin{aligned}
   {\boldsymbol{b}}_0\cdot\nabla\phi\times\nabla A =
   \frac{1}{J^2B}\left[\left(g_{yy}{\frac{\partial \phi}{\partial z}} -
   g_{yz}{\frac{\partial \phi}{\partial y}}\right){\frac{\partial A}{\partial x}} + \left(g_{yz}{\frac{\partial \phi}{\partial x}} -
   g_{xy}{\frac{\partial \phi}{\partial z}}\right){\frac{\partial A}{\partial y}} + \left(g_{xy}{\frac{\partial \phi}{\partial y}} -
   g_{yy}{\frac{\partial \phi}{\partial x}}\right){\frac{\partial A}{\partial z}}\right]\end{aligned}

The perpendicular nabla operator:

.. math::

   \begin{aligned}
   \nabla_\perp \equiv& \nabla - {\boldsymbol{b}}\left({\boldsymbol{b}}\cdot\nabla\right) \\ =& \nabla
       x\left({\frac{\partial }{\partial x}} - \frac{g_{xy}}{\left(JB\right)^2}{\frac{\partial }{\partial y}}\right) + \nabla
       z\left({\frac{\partial }{\partial z}} - \frac{g_{yz}}{\left(JB\right)^2}{\frac{\partial }{\partial y}}\right)\end{aligned}

.. _sec:jxb_fac:

J x B in field-aligned coordinates
----------------------------------

Components of the magnetic field in field-aligned coordinates:

.. math::

   \begin{aligned}
   B^y = \frac{{B_{\text{pol}}}}{{h_\theta}} \qquad B^x = B^z = 0\end{aligned}

and

.. math::

   \begin{aligned}
   B_x = {\sigma_{B\theta}}{B_{\text{tor}}}I R \qquad B_y = \frac{B^2{h_\theta}}{{B_{\text{pol}}}} \qquad B_z = {\sigma_{B\theta}}{B_{\text{tor}}}R\end{aligned}

Calculate current `{\boldsymbol{J}}= \frac{1}{\mu}{\nabla\times
{\boldsymbol{B}} }`

.. math::

   \begin{aligned}
   \left({\nabla\times {\boldsymbol{B}} }\right)^x = \frac{1}{J}\left({\frac{\partial B_z}{\partial y}} - {\frac{\partial B_y}{\partial z}}\right) = 0\end{aligned}

since `{B_{\text{tor}}}R` is a flux-surface quantity, and
`{\boldsymbol{B}}` is axisymmetric.

.. math::

   \begin{aligned}
   \left({\nabla\times {\boldsymbol{B}} }\right)^y =& -{\sigma_{B\theta}}\frac{{B_{\text{pol}}}}{{h_\theta}}{\frac{\partial }{\partial x}}\left({B_{\text{tor}}}R\right) \\
       \left({\nabla\times {\boldsymbol{B}} }\right)^z =&
       \frac{{B_{\text{pol}}}}{{h_\theta}}\left[{\frac{\partial }{\partial x}}\left(\frac{B^2{h_\theta}}{{B_{\text{pol}}}}\right) -
       {\sigma_{B\theta}}{\frac{\partial }{\partial y}}\left({B_{\text{tor}}}I R\right)\right]\end{aligned}

The second term can be simplified, again using
`{B_{\text{tor}}}R` constant on flux-surfaces:

.. math::

   \begin{aligned}
   {\frac{\partial }{\partial y}}\left({B_{\text{tor}}}I R\right) = {\sigma_{B\theta}}{B_{\text{tor}}}R{\frac{\partial \nu}{\partial x}} \qquad \nu =
   \frac{{h_\theta}{B_{\text{tor}}}}{R{B_{\text{pol}}}}\end{aligned}

From these, calculate covariant components:

.. math::
   :label: eq:curlb_y

   \begin{aligned}
   \left({\nabla\times {\boldsymbol{B}} }\right)_x =& -{B_{\text{tor}}}I R {\frac{\partial }{\partial x}}\left({B_{\text{tor}}}R\right) +
       \frac{IR^2{B_{\text{pol}}}}{{h_\theta}}\left[{\frac{\partial }{\partial x}}\left(\frac{B^2{h_\theta}}{{B_{\text{pol}}}}\right) - {B_{\text{tor}}}
       R{\frac{\partial \nu}{\partial x}}\right] \nonumber\\
   %
   \left({\nabla\times {\boldsymbol{B}} }\right)_y =& -{\sigma_{B\theta}}\frac{B^2{h_\theta}}{{B_{\text{pol}}}}{\frac{\partial }{\partial x}}\left({B_{\text{tor}}}R\right) +
       {\sigma_{B\theta}}{B_{\text{tor}}}R\left[{\frac{\partial }{\partial x}}\left(\frac{B^2{h_\theta}}{{B_{\text{pol}}}}\right) - {B_{\text{tor}}}R{\frac{\partial \nu}{\partial x}}\right]
       \\
   %
   \left({\nabla\times {\boldsymbol{B}} }\right)_z =& -{B_{\text{tor}}}R{\frac{\partial }{\partial x}}\left({B_{\text{tor}}}R\right) +
       \frac{R^2{B_{\text{pol}}}}{{h_\theta}}\left[{\frac{\partial }{\partial x}}\left(\frac{B^2{h_\theta}}{{B_{\text{pol}}}}\right) - {B_{\text{tor}}}
       R{\frac{\partial \nu}{\partial x}}\right] \nonumber\end{aligned}

Calculate `{\boldsymbol{J}}\times{\boldsymbol{B}}` using

.. math::

   \begin{aligned}
   {\boldsymbol{e}}^i = \frac{1}{J}\left({\boldsymbol{e}}_j \times {\boldsymbol{e}}_k\right) \qquad {\boldsymbol{e}}_i =
   J\left({\boldsymbol{e}}^j \times {\boldsymbol{e}}^k\right) \qquad i,j,k \texttt{ cyc } 1,2,3\end{aligned}

gives

.. math::

   \begin{aligned}
   \mu_0 \left({\boldsymbol{J}}\times{\boldsymbol{B}}\right)^x =& \frac{1}{J}\left[\left({\nabla\times {\boldsymbol{B}} }\right)_y B_z -
   \left({\nabla\times {\boldsymbol{B}} }\right)_z B_y \right]\\ =& -\frac{{B_{\text{pol}}}^3
   R^2}{{h_\theta}}\left[{\frac{\partial }{\partial x}}\left(\frac{B^2{h_\theta}}{{B_{\text{pol}}}}\right) - {B_{\text{tor}}}R{\frac{\partial \nu}{\partial x}}\right]\end{aligned}

Covariant components of `\nabla P`:

.. math::

   \begin{aligned}
   \left(\nabla P\right)_x = {\frac{\partial P}{\partial x}} \qquad \left(\nabla P\right)_y = \left(\nabla P\right)_z = 0\end{aligned}

and contravariant:

.. math::

   \begin{aligned}
   \left(\nabla P\right)^x = {\left({R{B_{\text{pol}}}}\right)^2}{\frac{\partial P}{\partial x}} \qquad \left(\nabla P\right)^y = 0 \qquad
   \left(\nabla P\right)^z = -I{\left({R{B_{\text{pol}}}}\right)^2}{\frac{\partial P}{\partial x}}\end{aligned}

Hence equating contravariant x components of
`{\boldsymbol{J}}\times{\boldsymbol{B}}= \nabla P`,

.. math::
   :label: eq:xbalance

   \begin{aligned}
   {\frac{\partial }{\partial x}}\left(\frac{B^2{h_\theta}}{{B_{\text{pol}}}}\right) - {B_{\text{tor}}}
   R{\frac{\partial }{\partial x}}\left(\frac{{B_{\text{tor}}}{h_\theta}}{R{B_{\text{pol}}}}\right) + \frac{\mu_0{h_\theta}}{{B_{\text{pol}}}}{\frac{\partial P}{\partial x}} =
   0
   \end{aligned}

Use this to calculate `{h_\theta}` profiles (need to fix
`{h_\theta}` at one radial location).

Close to x-points, the above expression becomes singular, so a better
way to write it is:

.. math::

   \begin{aligned}
   {\frac{\partial }{\partial x}}\left(B^2{h_\theta}\right) - {h_\theta}{B_{\text{pol}}}{\frac{\partial {B_{\text{pol}}}}{\partial x}} - {B_{\text{tor}}}
   R{\frac{\partial }{\partial x}}\left(\frac{{B_{\text{tor}}}{h_\theta}}{R}\right) + \mu_0{h_\theta}{\frac{\partial P}{\partial x}} = 0\end{aligned}

For solving force-balance by adjusting `P` and `f`
profiles, the form used is

.. math::

   \begin{aligned}
   {B_{\text{tor}}}{h_\theta}{\frac{\partial {B_{\text{tor}}}}{\partial x}} + \frac{{B_{\text{tor}}}^2{h_\theta}}{R}{\frac{\partial R}{\partial x}} +
   \mu_0{h_\theta}{\frac{\partial P}{\partial x}} = -{B_{\text{pol}}}{\frac{\partial }{\partial x}}\left({B_{\text{pol}}}{h_\theta}\right)\end{aligned}

A quick way to calculate f is to rearrange this to:

.. math::

   \begin{aligned}
   {\frac{\partial {B_{\text{tor}}}}{\partial x}} = {B_{\text{tor}}}\left[-\frac{1}{R}{\frac{\partial R}{\partial x}}\right] +
   \frac{1}{{B_{\text{tor}}}}\left[-\mu_0{\frac{\partial P}{\partial x}} -
   {\frac{\partial {B_{\text{pol}}}}{\partial {h_\theta}}}{\frac{\partial }{\partial x}}\left({B_{\text{pol}}}{h_\theta}\right)\right]\end{aligned}

and then integrate this using LSODE.

Parallel current
----------------

.. math::

   \begin{aligned}
   J_{||} = {\boldsymbol{b}}\cdot{\boldsymbol{J}}\qquad b^y = \frac{{B_{\text{pol}}}}{B{h_\theta}}\end{aligned}

and from equation :eq:`eq:curlb_y`:

.. math::

   \begin{aligned}
   J_y = \frac{{\sigma_{B\theta}}}{\mu_0}\left\{-\frac{B^2{h_\theta}}{{B_{\text{pol}}}}{\frac{\partial }{\partial x}}\left({B_{\text{tor}}}R\right) + {B_{\text{tor}}}
   R\left[{\frac{\partial }{\partial x}}\left(\frac{B^2{h_\theta}}{{B_{\text{pol}}}}\right) - {\sigma_{B\theta}}{B_{\text{tor}}}R{\frac{\partial \nu}{\partial x}}\right]\right\}\end{aligned}

since `J_{||} = b^yJ_y`,

.. math::

   \begin{aligned}
   \mu_0 J_{||} ={\sigma_{B\theta}}\frac{{B_{\text{pol}}}{B_{\text{tor}}}
   R}{B{h_\theta}}\left[{\frac{\partial }{\partial x}}\left(\frac{B^2{h_\theta}}{{B_{\text{pol}}}}\right) - {B_{\text{tor}}}R{\frac{\partial \nu}{\partial x}}\right] -
   {\sigma_{B\theta}}B{\frac{\partial }{\partial x}}\left({B_{\text{tor}}}R\right)\end{aligned}

Curvature
---------

For reduced MHD, need to calculate curvature term
`{\boldsymbol{b}}\times{\boldsymbol{\kappa}}`, where
`{\boldsymbol{\kappa}} =
\left({\boldsymbol{b}}\cdot\nabla\right){\boldsymbol{b}}=
-{\boldsymbol{b}}\times\left(\nabla\times{\boldsymbol{b}}\right)`.
Re-arranging, this becomes:

.. math::

   \begin{aligned}
   {\boldsymbol{b}}\times{\boldsymbol{\kappa}} = \nabla\times{\boldsymbol{b}}-
   {\boldsymbol{b}}\left({\boldsymbol{b}}\cdot\left(\nabla\times{\boldsymbol{b}}\right)\right)\end{aligned}

Components of `\nabla\times{\boldsymbol{b}}` are:

.. math::

   \begin{aligned}
   \left(\nabla\times{\boldsymbol{b}}\right)^x =& {\sigma_{B\theta}}\frac{{B_{\text{pol}}}}{{h_\theta}}{\frac{\partial }{\partial y}}\left(\frac{{B_{\text{tor}}}
   R}{B}\right) \\ \left(\nabla\times{\boldsymbol{b}}\right)^y =&
       -{\sigma_{B\theta}}\frac{{B_{\text{pol}}}}{{h_\theta}}{\frac{\partial }{\partial x}}\left(\frac{{B_{\text{tor}}}R}{B}\right) \\
       \left(\nabla\times{\boldsymbol{b}}\right)^z =&
       \frac{{B_{\text{pol}}}}{{h_\theta}}{\frac{\partial }{\partial x}}\left(\frac{B{h_\theta}}{{B_{\text{pol}}}}\right) - {\sigma_{B\theta}}\frac{{B_{\text{pol}}}{B_{\text{tor}}}
       R}{{h_\theta}B}{\frac{\partial \nu}{\partial x}} - {\sigma_{B\theta}}I\frac{{B_{\text{pol}}}}{{h_\theta}}{\frac{\partial }{\partial y}}\left(\frac{{B_{\text{tor}}}
       R}{B}\right) \\\end{aligned}

giving:

.. math::
   :label: eq:curvature

   \begin{aligned}
   {\boldsymbol{\kappa}} =& -\frac{{B_{\text{pol}}}}{B h_\theta}\left[{\frac{\partial }{\partial x}}\left(\frac{B
   h_\theta}{{B_{\text{pol}}}}\right) - {\sigma_{B\theta}}{\frac{\partial }{\partial y}}\left(\frac{{B_{\text{tor}}}I R}{B}\right)\right]\nabla x \nonumber
   \\ &+ {\sigma_{B\theta}}\frac{{B_{\text{pol}}}}{B h_\theta}{\frac{\partial }{\partial y}}\left(\frac{{B_{\text{tor}}}R}{B}\right)\nabla z
   \end{aligned}

.. math::

   \begin{aligned}
   {\boldsymbol{b}}\cdot\left(\nabla\times{\boldsymbol{b}}\right) = -{\sigma_{B\theta}}B{\frac{\partial }{\partial x}}\left(\frac{{B_{\text{tor}}}R}{B}\right) +
   {\sigma_{B\theta}}\frac{{B_{\text{tor}}}{B_{\text{pol}}}R}{B{h_\theta}}{\frac{\partial }{\partial x}}\left(\frac{B{h_\theta}}{{B_{\text{pol}}}}\right) -
   \frac{{B_{\text{pol}}}{B_{\text{tor}}}^2R^2}{{h_\theta}B^2}{\frac{\partial \nu}{\partial x}}\end{aligned}

therefore,

.. math::

   \begin{aligned}
   \left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\right)^x =& {\sigma_{B\theta}}\frac{{B_{\text{pol}}}}{{h_\theta}}{\frac{\partial }{\partial y}}\left(\frac{{B_{\text{tor}}}
   R}{B}\right) = -{\sigma_{B\theta}}\frac{{B_{\text{pol}}}{B_{\text{tor}}}R}{{h_\theta}B^2}{\frac{\partial B}{\partial y}} \\
   \left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\right)^y =& \frac{{B_{\text{pol}}}^2{B_{\text{tor}}}^2
   R^2}{B^3{h_\theta}^2}{\frac{\partial \nu}{\partial x}} - {\sigma_{B\theta}}\frac{{B_{\text{pol}}}^2{B_{\text{tor}}}
   R}{B^2{h_\theta}^2}{\frac{\partial }{\partial x}}\left(\frac{B{h_\theta}}{{B_{\text{pol}}}}\right) \\
   \left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\right)^z =&
   \frac{{B_{\text{pol}}}}{{h_\theta}}{\frac{\partial }{\partial x}}\left(\frac{B{h_\theta}}{{B_{\text{pol}}}}\right) - {\sigma_{B\theta}}\frac{{B_{\text{pol}}}{B_{\text{tor}}}
   R}{{h_\theta}B}{\frac{\partial \nu}{\partial x}} - I\left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\right)^x\end{aligned}

Using equation :eq:`eq:xbalance`:

.. math::

   \begin{aligned}
   B{\frac{\partial }{\partial x}}\left(\frac{B{h_\theta}}{{B_{\text{pol}}}}\right) + \frac{B{h_\theta}}{{B_{\text{pol}}}}{\frac{\partial B}{\partial x}} - {\sigma_{B\theta}}{B_{\text{tor}}}
   R{\frac{\partial }{\partial x}}\left(\frac{{B_{\text{tor}}}{h_\theta}}{R{B_{\text{pol}}}}\right) + \frac{\mu_0{h_\theta}}{{B_{\text{pol}}}}{\frac{\partial P}{\partial x}} =
   0\end{aligned}

we can re-write the above components as:

.. math::

   \begin{aligned}
   \left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\right)^y =& {\sigma_{B\theta}}\frac{{B_{\text{pol}}}{B_{\text{tor}}}
   R}{B^2{h_\theta}}\left[\frac{\mu_0}{B}{\frac{\partial P}{\partial x}} + {\frac{\partial B}{\partial x}}\right] \\
   \left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\right)^z =& -\frac{\mu_0}{B}{\frac{\partial P}{\partial x}} - {\frac{\partial B}{\partial x}} -
   I\left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\right)^x\end{aligned}

Curvature from `{\boldsymbol{b}}\times{\boldsymbol{\kappa}}`
------------------------

The vector `{\boldsymbol{b}}\times{\boldsymbol{\kappa}}` is an
approximation of

.. math::

   \begin{aligned}
   \frac{B}{2}\nabla\times\left(\frac{{\boldsymbol{b}}}{B}\right) \simeq {\boldsymbol{b}}\times{\boldsymbol{\kappa}}\end{aligned}

so can just derive from the original expression. Using the
contravariant components of `{\boldsymbol{b}}`, and the curl
operator in curvilinear coordinates (see appendix):

.. math::

   \begin{aligned}
   \nabla\times\left(\frac{{\boldsymbol{b}}}{B}\right) =&
       \frac{{B_{\text{pol}}}}{{h_\theta}}\left[\left({\frac{\partial }{\partial x}}\left(\frac{{h_\theta}}{{B_{\text{pol}}}}\right) -
       {\frac{\partial }{\partial y}}\left(\frac{{\sigma_{B\theta}}{B_{\text{tor}}}IR}{B^2}\right)\right){\boldsymbol{e}}_z \right.  \\ &+
       {\frac{\partial }{\partial y}}\left(\frac{{\sigma_{B\theta}}{B_{\text{tor}}}R}{B^2}\right){\boldsymbol{e}}_x \\ &+
       \left.{\frac{\partial }{\partial x}}\left(\frac{{\sigma_{B\theta}}{B_{\text{tor}}}R}{B^2}\right){\boldsymbol{e}}_y\right]\end{aligned}

This can be simplified using

.. math::

   \begin{aligned}
   {\frac{\partial }{\partial y}}\left(\frac{{\sigma_{B\theta}}{B_{\text{tor}}}IR}{B^2}\right) = I{\sigma_{B\theta}}{B_{\text{tor}}}
   R{\frac{\partial }{\partial y}}\left(\frac{1}{B^2}\right) + \frac{{B_{\text{tor}}}R}{B^2}{\frac{\partial \nu}{\partial x}}\end{aligned}

to give

.. math::

   \begin{aligned}
     \left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\right)^x =& -{\sigma_{B\theta}}\frac{{B_{\text{pol}}}{B_{\text{tor}}}R}{{h_\theta}B^2}{\frac{\partial B}{\partial y}} \\
       \left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\right)^y =& -{\sigma_{B\theta}}\frac{B{B_{\text{pol}}}}{2{h_\theta}}{\frac{\partial }{\partial x}}\left(\frac{{B_{\text{tor}}}
   R}{B^2}\right) \\ \left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\right)^z =&
       \frac{B{B_{\text{pol}}}}{2{h_\theta}}{\frac{\partial }{\partial x}}\left(\frac{{h_\theta}}{{B_{\text{pol}}}}\right) - \frac{{B_{\text{pol}}}{B_{\text{tor}}}
       R}{2{h_\theta}B}{\frac{\partial \nu}{\partial x}} - I\left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\cdot\nabla\right)^x\end{aligned}

The first and second terms in
`\left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\cdot\nabla\right)^z`
almost cancel, so by expanding out `\nu` a better expression is

.. math::

   \begin{aligned}
   \left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\right)^z = \frac{{B_{\text{pol}}}^3}{2{h_\theta}
   B}{\frac{\partial }{\partial x}}\left(\frac{{h_\theta}}{{B_{\text{pol}}}}\right) - \frac{{B_{\text{tor}}}
   R}{2B}{\frac{\partial }{\partial x}}\left(\frac{{h_\theta}}{{B_{\text{pol}}}}\right)\end{aligned}

Curvature of a single line
--------------------------

The curvature vector can be calculated from the field-line toroidal
coordinates `\left(R,Z,\phi\right)` as follows. The line element
is given by

.. math::

   \begin{aligned}
   d{\boldsymbol{r}} = dR{\hat{{\boldsymbol{R}}}}+ dZ{\hat{{\boldsymbol{Z}}}}+ Rd\phi{\hat{{\boldsymbol{\phi}}}}\end{aligned}

Hence the tangent vector is

.. math::

   \begin{aligned}
   \hat{{\boldsymbol{T}}} \equiv {\frac{d {\boldsymbol{r}}}{d s}} = {\frac{d R}{d s}}{\hat{{\boldsymbol{R}}}}+ {\frac{d Z}{d s}}{\hat{{\boldsymbol{Z}}}}+
   R{\frac{d \phi}{d s}}{\hat{{\boldsymbol{\phi}}}}\end{aligned}

where `s` is the distance along the field-line. From this, the
curvature vector is given by

.. math::

   \begin{aligned}
   {\boldsymbol{\kappa}}\equiv {\frac{d {\boldsymbol{T}}}{d s}} =& {\frac{d^2 R}{d s^2}}{\hat{{\boldsymbol{R}}}}+ {\frac{d R}{d s}}{\frac{d \phi}{d s}}{\hat{{\boldsymbol{\phi}}}}
       \\ &+ {\frac{d^2 Z}{d s^2}}{\hat{{\boldsymbol{Z}}}}\\ &+ {\frac{d R}{d s}}{\frac{d \phi}{d s}}{\hat{{\boldsymbol{\phi}}}}+
       R{\frac{d^2 \phi}{d s^2}}{\hat{{\boldsymbol{\phi}}}}- R\left({\frac{d \phi}{d s}}\right)^2 {\hat{{\boldsymbol{R}}}}\end{aligned}

i.e.

.. math::
   :label: eq:kappaline

   \begin{aligned}
   {\boldsymbol{\kappa}}= \left[{\frac{d^2 R}{d s^2}} - R\left({\frac{d \phi}{d s}}\right)^2\right]{\hat{{\boldsymbol{R}}}}+ {\frac{d^2 Z}{d s^2}}{\hat{{\boldsymbol{Z}}}}+
   \left[2{\frac{d R}{d s}}{\frac{d \phi}{d s}} + R{\frac{d^2 \phi}{d s^2}}\right]{\hat{{\boldsymbol{\phi}}}}
   \end{aligned}

Want the components of
`{\boldsymbol{b}}\times{\boldsymbol{\kappa}}`,
and since the vector `{\boldsymbol{b}}` is just the
tangent vector `{\boldsymbol{T}}` above, this can be
written using the cross-products

.. math::

   \begin{aligned}
   {\hat{{\boldsymbol{R}}}}\times{\hat{{\boldsymbol{Z}}}}= -{\hat{{\boldsymbol{\phi}}}}\qquad {\hat{{\boldsymbol{\phi}}}}\times{\hat{{\boldsymbol{Z}}}}= {\hat{{\boldsymbol{R}}}}\qquad
   {\hat{{\boldsymbol{R}}}}\times{\hat{{\boldsymbol{\phi}}}}= {\hat{{\boldsymbol{Z}}}}\end{aligned}

This vector must then be dotted with `\nabla\psi`,
`\nabla\theta`, and `\nabla\phi`. This is done by writing
these vectors in cylindrical coordinates:

.. math::

   \begin{aligned}
   \nabla\psi =& {\frac{\partial \psi}{\partial R}}\hat{{\boldsymbol{R}}} + {\frac{\partial \psi}{\partial Z}}\hat{{\boldsymbol{Z}}} \\ \nabla\theta =&
       \frac{1}{{B_{\text{pol}}}{h_\theta}}\nabla\phi\times\nabla\psi =
       \frac{1}{R{B_{\text{pol}}}{h_\theta}}\left({\frac{\partial \psi}{\partial Z}}\hat{{\boldsymbol{R}}} - {\frac{\partial \psi}{\partial R}}\hat{{\boldsymbol{Z}}}\right) \\\end{aligned}

An alternative is to use

.. math::

   \begin{aligned}
   {\boldsymbol{b}}\times \nabla\phi = \frac{{\sigma_{B\theta}}}{BR^2}\nabla\psi\end{aligned}

and that the tangent vector `{\boldsymbol{T}} =
{\boldsymbol{b}}`. This gives

.. math::
   :label: eq:flinenablapsi

   \begin{aligned}
   \nabla\psi = {\sigma_{B\theta}}BR\left[\frac{dR}{ds}{\boldsymbol{Z}} - \frac{dZ}{ds}{\boldsymbol{R}}\right]
   \end{aligned}

and so because
`d\phi / ds = {B_{\text{tor}}}/ \left(RB\right)`

.. math::
   :label: eq:flinekappsi

   \begin{aligned}
   {\boldsymbol{\kappa}}\cdot\nabla\psi = {\sigma_{B\theta}}BR\left[ \left( \frac{{B_{\text{tor}}}^2}{RB^2} -
   {\frac{d^2 R}{d s^2}}\right){\frac{d Z}{d s}} + {\frac{d^2 Z}{d s^2}}\frac{dR}{ds} \right]
   \end{aligned}

Taking the cross-product of the tangent vector with the curvature in
equation :eq:`eq:kappaline` above gives

.. math::

   \begin{aligned}
     {\boldsymbol{b}}\times{\boldsymbol{\kappa}}=& \left[\frac{{B_{\text{tor}}}}{B}{\frac{d^2 Z}{d s^2}} -
   {\frac{d Z}{d s}}\left(2{\frac{d R}{d s}}{\frac{d \phi}{d s}} + R{\frac{d^2 \phi}{d s^2}}\right)\right]{\boldsymbol{R}} \\ &+
       \left[{\frac{d R}{d s}}\left(2{\frac{d R}{d s}}{\frac{d \phi}{d s}} + R{\frac{d^2 \phi}{d s^2}}\right) -
       \frac{{B_{\text{tor}}}}{B}\left({\frac{d^2 R}{d s^2}} - R\left({\frac{d \phi}{d s}}\right)^2\right)\right]{\boldsymbol{Z}} \\ &+
           \left[{\frac{d Z}{d s}}\left({\frac{d^2 R}{d s^2}} - R\left({\frac{d \phi}{d s}}\right)^2\right) -
           {\frac{d R}{d s}}{\frac{d^2 Z}{d s^2}}\right]{\hat{{\boldsymbol{\phi}}}}\end{aligned}

The components in field-aligned coordinates can then be calculated:

.. math::

   \begin{aligned}
   \left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\right)^x =& {\sigma_{B\theta}}\left({\boldsymbol{b}}\times{\boldsymbol{\kappa}}\right)\cdot\nabla\psi \\ =&
       \frac{R{B_{\text{pol}}}^2}{B}\left(2{\frac{d R}{d s}}{\frac{d \phi}{d s}} + R{\frac{d^2 \phi}{d s^2}}\right) -
       R{B_{\text{tor}}}\left({\frac{d R}{d s}}{\frac{d^2 R}{d s^2}} + {\frac{d Z}{d s}}{\frac{d^2 Z}{d s^2}}\right) +
       \frac{{B_{\text{tor}}}^3}{B^2}{\frac{d R}{d s}}\end{aligned}

Curvature in toroidal coordinates
---------------------------------

In toroidal coordinates `\left(\psi,\theta,\phi\right)`, the
`{\boldsymbol{b}}` vector is

.. math::

   \begin{aligned}
   {\boldsymbol{b}}=& \frac{{B_{\text{pol}}}}{B}{\hat{{\boldsymbol{e}}}}_\theta + \frac{{B_{\text{tor}}}}{B}{\hat{{\boldsymbol{e}}}}_\phi \\ =&
       \frac{{B_{\text{pol}}}{h_\theta}}{B}\nabla\theta + \frac{R{B_{\text{tor}}}}{B}\nabla\phi\end{aligned}

The curl of this vector is

.. math::

   \begin{aligned}
   \left(\nabla\times{\boldsymbol{b}}\right)^\psi =& \frac{1}{\sqrt{g}}\left({\frac{\partial b_\phi}{\partial \theta}} -
       {\frac{\partial b_\theta}{\partial \phi}}\right) \\ \left(\nabla\times{\boldsymbol{b}}\right)^\theta =&
       \frac{1}{\sqrt{g}}\left({\frac{\partial b_\psi}{\partial \phi}} - {\frac{\partial b_\phi}{\partial \psi}}\right) \\
       \left(\nabla\times{\boldsymbol{b}}\right)^\phi =& \frac{1}{\sqrt{g}}\left({\frac{\partial b_\theta}{\partial \psi}}
       - {\frac{\partial b_\psi}{\partial \theta}}\right)\end{aligned}

where
`1/\sqrt{g} = {B_{\text{pol}}}/{h_\theta}`.
Therefore, in terms of unit vectors:

.. math::

   \begin{aligned}
   \nabla\times{\boldsymbol{b}}=
   \frac{1}{R{h_\theta}}{\frac{\partial }{\partial \theta}}\left(\frac{R{B_{\text{tor}}}}{B}\right){\hat{{\boldsymbol{e}}}}_\psi -
   {B_{\text{pol}}}{\frac{\partial }{\partial \psi}}\left(\frac{R{B_{\text{tor}}}}{B}\right){\hat{{\boldsymbol{e}}}}_\theta + \frac{{B_{\text{pol}}}
   R}{{h_\theta}}{\frac{\partial }{\partial \psi}}\left(\frac{{h_\theta}{B_{\text{pol}}}}{B}\right){\hat{{\boldsymbol{e}}}}_\phi\end{aligned}

psi derivative of the B field
-----------------------------

Needed to calculate magnetic shear, and one way to get the curvature.
The simplest way is to use finite differencing, but there is another way
using local derivatives (implemented using DCT).

.. math::

   \begin{aligned}
   {B_{\text{pol}}}= \frac{\left|\nabla\psi\right|}{R} = \frac{1}{R}\sqrt{\left({\frac{\partial \psi}{\partial R}}\right)^2 +
   \left({\frac{\partial \psi}{\partial R}}\right)^2}\end{aligned}

Using

.. math::

   \begin{aligned}
   \nabla{B_{\text{pol}}}= {\frac{\partial {B_{\text{pol}}}}{\partial \psi}}\nabla\psi + {\frac{\partial {B_{\text{pol}}}}{\partial \theta}}\nabla\theta +
   {\frac{\partial {B_{\text{pol}}}}{\partial \phi}}\nabla\phi\end{aligned}

we get

.. math::

   \begin{aligned}
   \nabla{B_{\text{pol}}}\cdot\nabla\psi = {\frac{\partial {B_{\text{pol}}}}{\partial \psi}}\left|\nabla\psi\right|^2\end{aligned}

and so

.. math::

   \begin{aligned}
   {\frac{\partial {B_{\text{pol}}}}{\partial \psi}} = \nabla{B_{\text{pol}}}\cdot\nabla\psi / \left(R{B_{\text{pol}}}\right)^2\end{aligned}

The derivatives of `{B_{\text{pol}}}` in `R` and
`Z` are:

.. math::

   \begin{aligned}
   {\frac{\partial {B_{\text{pol}}}}{\partial R}} =& -\frac{{B_{\text{pol}}}}{R} + \frac{1}{{B_{\text{pol}}}
   R^2}\left[{\frac{\partial \psi}{\partial R}}{\frac{\partial^2 \psi}{\partial {R}^2}} +
   {\frac{\partial \psi}{\partial Z}}\frac{\partial^2\psi}{\partial R\partial Z}\right] \\ {\frac{\partial {B_{\text{pol}}}}{\partial Z}}
   =& \frac{1}{{B_{\text{pol}}}R^2}\left[{\frac{\partial \psi}{\partial Z}}{\frac{\partial^2 \psi}{\partial {Z}^2}} +
   {\frac{\partial \psi}{\partial R}}\frac{\partial^2\psi}{\partial R\partial Z}\right]\end{aligned}

For the toroidal field, `{B_{\text{tor}}}= f/R`

.. math::

   \begin{aligned}
   {\frac{\partial {B_{\text{tor}}}}{\partial \psi}} = \frac{1}{R}{\frac{\partial f}{\partial \psi}} - \frac{f}{R^2}{\frac{\partial R}{\partial \psi}}\end{aligned}

As above,
`{\frac{\partial R}{\partial \psi}} = \nabla R \cdot\nabla\psi / \left(R{B_{\text{pol}}}\right)^2`,
and since `\nabla R\cdot\nabla R = 1`,

.. math::

   \begin{aligned}
   {\frac{\partial R}{\partial \psi}} = {\frac{\partial \psi}{\partial R}} / \left(R{B_{\text{pol}}}\right)^2\end{aligned}

similarly,

.. math::

   \begin{aligned}
   {\frac{\partial Z}{\partial \psi}} = {\frac{\partial \psi}{\partial Z}} / \left(R{B_{\text{pol}}}\right)^2\end{aligned}

and so the variation of toroidal field with `\psi` is

.. math::

   \begin{aligned}
   {\frac{\partial {B_{\text{tor}}}}{\partial \psi}} = \frac{1}{R}{\frac{\partial f}{\partial \psi}} -
   \frac{{B_{\text{tor}}}}{R^3{B_{\text{pol}}}^2}{\frac{\partial \psi}{\partial R}}\end{aligned}

From the definition
`B=\sqrt{{B_{\text{tor}}}^2 + {B_{\text{pol}}}^2}`,

.. math::

   \begin{aligned}
   {\frac{\partial B}{\partial \psi}} = \frac{1}{B}\left({B_{\text{tor}}}{\frac{\partial {B_{\text{tor}}}}{\partial \psi}} + {B_{\text{pol}}}{\frac{\partial {B_{\text{pol}}}}{\partial \psi}}\right)\end{aligned}

Parallel derivative of the B field
----------------------------------

To get the parallel nablaients of the `B` field components, start
with

.. math::

   \begin{aligned}
   {\frac{\partial }{\partial s}}\left(B^2\right) = {\frac{\partial }{\partial s}}\left({B_{\text{tor}}}^2\right) + {\frac{\partial }{\partial s}}\left({B_{\text{pol}}}^2\right)\end{aligned}

Using the fact that `R{B_{\text{tor}}}` is constant
along `s`,

.. math::

   \begin{aligned}
   {\frac{\partial }{\partial s}}\left(R^2{B_{\text{tor}}}^2\right) = R^2{\frac{\partial }{\partial s}}\left({B_{\text{tor}}}^2\right) +
   {B_{\text{tor}}}^2{\frac{\partial }{\partial s}}\left(R^2\right) = 0\end{aligned}

which gives

.. math::

   \begin{aligned}
     {\frac{\partial }{\partial s}}\left({B_{\text{tor}}}^2\right) = -\frac{{B_{\text{tor}}}^2}{R^2}{\frac{\partial }{\partial s}}\left(R^2\right)\end{aligned}

The poloidal field can be calculated from

.. math::

   \begin{aligned}
   {\frac{\partial }{\partial s}}\left(\nabla\psi \cdot \nabla\psi\right) = {\frac{\partial }{\partial s}}\left(R^2{B_{\text{pol}}}^2\right) =
   R^2{\frac{\partial }{\partial s}}\left({B_{\text{pol}}}^2\right) + {B_{\text{pol}}}^2{\frac{\partial }{\partial s}}\left(R^2\right)\end{aligned}

Using equation :eq:`eq:flinenablapsi`,
`\nabla\psi \cdot \nabla\psi` can also be written as

.. math::

   \begin{aligned}
   \nabla\psi \cdot \nabla\psi = B^2R^2\left[\left({\frac{\partial R}{\partial s}}\right)^2 +
   \left({\frac{\partial Z}{\partial s}}\right)^2\right]\end{aligned}

and so (unsurprisingly)

.. math::

   \begin{aligned}
   \frac{{B_{\text{pol}}}^2}{B^2} = \left[\left({\frac{\partial R}{\partial s}}\right)^2 + \left({\frac{\partial Z}{\partial s}}\right)^2\right]\end{aligned}

Hence

.. math::

   \begin{aligned}
   {\frac{\partial }{\partial s}}\left({B_{\text{pol}}}^2\right) = B^2{\frac{\partial }{\partial s}}\left[\left({\frac{\partial R}{\partial s}}\right)^2 +
   \left({\frac{\partial Z}{\partial s}}\right)^2\right] + \frac{{B_{\text{pol}}}^2}{B^2}{\frac{\partial }{\partial s}}\left(B^2\right)\end{aligned}

Which gives

.. math::

   \begin{aligned}
   {\frac{\partial }{\partial s}}\left(B^2\right) = -\frac{B^2}{R^2}{\frac{\partial }{\partial s}}\left(R^2\right) +
   \frac{B^4}{{B_{\text{tor}}}^2}{\frac{\partial }{\partial s}}\left[\left({\frac{\partial R}{\partial s}}\right)^2 + \left({\frac{\partial Z}{\partial s}}\right)^2\right]\end{aligned}

.. math::

   \begin{aligned}
   {\frac{\partial }{\partial s}}\left({B_{\text{pol}}}^2\right) = \left(1 +
   \frac{{B_{\text{pol}}}^2}{{B_{\text{tor}}}^2}\right)B^2{\frac{\partial }{\partial s}}\left[\left({\frac{\partial R}{\partial s}}\right)^2 +
   \left({\frac{\partial Z}{\partial s}}\right)^2\right] - \frac{{B_{\text{pol}}}^2}{R^2}{\frac{\partial }{\partial s}}\left(R^2\right)\end{aligned}

Magnetic shear from J x B
-------------------------

Re-arranging the radial force balance
equation :eq:`eq:xbalance` gives

.. math::

   \begin{aligned}
   \frac{{B_{\text{pol}}}^2R}{{B_{\text{tor}}}}{\frac{\partial \nu}{\partial \psi}} + \nu\left(\frac{2RB}{{B_{\text{tor}}}}{\frac{\partial B}{\partial \psi}} +
   \frac{B^2}{{B_{\text{tor}}}}{\frac{\partial R}{\partial \psi}} - \frac{B^2R}{{B_{\text{tor}}}^2}{\frac{\partial {B_{\text{tor}}}}{\partial \psi}}\right) +
   \frac{\mu_0{h_\theta}}{{B_{\text{pol}}}}{\frac{\partial P}{\partial \psi}} = 0\end{aligned}

Magnetic shear
--------------

The field-line pitch is given by

.. math::

   \begin{aligned}
   \nu = \frac{{h_\theta}{B_{\text{tor}}}}{{B_{\text{pol}}}R}\end{aligned}

and so

.. math::

   \begin{aligned}
   {\frac{\partial \nu}{\partial \psi}} = \frac{\nu}{{h_\theta}}{\frac{\partial {h_\theta}}{\partial \psi}} +
   \frac{\nu}{{B_{\text{tor}}}}{\frac{\partial {B_{\text{tor}}}}{\partial \psi}} - \frac{\nu}{{B_{\text{pol}}}}{\frac{\partial {B_{\text{pol}}}}{\partial \psi}} -
   \frac{\nu}{R}{\frac{\partial R}{\partial \psi}}\end{aligned}

The last three terms are given in the previous section, but
`\partial{h_\theta}/\partial\psi` needs to be evaluated

psi derivative of h
-------------------

From the expression for curvature (equation :eq:`eq:curvature`),
and using
`\nabla x \cdot \nabla \psi = {\sigma_{B\theta}}\left(R{B_{\text{pol}}}\right)^2`
and
`\nabla z\cdot\nabla \psi = -{\sigma_{B\theta}}I \left(R{B_{\text{pol}}}\right)^2`

.. math::

   \begin{aligned}
   {\boldsymbol{\kappa}}\cdot\nabla\psi =& -{\sigma_{B\theta}}
       \frac{{B_{\text{pol}}}}{B{h_\theta}}{\left({R{B_{\text{pol}}}}\right)^2}\left[{\frac{\partial }{\partial x}}\left(\frac{B{h_\theta}}{{B_{\text{pol}}}}\right) -
       {\sigma_{B\theta}}{\frac{\partial }{\partial y}}\left(\frac{{B_{\text{tor}}}IR}{B}\right)\right] \\ &- I{\left({R{B_{\text{pol}}}}\right)^2}
           \frac{{B_{\text{pol}}}}{B{h_\theta}}{\frac{\partial }{\partial y}}\left(\frac{{B_{\text{tor}}}R}{B}\right)\end{aligned}

The second and third terms partly cancel, and using
`{\frac{\partial I}{\partial y}} = {\sigma_{B\theta}}
{\frac{\partial \nu}{\partial x}}`

.. math::

   \begin{aligned}
     \frac{{\boldsymbol{\kappa}}\cdot\nabla\psi}{{\left({R{B_{\text{pol}}}}\right)^2}} =&
       -{\sigma_{B\theta}}\frac{{B_{\text{pol}}}}{B{h_\theta}}{\frac{\partial }{\partial x}}\left(\frac{B{h_\theta}}{{B_{\text{pol}}}}\right) +
       {\sigma_{B\theta}}\frac{{B_{\text{pol}}}}{B{h_\theta}}\frac{{B_{\text{tor}}}R}{B}{\frac{\partial \nu}{\partial x}} \\ =&
       -{\sigma_{B\theta}}\frac{{B_{\text{pol}}}}{B{h_\theta}}\left[{\frac{\partial }{\partial x}}\left(\frac{B{h_\theta}}{{B_{\text{pol}}}}\right) - \frac{{B_{\text{tor}}}
       R}{B}{\frac{\partial }{\partial x}}\left(\frac{{B_{\text{tor}}}{h_\theta}}{{B_{\text{pol}}}R}\right)\right] \\ =&
               -{\sigma_{B\theta}}\frac{{B_{\text{pol}}}}{B{h_\theta}}\left[{h_\theta}{\frac{\partial }{\partial x}}\left(\frac{B}{{B_{\text{pol}}}}\right) -
               {h_\theta}\frac{{B_{\text{tor}}}R}{B}{\frac{\partial }{\partial x}}\left(\frac{{B_{\text{tor}}}}{{B_{\text{pol}}}R}\right) +
           \frac{B^2}{B{B_{\text{pol}}}}{\frac{\partial {h_\theta}}{\partial x}} -
       \frac{{B_{\text{tor}}}^2}{B{B_{\text{pol}}}}{\frac{\partial {h_\theta}}{\partial x}}\right] \\ =& -{\sigma_{B\theta}}
           \frac{{B_{\text{pol}}}}{B^2{h_\theta}}{\frac{\partial {h_\theta}}{\partial x}} -
           {\sigma_{B\theta}}\frac{{B_{\text{pol}}}}{B^2}\left[B{\frac{\partial }{\partial x}}\left(\frac{B}{{B_{\text{pol}}}}\right) - {B_{\text{tor}}}
           R{\frac{\partial }{\partial x}}\left(\frac{{B_{\text{tor}}}}{{B_{\text{pol}}}R}\right)\right]\end{aligned}

Writing

.. math::

   \begin{aligned}
   B{\frac{\partial }{\partial x}}\left(\frac{B}{{B_{\text{pol}}}}\right) =& {\frac{\partial }{\partial x}}\left(\frac{B^2}{{B_{\text{pol}}}}\right) -
       \frac{B}{{B_{\text{pol}}}}{\frac{\partial B}{\partial x}} \\ {B_{\text{tor}}}R{\frac{\partial }{\partial x}}\left(\frac{{B_{\text{tor}}}}{{B_{\text{pol}}}R}\right) =&
       {\frac{\partial }{\partial x}}\left(\frac{{B_{\text{tor}}}^2}{{B_{\text{pol}}}}\right) - \frac{{B_{\text{tor}}}}{{B_{\text{pol}}}R}{\frac{\partial }{\partial x}}\left({B_{\text{tor}}}
       R\right)\end{aligned}

and using
`B{\frac{\partial B}{\partial x}} = {B_{\text{tor}}}{\frac{\partial {B_{\text{tor}}}}{\partial x}} + {B_{\text{pol}}}{\frac{\partial {B_{\text{pol}}}}{\partial x}}`,
this simplifies to give

.. math::
   :label: eq:dhdpsi

   \begin{aligned}
   \frac{{\boldsymbol{\kappa}}\cdot\nabla\psi}{{\left({R{B_{\text{pol}}}}\right)^2}} =
   -{\sigma_{B\theta}}\frac{{B_{\text{pol}}}^2}{B^2{h_\theta}}{\frac{\partial {h_\theta}}{\partial x}} - {\sigma_{B\theta}}\frac{{B_{\text{tor}}}^2}{B^2
   R}{\frac{\partial R}{\partial x}}
   \end{aligned}

This can be transformed into an expression for
`{\frac{\partial {h_\theta}}{\partial x}}`
involving only derivatives along field-lines. Writing `\nabla R =
{\frac{\partial R}{\partial \psi}}\nabla\psi + {\frac{\partial R}{\partial \theta}}\nabla\theta`,

.. math::

   \begin{aligned}
   \nabla R \cdot \nabla\psi = {\frac{\partial R}{\partial \psi}}{\left({R{B_{\text{pol}}}}\right)^2}\end{aligned}

Using :eq:`eq:flinenablapsi`,

.. math::

   \begin{aligned}
   \nabla\psi \cdot \nabla R = -{\sigma_{B\theta}}B R\frac{dZ}{ds}\end{aligned}

and so

.. math::

   \begin{aligned}
   {\frac{\partial R}{\partial x}} = -\frac{BR}{{\left({R{B_{\text{pol}}}}\right)^2}}\frac{dZ}{ds}\end{aligned}

Substituting this and equation :eq:`eq:flinekappsi`
for `{\boldsymbol{\kappa}}\cdot\nabla\psi` into
equation :eq:`eq:dhdpsi` the
`{\frac{\partial R}{\partial x}}` term cancels with
part of the `{\boldsymbol{\kappa}}\cdot\nabla\psi`
term, simplifying to

.. math::

   \begin{aligned}
   {\frac{\partial {h_\theta}}{\partial x}} =
   -{h_\theta}\frac{B^3R}{{B_{\text{pol}}}^2{\left({R{B_{\text{pol}}}}\right)^2}}\left[\frac{d^2Z}{ds^2}\frac{dR}{ds} -
   \frac{d^2R}{ds^2}\frac{dZ}{ds}\right]\end{aligned}

.. _sec:shiftcoords:

Shifted radial derivatives
==========================

The coordinate system given by
equation :eq:`eq:coordtransform` and used in the
above sections has a problem: There is a special poloidal location
`\theta_0` where the radial basis vector
`{\boldsymbol{e}}_x` is purely in the
`\nabla\psi` direction. Moving away from this location, the
coordinate system becomes sheared in the toroidal direction.

Making the substitution

.. math::

   \begin{aligned}
   {\frac{\partial }{\partial x}} = {\frac{\partial }{\partial \psi}} + I{\frac{\partial }{\partial z}}\end{aligned}

we also get the mixed derivative

.. math::

   \begin{aligned}
   \frac{\partial}{\partial z\partial x} =& {\frac{\partial }{\partial z}}{\frac{\partial }{\partial \psi}} +
       {\frac{\partial I}{\partial z}}{\frac{\partial }{\partial z}} + I\frac{\partial^2}{\partial z^2} \nonumber \\ =&
       \frac{\partial^2}{\partial z\partial \psi} + I\frac{\partial^2}{\partial
       z^2}\end{aligned}

and second-order `x` derivative

.. math::

   \begin{aligned}
   \frac{\partial^2}{\partial x^2} =& \frac{\partial^2}{\partial \psi^2} +
       {\frac{\partial }{\partial \psi}}\left(I{\frac{\partial }{\partial z}}\right) + I{\frac{\partial }{\partial z}}\left({\frac{\partial }{\partial \psi}} +
       I{\frac{\partial }{\partial z}}\right) \nonumber \\ =& \frac{\partial^2}{\partial \psi^2} +
       I^2\frac{\partial^2}{\partial z^2} + 2I\frac{\partial^2}{\partial z\partial
       \psi} + {\frac{\partial I}{\partial \psi}}{\frac{\partial }{\partial z}}\end{aligned}

Perpendicular Laplacian
-----------------------

.. math::

   \begin{aligned}
   \nabla_\perp^2= {\left({R{B_{\text{pol}}}}\right)^2}\left[{\frac{\partial^2 }{\partial {x}^2}} - 2I\frac{\partial^2}{\partial z\partial x} +
   \left(I^2 + \frac{B^2}{\left({R{B_{\text{pol}}}}\right)^4}\right){\frac{\partial^2 }{\partial {z}^2}}\right]\end{aligned}

transforms to

.. math::
   :label: eq:delp

   \begin{aligned}
   \nabla_\perp^2= {\left({R{B_{\text{pol}}}}\right)^2}\left[{\frac{\partial^2 }{\partial {\psi}^2}} + {\frac{\partial I}{\partial \psi}}{\frac{\partial }{\partial z}} +
   \frac{B^2}{\left({R{B_{\text{pol}}}}\right)^4}{\frac{\partial^2 }{\partial {z}^2}}\right]
   \end{aligned}

The extra term involving `I` disappears, but only if both the
`x` and `z` first derivatives are taken into account:

.. math::

   \begin{aligned}
   \nabla_\perp^2= {\left({R{B_{\text{pol}}}}\right)^2}\left[{\frac{\partial^2 }{\partial {x}^2}} - 2I\frac{\partial^2}{\partial z\partial x} +
   \left(I^2 + \frac{B^2}{\left({R{B_{\text{pol}}}}\right)^4}\right){\frac{\partial^2 }{\partial {z}^2}}\right] + \nabla^2 x {\frac{\partial }{\partial x}} +
   \nabla^2 z{\frac{\partial }{\partial z}}\end{aligned}

with

.. math::

   \begin{aligned}
   \nabla^2 x = \frac{1}{J}{\frac{\partial }{\partial x}}\left[J{\left({R{B_{\text{pol}}}}\right)^2}\right]\end{aligned}

.. math::
   :label: eq:delpz

   \begin{aligned}
   \nabla^2 z =& \frac{1}{J}\left[-{\frac{\partial }{\partial x}}\left(JI{\left({R{B_{\text{pol}}}}\right)^2}\right) -
   {\frac{\partial }{\partial y}}\left(\frac{{B_{\text{tor}}}}{{B_{\text{pol}}}^2R}\right)\right] \nonumber \\ =&
       \frac{1}{J}\left[-I{\frac{\partial }{\partial x}}\left(J{\left({R{B_{\text{pol}}}}\right)^2}\right) - {\frac{\partial I}{\partial x}}J{\left({R{B_{\text{pol}}}}\right)^2}-
       {\frac{\partial }{\partial y}}\left(\frac{{B_{\text{tor}}}}{{B_{\text{pol}}}^2R}\right)\right] \end{aligned}

where `J={h_\theta}/ {B_{\text{pol}}}` is
the Jacobian. Transforming into `\psi` derivatives, the middle
term of equation :eq:`eq:delpz` cancels the `I` term
in equation :eq:`eq:delp`, but introduces another `I`
term (first term in equation :eq:`eq:delpz`). This term
cancels with the `\nabla^2 x` term when
`{\frac{\partial }{\partial x}}` is expanded, so the
full expression for `\nabla_\perp^2` using `\psi`
derivatives is:

.. math::
   :label: eq:delp_shift

   \begin{aligned}
   \nabla_\perp^2=& {\left({R{B_{\text{pol}}}}\right)^2}\left[{\frac{\partial^2 }{\partial {\psi}^2}} + \frac{B^2}{\left({R{B_{\text{pol}}}}\right)^4}{\frac{\partial^2 }{\partial {z}^2}}\right]
       \nonumber \\ &+ \frac{1}{J}{\frac{\partial }{\partial \psi}}\left[J{\left({R{B_{\text{pol}}}}\right)^2}\right]{\frac{\partial }{\partial \psi}} -
       \frac{1}{J}{\frac{\partial }{\partial y}}\left(\frac{{B_{\text{tor}}}}{{B_{\text{pol}}}^2R}\right){\frac{\partial }{\partial z}}
   \end{aligned}

In orthogonal (psi, theta, zeta) flux coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For comparison, the perpendicular Laplacian can be derived in orthogonal
“flux” coordinates

.. math::

   \begin{aligned}
   \left|\nabla\psi\right| = {R{B_{\text{pol}}}}\qquad \left|\nabla\theta\right| = 1/{h_\theta}\qquad
   \left|\nabla\zeta\right| = 1/R\end{aligned}

The Laplacian operator is given by

.. math::

   \begin{aligned}
   \nabla^2 A =& {\left({R{B_{\text{pol}}}}\right)^2}{\frac{\partial^2 A}{\partial {\psi}^2}} + \frac{1}{{h_\theta}^2}{\frac{\partial^2 A}{\partial {\theta}^2}} +
       \frac{1}{R^2}{\frac{\partial^2 A}{\partial {\zeta}^2}} \nonumber \\ &+
       \frac{1}{J}{\frac{\partial }{\partial \psi}}\left[J{\left({R{B_{\text{pol}}}}\right)^2}\right]{\frac{\partial A}{\partial \psi}} +
       \frac{1}{J}{\frac{\partial }{\partial \theta}}\left(J/{h_\theta}^2\right){\frac{\partial A}{\partial \theta}}\end{aligned}

parallel derivative by

.. math::

   \begin{aligned}
   \partial_{||} \equiv {\boldsymbol{b}}\cdot\nabla = \frac{{B_{\text{pol}}}}{B{h_\theta}}{\frac{\partial }{\partial \theta}} +
   \frac{{B_{\text{tor}}}}{RB}{\frac{\partial }{\partial \zeta}}\end{aligned}

and so

.. math::

   \begin{aligned}
   \partial^2_{||}A \equiv \partial_{||}\left(\partial_{||}A\right) =&
       \left(\frac{{B_{\text{pol}}}}{B{h_\theta}}\right)^2{\frac{\partial^2 A}{\partial {\theta}^2}} +
       \left(\frac{{B_{\text{tor}}}}{RB}\right)^2{\frac{\partial^2 A}{\partial {\zeta}^2}} \nonumber \\ &+
       2\frac{{B_{\text{pol}}}{B_{\text{tor}}}}{B^2{h_\theta}R}\frac{\partial^2 A}{\partial\theta\partial\zeta}
       \nonumber \\ &+ {\frac{\partial }{\partial \theta}}\left(\frac{{B_{\text{pol}}}}{B{h_\theta}}\right){\frac{\partial A}{\partial \theta}} +
       {\frac{\partial }{\partial \theta}}\left(\frac{{B_{\text{tor}}}}{RB}\right){\frac{\partial A}{\partial \zeta}}\end{aligned}

Hence in orthogonal flux coordinates, the perpendicular Laplacian is:

.. math::
   :label: eq:delp_flux

   \begin{aligned}
   \nabla_\perp^2\equiv \nabla^2 - \partial_{||}^2 = {\left({R{B_{\text{pol}}}}\right)^2}\left[{\frac{\partial^2 }{\partial {\psi}^2}} +
   \frac{1}{R^4B^2}{\frac{\partial^2 }{\partial {\zeta^2}^2}}\right] +
   \frac{{B_{\text{tor}}}^2}{{h_\theta}^2B^2}{\frac{\partial^2 }{\partial {\theta}^2}} + \cdots
   \end{aligned}

where the neglected terms are first-order derivatives. The coefficient
for the second-order `z` derivative differs from
equation :eq:`eq:delp_shift`, and
equation :eq:`eq:delp_flux` still contains a
derivative in `\theta`. This shows that the transformation made to
get equation :eq:`eq:delp_shift` doesn’t result in
the same answer as orthogonal flux coordinates:
equation :eq:`eq:delp_shift` is in field-aligned
coordinates.

Note that in the limit of `{B_{\text{pol}}}= B`, both equations
:eq:`eq:delp_shift` and :eq:`eq:delp_flux` are the same, as they
should be.

Operator B x Nabla Phi Dot Nabla A
----------------------------------

.. math::

   \begin{aligned}
   {\boldsymbol{B}}\times\nabla\phi\cdot\nabla A =& \left({\frac{\partial \phi}{\partial x}}{\frac{\partial A}{\partial y}} -
       {\frac{\partial \phi}{\partial y}}{\frac{\partial A}{\partial x}}\right)\left(-{B_{\text{tor}}}\frac{{R{B_{\text{pol}}}}}{{h_\theta}}\right) \\ &+
       \left({\frac{\partial \phi}{\partial x}}{\frac{\partial A}{\partial z}} - {\frac{\partial \phi}{\partial z}}{\frac{\partial A}{\partial x}}\right)\left(-B^2\right)
       \\ &- \left({\frac{\partial \phi}{\partial y}}{\frac{\partial A}{\partial z}} -
       {\frac{\partial \phi}{\partial z}}{\frac{\partial A}{\partial y}}\right)\left(I{B_{\text{tor}}}\frac{{R{B_{\text{pol}}}}}{{h_\theta}}\right)\end{aligned}

.. math::

   \begin{aligned}
   {\boldsymbol{B}}\times\nabla\phi\cdot\nabla A =& \left({\frac{\partial \phi}{\partial \psi}}{\frac{\partial A}{\partial y}} + I
       {\frac{\partial \phi}{\partial z}}{\frac{\partial A}{\partial y}} - {\frac{\partial \phi}{\partial y}}{\frac{\partial A}{\partial \psi}} -
       I{\frac{\partial \phi}{\partial y}}{\frac{\partial A}{\partial z}}\right)\left(-{B_{\text{tor}}}\frac{{R{B_{\text{pol}}}}}{{h_\theta}}\right) \\ &+
       \left({\frac{\partial \phi}{\partial \psi}}{\frac{\partial A}{\partial z}} + I{\frac{\partial \phi}{\partial z}}{\frac{\partial A}{\partial z}} -
       {\frac{\partial \phi}{\partial z}}{\frac{\partial A}{\partial \psi}} - I{\frac{\partial \phi}{\partial z}}{\frac{\partial A}{\partial z}}\right)\left(-B^2\right)
       \\ &- \left({\frac{\partial \phi}{\partial y}}{\frac{\partial A}{\partial z}} -
       {\frac{\partial \phi}{\partial z}}{\frac{\partial A}{\partial y}}\right)\left(I{B_{\text{tor}}}\frac{{R{B_{\text{pol}}}}}{{h_\theta}}\right)\end{aligned}

.. math::

   \begin{aligned}
   {\boldsymbol{B}}\times\nabla\phi\cdot\nabla A =& \left({\frac{\partial \phi}{\partial \psi}}{\frac{\partial A}{\partial y}} -
       {\frac{\partial \phi}{\partial y}}{\frac{\partial A}{\partial \psi}}\right)\left(-{B_{\text{tor}}}\frac{{R{B_{\text{pol}}}}}{{h_\theta}}\right) \nonumber \\
       &+ \left({\frac{\partial \phi}{\partial \psi}}{\frac{\partial A}{\partial z}} - {\frac{\partial \phi}{\partial z}}{\frac{\partial A}{\partial \psi}}
       \right)\left(-B^2\right)\end{aligned}

Useful identities
=================

`\mathbf{b}\times\mathbf{\kappa}\cdot\nabla\psi \simeq -RB_\zeta\partial_{||}\ln B`
-----------------------------------------------------------------------------------------

Using
`\mathbf{b}\times\mathbf{\kappa} \simeq \frac{B}{2}\nabla\times\frac{\mathbf{b}}{B}`,
and working in orthogonal `\left(\psi, \theta, \zeta\right)`
coordinates. The magnetic field unit vector is:

.. math:: \mathbf{b} = \frac{B_\theta h_\theta}{B}\nabla\theta + \frac{B_\zeta R}{B}\nabla\zeta

and using the definition of curl (equation :eq:`eq:curlcurvilinear`)
we can write

.. math:: \mathbf{b}\times\mathbf{\kappa} \simeq \frac{B}{2}\nabla\times\frac{\mathbf{b}}{B} = \frac{B}{2}\frac{B_\theta}{h_\theta}\left[\frac{\partial}{\partial\theta}\left(\frac{B_\zeta R}{B^2}\right) - \frac{\partial}{\partial\zeta}\left(\frac{B_\theta h_\theta}{B^2}\right)\right]\mathbf{e}_\psi + \left[\cdot\right]\mathbf{e}_\theta + \left[\cdot\right]\mathbf{e}_\zeta

so that when dotted with `\nabla\psi`, only the first bracket
survives. The parallel gradient is

.. math:: \partial_{||} = \mathbf{b}\cdot\nabla = \frac{B_\theta}{Bh_\theta}\frac{\partial}{\partial\theta} + \frac{B_\theta}{BR}\frac{\partial}{\partial\zeta}

Neglecting derivatives for axisymmetric equilibrium

.. math:: \frac{B}{2}\nabla\times\frac{\mathbf{b}}{B}\cdot\nabla\psi = \frac{B}{2}B\partial_{||}\left(\frac{B_\zeta R}{B^2}\right)

Since `B_\zeta R` is a flux function, this can be written as

.. math:: \frac{B}{2}\nabla\times\frac{\mathbf{b}}{B}\cdot\nabla\psi = -B_\zeta R\frac{1}{B}\partial_{||} B

and so

.. math:: \mathbf{b}\times\mathbf{\kappa}\cdot\nabla\psi \simeq -RB_\zeta\partial_{||}\ln B

.. raw:: latex

   \bibliographystyle{unsrt}

.. raw:: latex

   \appendix

Differential geometry
=====================

.. warning:: Several mistakes have been found (and is now corrected)
  in this section, so it should be proof read before removing this
  warning!  The following are notes from [haeseler]_.

Sets of vectors `\left\{\mathbf{A, B, C}\right\}` and
`\left\{\mathbf{a, b, c}\right\}` are reciprocal if

.. math::

   \begin{aligned}
   \mathbf{A\cdot a} = \mathbf{B\cdot b} = \mathbf{C\cdot c} = 1\\ \mathbf{A\cdot
   b} = \mathbf{A\cdot c} = \mathbf{B\cdot a} = \mathbf{B\cdot c} = \mathbf{C\cdot
   a} = \mathbf{C\cdot b} = 0 \\\end{aligned}

which implies that `\left\{\mathbf{A, B, C}\right\}` and
`\left\{\mathbf{a, b, c}\right\}` are each linearly independent.
Equivalently,

.. math::

   \begin{aligned}
   \mathbf{a} = \frac{\mathbf{B\times C}}{\mathbf{A\cdot\left(B\times C\right)}}\qquad
   {\boldsymbol{b}}= \frac{\mathbf{C\times A}}{\mathbf{B\cdot\left(C\times A\right)}}\qquad
   \mathbf{c} = \frac{\mathbf{A\times B}}{\mathbf{C\cdot\left(A\times B\right)}}\end{aligned}

Either of these sets can be used as a basis, and any vector
`\mathbf{w}` can be represented as
`\mathbf{w} = \left(\mathbf{w\cdot a}\right)\mathbf{A} +
\left(\mathbf{w\cdot b}\right){\boldsymbol{B}}+ \left(\mathbf{w\cdot c}\right)\mathbf{C}`
or as
`\mathbf{w} = \left(\mathbf{w\cdot A}\right)\mathbf{a} + \left(\mathbf{w\cdot B}\right){\boldsymbol{b}}
+ \left(\mathbf{w\cdot C}\right)\mathbf{c}`. In the Cartesian coordinate
system, the basis vectors are reciprocal to themselves so this
distinction is not needed. For a general set of coordinates
`\left\{u^1, u^2, u^3\right\}`, tangent basis vectors can be
defined. If the Cartesian coordinates of a point are given by
`\left(x, y, z\right) = \mathbf{R}\left(u^1, u^2, u^3\right)` then
the tangent basis vectors are:

.. math::

   \begin{aligned}
   {\boldsymbol{e}}_i = \frac{\partial\mathbf{R}}{\partial u^i}\end{aligned}

and in general these will vary from point to point. The scale factor or
metric coefficient
`h_i =\left|{\boldsymbol{e}}_i\right|` is the distance
moved for a unit change in `u^i`. The unit vector
`\hat{{\boldsymbol{e}}}_i = {\boldsymbol{e}}_i/h_i`.
Definition of nabla operator:

.. raw:: latex

   \framebox{$\nabla\Phi$ of a function $\Phi$ is defined so that $d\Phi =
   \nabla\Phi\cdot d{\mathbf{R}}$}

From the chain rule,
`d\mathbf{R} = \frac{\partial\mathbf{R}}{\partial u^i}du^i
= {\boldsymbol{e}}_idu^i` and substituting `\Phi = u^i`

.. math::

   \begin{aligned}
   du^i = \nabla u^i\cdot{\boldsymbol{e}}_jdu^j\end{aligned}

which can only be true if
`\nabla u^i\cdot{\boldsymbol{e}}_j = \delta^i_j` i.e.
if

.. math::

   \text{Sets of vectors $\boldsymbol{e}^i\equiv\nabla u^i$ and
   $\boldsymbol{e}_j$ are reciprocal}

Since the sets of vectors
`\left\{{\boldsymbol{e}}^i\right\}` and
`\left\{{\boldsymbol{e}}_i\right\}` are reciprocal, any
vector `\mathbf{D}` can be written as
`\mathbf{D} = D_i{\boldsymbol{e}}^i
= D^i{\boldsymbol{e}}_i` where
`D_i = \mathbf{D\cdot e}_i` are the covariant components and
`D^i = \mathbf{D\cdot e}^i` are the contravariant components. To
convert between covariant and contravariant components, define the
metric coefficients `g_{ij} = \mathbf{e_i\cdot e_j}` and
`g^{ij} =
\mathbf{e^i\cdot e^j}` so that
`{\boldsymbol{e}}_i = g_{ij}{\boldsymbol{e}}^j`.
`g_{ij}` and `g^{ij}` are symmetric and if the basis is
orthogonal then `g_{ij}=g^{ij} = 0` for `i\neq j` i.e. the
metric is diagonal.

.. math::

   \text{$g_{ij} = h_ih_j\boldsymbol{e}_i\cdot\boldsymbol{e}_j$ and so $g_{ii} = h_i^2$}

For a general set of coordinates, the nabla operator can be expressed as

.. math::

   \begin{aligned}
   \nabla = \nabla u^i\frac{\partial}{\partial u^i} =
   {\boldsymbol{e}}^i\frac{\partial}{\partial u^i}\end{aligned}

and for a general set of (differentiable) coordinates
`\left\{u^i\right\}`, the Laplacian is given by

.. math::
   :label: eq:laplacegen

   \begin{aligned}
   \nabla^2\phi = \frac{1}{J}\frac{\partial}{\partial
   u^i}\left(Jg^{ij}\frac{\partial\phi}{\partial u^j}\right)
   \end{aligned}

which can be expanded as

.. math::
   :label: eq:laplace_expand

   \begin{aligned}
   \nabla^2\phi = g^{ij}\frac{\partial^2\phi}{\partial u^i\partial u^j} +
   \underbrace{\frac{1}{J}\frac{\partial}{\partial
   u^i}\left(Jg^{ij}\right)}_{G^j}\frac{\partial\phi}{\partial u^j}
   \end{aligned}

where `G^j` must **not** be mistaken as the so called connection
coefficients (i.e. the Christoffel symbols of second kind). Setting
`\phi =
u^k` in equation :eq:`eq:laplacegen` gives
`\nabla^2u^k = G^k`. Expanding
:eq:`eq:laplacegen` and setting
`\left\{u^i\right\} = \left\{x, y, z\right\}` gives

.. math::
   :label: eq:general_laplacian

   \begin{aligned}
   \nabla^2f &= \nabla\cdot\nabla f = \nabla\cdot\left(\frac{\partial}{\partial
   x}\nabla x + \frac{\partial}{\partial y}\nabla y + \frac{\partial}{\partial
   z}\nabla z\right) \nonumber \\
   &= \frac{\partial^2 f}{\partial x^2}\left|\nabla x\right|^2 + \frac{\partial^2
   f}{\partial y^2}\left|\nabla y\right|^2 + \frac{\partial^2 f}{\partial z^2}\left|\nabla
   z\right|^2 \\ +2\frac{\partial^2 f}{\partial x\partial y}\left(\nabla x\cdot\nabla
   y\right) +2\frac{\partial^2 f}{\partial x\partial z}\left(\nabla x\cdot\nabla z\right)
   +2\frac{\partial^2 f}{\partial y\partial z}\left(\nabla y\cdot\nabla z\right)
   \nonumber \\ +\nabla^2x\frac{\partial f}{\partial x} +\nabla^2y\frac{\partial
   f}{\partial y} + \nabla^2z\frac{\partial f}{\partial z} \nonumber
   \end{aligned}

Curl defined as:

.. math::
   :label: eq:curlcurvilinear

   \begin{aligned}
   \nabla\times\mathbf{A} = \frac{1}{\sqrt{g}}\sum_k\left(\frac{\partial
   A_j}{\partial u_i} - \frac{\partial A_i}{\partial u_j}\right){\boldsymbol{e}}_k \qquad i,j,k
   \texttt{ cyc } 1,2,3 \end{aligned}

Cross-product relation between contravariant and covariant vectors:

.. math::

   \begin{aligned}
   {\boldsymbol{e}}^i = \frac{1}{J}\left({\boldsymbol{e}}_j \times {\boldsymbol{e}}_k\right) \qquad {\boldsymbol{e}}_i =
   J\left({\boldsymbol{e}}^j \times {\boldsymbol{e}}^k\right) \qquad i,j,k \texttt{ cyc } 1,2,3\end{aligned}

Derivation of operators in the BOUT++ Clebsch system
====================================================

The Clebsch system in BOUT++ goes like this

.. math::

   \begin{aligned}
       {\boldsymbol{B}}=&\nabla z \times \nabla x\\ =&{\boldsymbol{e}}^z \times {\boldsymbol{e}}^x\\
       J^{-1}{\boldsymbol{e}}_y=&{\boldsymbol{e}}^z \times {\boldsymbol{e}}^x\end{aligned}

We have

.. math::

   \begin{aligned}
       B{\overset{\text{def}}{=}}& \sqrt{{\boldsymbol{B}}\cdot{\boldsymbol{B}}} = \sqrt{J^{-1}{\boldsymbol{e}}_y\cdot
   J^{-1}{\boldsymbol{e}}_y} = \sqrt{J^{-2}g_{yy}} = J^{-1}\sqrt{g_{yy}}\end{aligned}

Further on

.. math::

   \begin{aligned}
       {\boldsymbol{B}}=&B{\boldsymbol{b}}\\ {\boldsymbol{b}}=&\frac{{\boldsymbol{B}}}{B}
       =\frac{J^{-1}{\boldsymbol{e}}_y}{J^{-1}\sqrt{g_{yy}}} =\frac{{\boldsymbol{e}}_y}{\sqrt{g_{yy}}}\end{aligned}

The parallel and perpendicular gradients
----------------------------------------

We have that

.. math::

   \begin{aligned}
       {\nabla}=& {\boldsymbol{e}}^i \partial_i = {\boldsymbol{e}}^x \partial_x + {\boldsymbol{e}}^y \partial_y +
       {\boldsymbol{e}}^z \partial_z\end{aligned}

and that

.. math::

   \begin{aligned}
       {\nabla}_\| =& \left({\boldsymbol{b}} \cdot {\nabla}\right) {\boldsymbol{b}} = {\boldsymbol{b}} {\boldsymbol{b}} \cdot {\nabla}=
       \frac{{\boldsymbol{e}}_y {\boldsymbol{e}}_y}{g_{yy}} \cdot {\nabla}= \frac{{\boldsymbol{e}}_y
       {\boldsymbol{e}}_y}{g_{yy}} \cdot {\boldsymbol{e}}^i \partial_i = \frac{{\boldsymbol{e}}_y}{g_{yy}}
       \partial_y\end{aligned}

so that

.. math::

   \begin{aligned}
       {\nabla}_\perp =& {\nabla}- {\nabla}_\|\\
   %
                   =& {\boldsymbol{e}}^x \partial_x + {\boldsymbol{e}}^y \partial_y + {\boldsymbol{e}}^z
       \partial_z - \frac{{\boldsymbol{e}}_y}{g_{yy}} \partial_y\\
   %
                   =& {\boldsymbol{e}}^x \partial_x + {\boldsymbol{e}}^y \partial_y + {\boldsymbol{e}}^z
       \partial_z - \frac{g_{yi}{\boldsymbol{e}}^i}{g_{yy}} \partial_y\\
   %
                   =& {\boldsymbol{e}}^x \partial_x + {\boldsymbol{e}}^y \partial_y + {\boldsymbol{e}}^z
       \partial_z - \frac{g_{yx}{\boldsymbol{e}}^x +g_{yy}{\boldsymbol{e}}^y +g_{yz}{\boldsymbol{e}}^z
       }{g_{yy}}\partial_y\\
   %
                   =& {\boldsymbol{e}}^x \left(\partial_x - \frac{g_{yx}}{g_{yy}}\partial_y\right)
       +  {\boldsymbol{e}}^z \left(\partial_z - \frac{g_{yz}}{g_{yy}}\partial_y\right)\end{aligned}

The perpendicular gradients in Laplacian inversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the Laplacian inversion BOUT++ currently neglects the parallel
`y` derivatives if `g_{xy}` and `g_{yz}` are non-zero,
thus

.. math::
   :label: eq:reduced_grad_perp

   \begin{aligned}
       {\nabla}_\perp \simeq& {\boldsymbol{e}}^x \partial_x +  {\boldsymbol{e}}^z \partial_z
       \end{aligned}

The Laplacian
-------------

We would here like to find an expression for the Laplacian

.. math::

   \begin{aligned}
       {\nabla}^2 = {\nabla\cdot}{\nabla}\end{aligned}

In general we have (using equation (2.6.39) in D’Haeseleer [haeseler]_)

.. math::
   :label: eq:divA

   \begin{aligned}
       {\nabla\cdot}{\boldsymbol{A}} = \frac{1}{J} \partial_i \left(JA^i\right)
       \end{aligned}

and that

.. math::

   \begin{aligned}
       A^i = {\boldsymbol{A}}\cdot {\boldsymbol{e}}^i\end{aligned}

In our case `A \to {\nabla}`, so that

.. math::

   \begin{aligned}
       {\nabla}^i = \left({\nabla}\right)\cdot {\boldsymbol{e}}^i = {\boldsymbol{e}}^i \cdot \left({\nabla}\right) = {\boldsymbol{e}}^i
       \cdot \left({\boldsymbol{e}}^j \partial_j\right) = g^{ij} \partial_j\end{aligned}

Thus

.. math::

   \begin{aligned}
       {\nabla}^2 =& \frac{1}{J} \partial_i \left(J g^{ij} \partial_j\right)\\ =&
       \frac{1}{J} g^{ij} J \partial_i \partial_j + \frac{1}{J} \partial_i \left(J
       g^{ij} \right) \partial_j\\ =& g^{ij} \partial_i \partial_j + G^j \partial_j\\\end{aligned}

where we have defined [#f1]_

.. math::

   \begin{aligned}
       G^j =& \frac{1}{J} \partial_i \left(J g^{ij} \right)\\ =& \frac{1}{J} \left(
       \partial_x \left[J g^{xj} \right] + \partial_y \left[J g^{yj} \right] + \partial_z \left[J
       g^{zj} \right] \right)\end{aligned}

By writing the terms out, we get

.. math::

   \begin{aligned}
       {\nabla}^2 =& g^{ij} \partial_i \partial_j + G^j \partial_j\\
   %
               =& \left(  g^{xj} \partial_x \partial_j + g^{yj} \partial_y \partial_j
       + g^{zj} \partial_z \partial_j\right) + \left(G^j \partial_j\right)\\
   %
               =& \quad \, \left(  g^{xx} \partial_x^2 + g^{yx} \partial_y \partial_x
       + g^{zx} \partial_z \partial_x\right) + \left(G^x \partial_x\right)\\ &+ \left(  g^{xy}
       \partial_x \partial_y + g^{yy} \partial_y^2 + g^{zy} \partial_z
       \partial_y\right) + \left(G^y \partial_y\right)\\ &+ \left(  g^{xz} \partial_x \partial_z
       + g^{yz} \partial_y \partial_z + g^{zz} \partial_z^y\right) + \left(G^z
       \partial_z\right)\end{aligned}

We now use that the metric tensor is symmetric (by definition), so that
`g^{ij}=g^{ji}`, and `g_{ij}=g_{ji}`, and that the partial
derivatives commutes for smooth functions
`\partial_i\partial_j=\partial_j\partial_i`. This gives

.. math::

   \begin{aligned}
       {\nabla}^2 =&\quad \, \left(g^{xx} \partial_x^2 \right) + \left(G^x \partial_x\right)\\ &+
       \left(g^{yy} \partial_y^2 \right) + \left(G^y \partial_y\right)\\ &+ \left(g^{zz}
       \partial_z^2\right) + \left(G^z \partial_z\right)\\ &+ 2\left( g^{xy} \partial_x
       \partial_y + g^{xz} \partial_x \partial_z + g^{yz} \partial_y \partial_z
       \right)\\
   %
              =&\quad \, \left(g^{xx} \partial_x^2\right) + \left( \frac{1}{J} \left[
   \partial_x \left\{J g^{xx} \right\} + \partial_y \left\{J g^{yx} \right\} + \partial_z \left\{J
   g^{zx} \right\} \right] \partial_x\right)\\ &+ \left(g^{yy} \partial_y^2\right) + \left( \frac{1}{J}
       \left[ \partial_x \left\{J g^{xy} \right\} + \partial_y \left\{J g^{yy} \right\} +
       \partial_z \left\{J g^{zy} \right\} \right] \partial_y\right)\\ &+ \left(g^{zz}
           \partial_z^2\right) + \left( \frac{1}{J} \left[ \partial_x \left\{J g^{xz} \right\} +
           \partial_y \left\{J g^{yz} \right\} + \partial_z \left\{J g^{zz} \right\} \right]
           \partial_z\right)\\ &+ 2\left( g^{xy} \partial_x \partial_y + g^{xz}
           \partial_x \partial_z + g^{yz} \partial_y \partial_z \right)\end{aligned}

Notice that `G^i` does not operate on `\partial_i`, but
rather that the two are multiplied together.

The parallel Laplacian
----------------------

We have that

.. math::

   \begin{aligned}
       {\nabla}_\| =& \left({\boldsymbol{b}} \cdot {\nabla}\right) {\boldsymbol{b}}\ = {\boldsymbol{b}} {\boldsymbol{b}} \cdot {\nabla}=
       \frac{{\boldsymbol{e}}_y {\boldsymbol{e}}_y}{g_{yy}} \cdot {\nabla}= \frac{{\boldsymbol{e}}_y
       {\boldsymbol{e}}_y}{g_{yy}} \cdot {\boldsymbol{e}}^i \partial_i = \frac{{\boldsymbol{e}}_y}{g_{yy}}
       \partial_y\end{aligned}

we have that

.. math::

   \begin{aligned}
       {\nabla}_\|^i =& \left(\frac{{\boldsymbol{e}}_y}{g_{yy}} \partial_y\right)\cdot {\boldsymbol{e}}^i =
       {\boldsymbol{e}}^i \cdot \left(\frac{{\boldsymbol{e}}_y}{g_{yy}} \partial_y\right)\end{aligned}

so that by equation :eq:`eq:divA`,

.. math::

   \begin{aligned}
       {\nabla}_\|^2 =& {\nabla\cdot}\left({\boldsymbol{b}} {\boldsymbol{b}} \cdot {\nabla}\right)\\ =&
       {\nabla\cdot}\left(\frac{{\boldsymbol{e}}_y}{g_{yy}} \cdot \partial_y\right)\\ =& \frac{1}{J}
       \partial_i \left( J{\boldsymbol{e}}^i \cdot \left[\frac{{\boldsymbol{e}}_y}{g_{yy}} \partial_y\right]
       \right)\\ =& \frac{1}{J} \partial_y \left(\frac{J}{g_{yy}} \partial_y\right)\end{aligned}

The perpendicular Laplacian
---------------------------

For the perpendicular Laplacian, we have that

.. math::

   \begin{aligned}
       {\nabla}_\perp^2 =& {\nabla}^2 - {\nabla}_\|^2\\ =& g^{ij} \partial_i \partial_j +
       G^j \partial_j -\frac{1}{J} \partial_y \left(\frac{J}{g_{yy}} \partial_y\right)\\
   %
               =& \quad \, \left(g^{xx} \partial_x^2\right) + \left( \frac{1}{J} \left[
   \partial_x \left\{J g^{xx} \right\} + \partial_y \left\{J g^{yx} \right\} + \partial_z \left\{J
   g^{zx} \right\} \right] \partial_x\right)\\ &+ \left(g^{yy} \partial_y^2\right) + \left( \frac{1}{J}
       \left[ \partial_x \left\{J g^{xy} \right\} + \partial_y \left\{J g^{yy} \right\} +
       \partial_z \left\{J g^{zy} \right\} \right] \partial_y\right)\\ &+ \left(g^{zz}
           \partial_z^2\right) + \left( \frac{1}{J} \left[ \partial_x \left\{J g^{xz} \right\} +
           \partial_y \left\{J g^{yz} \right\} + \partial_z \left\{J g^{zz} \right\} \right]
           \partial_z\right)\\ &+ 2\left( g^{xy} \partial_x \partial_y + g^{xz}
           \partial_x \partial_z + g^{yz} \partial_y \partial_z \right)\\ &-
           \frac{1}{J} \partial_y \left(\frac{J}{g_{yy}} \partial_y\right)\end{aligned}

The perpendicular Laplacian in Laplacian inversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Notice that BOUT++ currently assumes small parallel gradients in the
dependent variable in Laplacian inversion if `g_{xy}` and
`g_{yz}` are non-zero (if these are zero, the derivation can be
done directly from equation
:eq:`eq:reduced_grad_perp` instead), so that

.. math::

   \begin{aligned}
       {\nabla}_\perp^2 \simeq& \quad \, \left(g^{xx} \partial_x^2\right) + \left( \frac{1}{J}
       \left[ \partial_x \left\{J g^{xx} \right\} + \partial_y \left\{J g^{yx} \right\} +
       \partial_z \left\{J g^{zx} \right\} \right] \partial_x\right)\\ &+ \left(g^{zz}
           \partial_z^2\right) + \left( \frac{1}{J} \left[ \partial_x \left\{J g^{xz} \right\} +
           \partial_y \left\{J g^{yz} \right\} + \partial_z \left\{J g^{zz} \right\} \right]
           \partial_z\right)\\ &+ 2\left(g^{xz} \partial_x \partial_z\right)\\
   %
              =& \left(g^{xx} \partial_x^2\right) + G^x\partial_x + \left(g^{zz}
           \partial_z^2\right) + G^z \partial_z + 2\left(g^{xz} \partial_x \partial_z\right)\end{aligned}

The Poisson bracket operator
----------------------------

We will here derive the bracket operators, as they are used in BOUT++.

The electrostatic ExB velocity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Under electrostatic conditions, we have that
`{\boldsymbol{v}}_E =
-\frac{\nabla\phi\times{\boldsymbol{b}}}{B}`, which is
similar to
`{\boldsymbol{v}}={\boldsymbol{k}}\times\nabla\psi`
found in incompressible fluid flow

.. math::

   \begin{aligned}
       {\boldsymbol{v}}_E =& -\frac{\nabla\phi\times{\boldsymbol{b}}}{B}\\
                %
                =&-\frac{\nabla\phi\times{\boldsymbol{e}}_y}{
   \sqrt{g_{yy}}J^{-1}\sqrt{g_{yy}}}\\
                %
                =&-\frac{J}{g_{yy}}\nabla\phi\times{\boldsymbol{e}}_y\\
                %
                =&\frac{J}{g_{yy}}{\boldsymbol{e}}_y\times\nabla\phi\\
                %
                =&\frac{J}{g_{yy}}{\boldsymbol{e}}_y\times \left({\boldsymbol{e}}^x\partial_x + {\boldsymbol{e}}^y\partial_y +
   {\boldsymbol{e}}^z\partial_z\right)\phi\\
                %
                =&\frac{J}{g_{yy}} \left(g_{yx}{\boldsymbol{e}}^x + g_{yy}{\boldsymbol{e}}^y +
   g_{yz}{\boldsymbol{e}}^z\right) \times \left({\boldsymbol{e}}^x\partial_x + {\boldsymbol{e}}^y\partial_y +
   {\boldsymbol{e}}^z\partial_z\right)\phi\\
                %
                =&\frac{J}{g_{yy}} \left( g_{yx}{\boldsymbol{e}}^x\times{\boldsymbol{e}}^x\partial_x +
   g_{yy}{\boldsymbol{e}}^y\times{\boldsymbol{e}}^x\partial_x + g_{yz}{\boldsymbol{e}}^z\times{\boldsymbol{e}}^x\partial_x
   \right.  \\ &\quad\; + g_{yx}{\boldsymbol{e}}^x\times{\boldsymbol{e}}^y\partial_y +
   g_{yy}{\boldsymbol{e}}^y\times{\boldsymbol{e}}^y\partial_y + g_{yz}{\boldsymbol{e}}^z\times{\boldsymbol{e}}^y\partial_y
   \\ &\quad\; \left.  + g_{yx}{\boldsymbol{e}}^x\times{\boldsymbol{e}}^z\partial_z +
   g_{yy}{\boldsymbol{e}}^y\times{\boldsymbol{e}}^z\partial_z + g_{yz}{\boldsymbol{e}}^z\times{\boldsymbol{e}}^z\partial_z
   \right) \phi\\
                %
                =&\frac{J}{g_{yy}} \left( - g_{yy}{\boldsymbol{e}}^y\times{\boldsymbol{e}}^x\partial_x +
   g_{yz}{\boldsymbol{e}}^z\times{\boldsymbol{e}}^x\partial_x \right.  \\ &\quad +
   g_{yx}{\boldsymbol{e}}^x\times{\boldsymbol{e}}^y\partial_y - g_{yz}{\boldsymbol{e}}^z\times{\boldsymbol{e}}^y\partial_y
   \\ &\quad \left.  - g_{yx}{\boldsymbol{e}}^x\times{\boldsymbol{e}}^z\partial_z +
   g_{yy}{\boldsymbol{e}}^y\times{\boldsymbol{e}}^z\partial_z \right) \phi\\
                %
                =&\frac{1}{g_{yy}} \left( - g_{yy}{\boldsymbol{e}}_z\partial_x +
   g_{yz}{\boldsymbol{e}}_y\partial_x + g_{yx}{\boldsymbol{e}}_z\partial_y - g_{yz}{\boldsymbol{e}}_x\partial_y
   - g_{yx}{\boldsymbol{e}}_y\partial_z + g_{yy}{\boldsymbol{e}}_x\partial_z \right) \phi\end{aligned}

The electrostatic ExB advection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The electrostatic `E\times B` advection operator thus becomes

.. math::

   \begin{aligned}
       {\boldsymbol{v}}_E\cdot\nabla =& -\frac{\nabla\phi\times{\boldsymbol{b}}}{B}\cdot\nabla\\
       %
       =&\frac{1}{g_{yy}} \left( - g_{yy}{\boldsymbol{e}}_z\partial_x +
       g_{yz}{\boldsymbol{e}}_y\partial_x + g_{yx}{\boldsymbol{e}}_z\partial_y -
       g_{yz}{\boldsymbol{e}}_x\partial_y - g_{yx}{\boldsymbol{e}}_y\partial_z +
       g_{yy}{\boldsymbol{e}}_x\partial_z \right) \phi \cdot\left({\boldsymbol{e}}^x\partial_x +
       {\boldsymbol{e}}^y\partial_y + {\boldsymbol{e}}^z\partial_z\right)\\
       %
       =& \frac{1}{g_{yy}} \left( - g_{yy}\partial_x\phi\partial_z +
       g_{yz}\partial_x\phi\partial_y + g_{yx}\partial_y\phi\partial_z -
       g_{yz}\partial_y\phi\partial_x - g_{yx}\partial_z\phi\partial_y +
       g_{yy}\partial_z\phi\partial_x \right)\\
       %
       =& \frac{1}{g_{yy}} \left( \left[ g_{yy}\partial_z\phi - g_{yz}\partial_y\phi
   \right]\partial_x + \left[ g_{yz}\partial_x\phi - g_{yx}\partial_z\phi \right]\partial_y +
   \left[ g_{yx}\partial_y\phi - g_{yy}\partial_x\phi \right]\partial_z \right)\\
       %
       =& \frac{1}{g_{yy}} \left( g_{yx}\{\phi, \cdot\}_{y,z} + g_{yy}\{\phi,
       \cdot\}_{z,x} + g_{yz}\{\phi, \cdot\}_{x,y} \right)\end{aligned}

Where we have used the definition of the Poisson bracket

.. math::

   \begin{aligned}
       \{a, b\}_{i,j} = \left(\partial_i a\right) \partial_j b - \left(\partial_j a\right)
       \partial_i b\end{aligned}

The pure solenoidal advection is thus

.. math::
   :label: eq:brackets

       B{\boldsymbol{v}}_E\cdot\nabla =& -\nabla\phi\times{\boldsymbol{b}}\cdot\nabla\\
       =& {\boldsymbol{b}} \times \nabla\phi\cdot\nabla\\
       =& \frac{\sqrt{g_{yy}}}{Jg_{yy}} \left( g_{yx}\{\phi, \cdot\}_{y,z} +
       g_{yy}\{\phi, \cdot\}_{z,x} + g_{yz}\{\phi, \cdot\}_{x,y} \right) \\
       =& \frac{1}{J\sqrt{g_{yy}}} \left( g_{yx}\{\phi, \cdot\}_{y,z} + g_{yy}\{\phi,
       \cdot\}_{z,x} + g_{yz}\{\phi, \cdot\}_{x,y} \right)


The brackets operator in BOUT++
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Notice that the (phi,f)@ operators in BOUT++ returns
`-\frac{\nabla\phi\times{\boldsymbol{b}}}{B}\cdot\nabla f`
rather than
`-\nabla\phi\times{\boldsymbol{b}}\cdot\nabla f`.

Notice also that the Arakawa brackets neglects the `\partial_y`
derivative terms (the `y`-derivative terms) if `g_{xy}` and
`g_{yz}` are non-zero, so for the Arakawa brackets, BOUT++ returns

.. math::

   \begin{aligned}
       {\boldsymbol{v}}_E\cdot\nabla =& -\frac{\nabla\phi\times{\boldsymbol{b}}}{B}\cdot\nabla\\
       %
       \simeq& \frac{1}{g_{yy}} \left( g_{yy}\{\phi, \cdot\}_{z,x} \right)\\
       %
       =& \partial_z\phi\partial_x - \partial_x\phi\partial_z\end{aligned}

Divergence of ExB velocity
==========================

.. math::

   \begin{aligned}
   {\boldsymbol{v}}_{ExB} = \frac{{\boldsymbol{b}}\times\nabla\phi}{B}\end{aligned}

Using

.. math::

   \begin{aligned}
   \nabla\cdot\left({\boldsymbol{F}}\times{\boldsymbol{G}}\right) = \left(\nabla\times{\boldsymbol{F}}\right)\cdot{\boldsymbol{G}} -
   {\boldsymbol{F}}\cdot\left(\nabla\times{\boldsymbol{G}}\right)\end{aligned}

the divergence of the
`{\boldsymbol{E}}\times{\boldsymbol{B}}`
velocity can be written as

.. math::
   :label: eq:exb1

   \begin{aligned}
   \nabla\cdot\left(\frac{1}{B}{\boldsymbol{b}}\times\nabla\phi\right) =
   \left[\nabla\times\left(\frac{1}{B}{\boldsymbol{b}}\right)\right]\cdot\nabla\phi -
   \frac{1}{B}{\boldsymbol{b}}\cdot\nabla\times\nabla\phi
   \end{aligned}

The second term on the right is identically zero (curl of a nablaient).
The first term on the right can be expanded as

.. math::

   \begin{aligned}
   \left[\nabla\times\left(\frac{1}{B}{\boldsymbol{b}}\right)\right]\cdot\nabla\phi =
   \left[\nabla\left(\frac{1}{B}\right)\times{\boldsymbol{b}} +
   \frac{1}{B}\nabla\times{\boldsymbol{b}}\right]\cdot\nabla\phi\end{aligned}

Using

.. math::

   \begin{aligned}
   {\boldsymbol{b}}\times{\boldsymbol{\kappa}} = \nabla\times{\boldsymbol{b}} -
   {\boldsymbol{b}}\left[{\boldsymbol{b}}\cdot\left(\nabla\times{\boldsymbol{b}}\right)\right]\end{aligned}

this becomes:

.. math::

   \begin{aligned}
     \nabla\cdot\left(\frac{1}{B}{\boldsymbol{b}}\times\nabla\phi\right) =
     &-{\boldsymbol{b}}\times\nabla\left(\frac{1}{B}\right)\cdot\nabla\phi \\ &+
     \frac{1}{B}{\boldsymbol{b}}\times{\boldsymbol{\kappa}}\cdot\nabla\phi \\ &+
     \left[{\boldsymbol{b}}\cdot\left(\nabla\times{\boldsymbol{b}}\right)\right]{\boldsymbol{b}}\cdot\nabla\phi\end{aligned}

Alternatively, equation :eq:`eq:exb1` can be expanded as

.. math::

   \begin{aligned}
     \nabla\cdot\left(\frac{1}{B}{\boldsymbol{b}}\times\nabla\phi\right) =&
       -B{\boldsymbol{b}}\times\nabla\left(\frac{1}{B^2}\right)\cdot\nabla\phi +
       \frac{1}{B^2}\nabla\times{\boldsymbol{B}}\cdot\nabla\phi \\ =&
       -B{\boldsymbol{b}}\times\nabla\left(\frac{1}{B^2}\right)\cdot\nabla\phi +
       \frac{1}{B^2}{\boldsymbol{J}}\cdot\nabla\phi\end{aligned}

.. math::

   \begin{aligned}
   \nabla\cdot\left(n\frac{\mathbf{b}\times\nabla\phi}{B}\right) &=& \frac{1}{J}\frac{\partial}{\partial\psi}\left(Jn\frac{\partial\phi}{\partial z} \right) - \frac{1}{J}\frac{\partial}{\partial z}\left(Jn\frac{\partial\phi}{\partial\psi}\right)  \\
                                                                 &+& \frac{1}{J}\frac{\partial}{\partial\psi}\left(Jn\frac{g^{\psi\psi}g^{yz}}{B^2}\frac{\partial\phi}{\partial y}\right) - \frac{1}{J}\frac{\partial}{\partial y}\left(Jn\frac{g^{\psi\psi}g^{yz}}{B^2}\frac{\partial\phi}{\partial\psi}\right)\end{aligned}

.. [haeseler] Haeseler, W. D.: Flux Coordinates and Magnetic Field Structure, Springer-Verlag, 1991, ISBN 3-540-52419-3

.. rubric:: Footnotes

.. [#f1] Notice that `G^i` is **not** the same as the
     *Christoffel symbols of second kind* (also known as the
     *connection coefficients* or
     `\Gamma^i_{jk}={\boldsymbol{e}}^i\cdot\partial_k
     {\boldsymbol{e}}_j`), although the derivation of the two are
     quite similar.  | We find that
     `\Gamma^i_{ji}={\boldsymbol{e}}^i\cdot\partial_i
     {\boldsymbol{e}}_j = {\nabla\cdot}{\boldsymbol{e}}_j`, whereas
     using equation :eq:`eq:divA` leads to
     `G^i={\boldsymbol{e}}^i\cdot\partial_i {\boldsymbol{e}}^j =
     {\nabla\cdot} {\boldsymbol{e}}^j`, since `g^{ji}=g^{ij}`
     due to symmetry.
