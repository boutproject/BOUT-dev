=========================
Field-aligned coordinates
=========================

:Author: B.Dudson§, M.V.Umansky, L.C.Wang, X.Q.Xu, L.L.LoDestro
§Department of Physics, University of York, UK
Lawrence Livermore National Laboratory, USA
IFTS, China

.. role:: raw-latex(raw)
   :format: latex
..

.. raw:: latex

   \maketitle

.. raw:: latex

   \tableofcontents

Introduction
============

This manual covers the field-aligned coordinate system used in many
BOUT++ tokamak models, and useful derivations and expressions.

.. _sec:coordinates:

Orthogonal toroidal coordinates
===============================

Starting with an orthogonal toroidal coordinate system
:math:`\left(\psi, \theta,
\zeta\right)`, where :math:`\psi` is the poloidal flux, :math:`\theta`
the poloidal angle (from :math:`0` to :math:`2\pi`), and :math:`\zeta`
the toroidal angle (also :math:`0` to :math:`2\pi`). We have that the
magnetic field :math:`\ensuremath{\boldsymbol{B}}` can be expressed as

.. math::

   \begin{aligned}
    \ensuremath{\boldsymbol{B}}=& B_\theta \nabla \theta + B_\zeta \nabla \zeta\\ =& B_\theta
       \ensuremath{\boldsymbol{e}}_\theta + B_\zeta \ensuremath{\boldsymbol{e}}_\zeta\\ =& \ensuremath{B_{\text{pol}}}h_\theta \ensuremath{\boldsymbol{e}}_\theta + \ensuremath{B_{\text{tor}}}
       R \ensuremath{\boldsymbol{e}}_\zeta\\ =& \ensuremath{B_{\text{pol}}}\hat{\ensuremath{\boldsymbol{e}}}_\theta + \ensuremath{B_{\text{tor}}}\hat{\ensuremath{\boldsymbol{e}}}_\zeta\end{aligned}

 The magnitudes of the unit vectors are

.. math::

   \begin{aligned}
   \left|\hat{\ensuremath{\boldsymbol{e}}}_\psi\right| = \frac{1}{R\left|\ensuremath{B_{\text{pol}}}\right|} \qquad \left|\hat{\ensuremath{\boldsymbol{e}}}_\theta\right| = \ensuremath{h_\theta}
   \qquad \left|\hat{\ensuremath{\boldsymbol{e}}}_\zeta\right| = R\end{aligned}

 where :math:`\ensuremath{h_\theta}` is the poloidal arc length per
radian. The coordinate system is right handed, so
:math:`\hat{\ensuremath{\boldsymbol{e}}}_\psi\times\hat{\ensuremath{\boldsymbol{e}}}_\theta = \hat{\ensuremath{\boldsymbol{e}}}_\zeta`,
:math:`\hat{\ensuremath{\boldsymbol{e}}}_\psi\times\hat{\ensuremath{\boldsymbol{e}}}_\zeta = -\hat{\ensuremath{\boldsymbol{e}}}_\theta`
and
:math:`\hat{\ensuremath{\boldsymbol{e}}}_\theta\times\hat{\ensuremath{\boldsymbol{e}}}_\zeta = \hat{\ensuremath{\boldsymbol{e}}}_\psi`.
The covariant metric coefficients are

.. math::

   \begin{aligned}
   g_{\psi\psi} = \frac{1}{\left(R\left|\ensuremath{B_{\text{pol}}}\right|\right)^2} \qquad g_{\theta\theta} =
   h_\theta^2 \qquad g_{\zeta\zeta} = R^2\end{aligned}

 and the magnitudes of the reciprocal vectors are therefore

.. math::

   \begin{aligned}
   \left|\nabla\psi\right| = R\left|\ensuremath{B_{\text{pol}}}\right| \qquad \left|\nabla\theta\right| = \frac{1}{h_\theta}
   \qquad \left|\nabla\zeta\right| = \frac{1}{R}\end{aligned}

 Because the coordinate system is orthogonal, :math:`g^{ii} = 1/g_{ii}`
and so the cross-products can be calculated as

.. math::

   \begin{aligned}
   \nabla\psi\times\nabla\theta = &\hat{\ensuremath{\boldsymbol{e}}}^\psi\times \hat{\ensuremath{\boldsymbol{e}}}^\theta =
       g^{\psi\psi}\ensuremath{\boldsymbol{e}}_\psi\times g^{\theta\theta}\ensuremath{\boldsymbol{e}}_\theta \nonumber \\ =
       & g^{\psi\psi}g^{\theta\theta}h_\psi h_\theta
       \hat{e}_\psi\times\hat{e}_\theta \nonumber \\ = &\frac{1}{h_\psi
   h_\theta}\hat{\ensuremath{\boldsymbol{e}}}_\zeta = \frac{R\left|\ensuremath{B_{\text{pol}}}\right|}{h_\theta}\hat{e}_\zeta\end{aligned}

 Similarly,

.. math::

   \begin{aligned}
   \nabla\psi\times\nabla\zeta = -\left|\ensuremath{B_{\text{pol}}}\right|\hat{\ensuremath{\boldsymbol{e}}}_\theta \qquad
   \nabla\theta\times\nabla\zeta = \frac{1}{Rh_\theta}\hat{\ensuremath{\boldsymbol{e}}}_\psi =
   \frac{1}{h_\theta R^2\left|\ensuremath{B_{\text{pol}}}\right|}\nabla \psi\end{aligned}

Field-aligned coordinates
=========================

In order to efficiently simulate (predominantly) field-aligned
structures, grid-points are placed in a field-aligned coordinate system.
We define
:math:`\sigma_{B\theta} \equiv \ensuremath{B_{\text{pol}}}/ \left|\ensuremath{B_{\text{pol}}}\right|`
i.e. the sign of the poloidal field. The new coordinates
:math:`\left(x,y,z\right)` are defined by:

.. math::

   \begin{aligned}
   x = \ensuremath{\sigma_{B\theta}}\left(\psi - \psi_0\right) \qquad y = \theta \qquad z = \sigma_{B\theta}
   \left(\zeta - \int_{\theta_0}^{\theta}\nu\left(\psi,\theta\right)d\theta\right)
   \label{eq:coordtransform}\end{aligned}

 Where :math:`\nu` is the local field-line pitch given by

.. math::

   \begin{aligned}
   \nu\left(\psi, \theta\right) = \frac{\ensuremath{\boldsymbol{B}}\cdot\nabla\zeta}{\ensuremath{\boldsymbol{B}}\cdot\nabla\theta} =
   \frac{\ensuremath{B_{\text{tor}}}\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}R} = \frac{\left(F/R\right)h_\theta}{\ensuremath{B_{\text{pol}}}R} = FJ/R^2\end{aligned}

 where :math:`F=\ensuremath{B_{\text{tor}}}R` is a function only of
:math:`\psi` (sometimes called the poloidal current function).

The coordinate system is chosen so that :math:`x` increases radially
outwards, from plasma to the wall. The sign of the toroidal field
:math:`\ensuremath{B_{\text{tor}}}` can then be either + or -.

The contravariant basis vectors are therefore

.. math::

   \begin{aligned}
   \nabla x = \ensuremath{\sigma_{B\theta}}\nabla \psi \qquad \nabla y = \nabla \theta \qquad \nabla z =
   \ensuremath{\sigma_{B\theta}}\left(\nabla\zeta - \left[\int_{\theta_0}^\theta\ensuremath{\frac{\partial \nu\left(\psi,
   \theta\right)}{\partial \psi}} d\theta\right] \nabla\psi - \nu\left(\psi, \theta\right)\nabla\theta\right)\end{aligned}

 The term in square brackets is the integrated local shear:

.. math::

   \begin{aligned}
   I = \int_{y_0}^y\frac{\partial\nu\left(x, y\right)}{\partial\psi}dy\end{aligned}

Magnetic field
--------------

Magnetic field is given in Clebsch form by:

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{B}}= \nabla z\times \nabla x = \frac{1}{J}\ensuremath{\boldsymbol{e}}_y\end{aligned}

 The contravariant components of this are then

.. math::

   \begin{aligned}
   B^y = \frac{\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}} \qquad B^x = B^z = 0\end{aligned}

 i.e. :math:`\ensuremath{\boldsymbol{B}}` can be written as

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{B}}= \frac{\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}}\ensuremath{\boldsymbol{e}}_y\end{aligned}

 and the covariant components calculated using :math:`g_{ij}` as

.. math::

   \begin{aligned}
   B_x = \ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}I R \qquad B_y = \frac{B^2 \ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}} \qquad B_z = \ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}R\end{aligned}

 The unit vector in the direction of equilibrium
:math:`\ensuremath{\boldsymbol{B}}` is therefore

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{b}} = \frac{1}{JB}\ensuremath{\boldsymbol{e}}_y = \frac{1}{JB}\left[g_{xy}\nabla x + g_{yy}\nabla y
   + g_{yz}\nabla z\right]\end{aligned}

Jacobian and metric tensors
---------------------------

The Jacobian of this coordinate system is

.. math::

   \begin{aligned}
   J^{-1} \equiv \left(\nabla x\times\nabla y\right)\cdot\nabla z = \ensuremath{B_{\text{pol}}}/ \ensuremath{h_\theta}\end{aligned}

 which can be either positive or negative, depending on the sign of
:math:`\ensuremath{B_{\text{pol}}}`. The contravariant metric tensor is
given by:

.. math::

   \begin{aligned}
   g^{ij} \equiv \ensuremath{\boldsymbol{e}}^i \cdot\ensuremath{\boldsymbol{e}}^j \equiv \nabla u^i \cdot \nabla u^j = \left(%
   \begin{array}{ccc}
   \left(R\ensuremath{B_{\text{pol}}}\right)^2 & 0 & -I\left(R\ensuremath{B_{\text{pol}}}\right)^2 \\
   0 & 1 / \ensuremath{h_\theta}^2 & -\ensuremath{\sigma_{B\theta}}\nu / \ensuremath{h_\theta}^2 \\
   -I\left(R\ensuremath{B_{\text{pol}}}\right)^2 & -\ensuremath{\sigma_{B\theta}}\nu / \ensuremath{h_\theta}^2 & I^2\left(R\ensuremath{B_{\text{pol}}}\right)^2 + B^2 /
   \left(R\ensuremath{B_{\text{pol}}}\right)^2
   \end{array}
   %
    \right)\end{aligned}

 and the covariant metric tensor:

.. math::

   \begin{aligned}
   g_{ij} \equiv \ensuremath{\boldsymbol{e}}_i \cdot\ensuremath{\boldsymbol{e}}_j = \left(%
   \begin{array}{ccc}
   I^2 R^2 + 1 / \ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}& \ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}\ensuremath{h_\theta}I R / \ensuremath{B_{\text{pol}}}& I R^2 \\
   \ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}\ensuremath{h_\theta}I R / \ensuremath{B_{\text{pol}}}& B^2\ensuremath{h_\theta}^2 / \ensuremath{B_{\text{pol}}}^2 & \ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}\ensuremath{h_\theta}R / \ensuremath{B_{\text{pol}}}\\
   I R^2 & \ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}\ensuremath{h_\theta}R / \ensuremath{B_{\text{pol}}}& R^2
   \end{array}
   %
    \right)\end{aligned}

Differential operators
----------------------

The derivative of a scalar field :math:`f` along the *unperturbed*
magnetic field :math:`\ensuremath{\boldsymbol{b}}_0` is given by

.. math::

   \begin{aligned}
   \partial^0_{||}f \equiv \ensuremath{\boldsymbol{b}}_0 \cdot\nabla f =
   \frac{1}{\sqrt{g_{yy}}}\ensuremath{\frac{\partial f}{\partial y}} = \frac{\ensuremath{B_{\text{pol}}}}{B\ensuremath{h_\theta}}\ensuremath{\frac{\partial f}{\partial y}}\end{aligned}

 whilst the parallel divergence is given by

.. math::

   \begin{aligned}
   \nabla^0_{||}f = B_0\partial^0_{||}\left(\frac{f}{B_0}\right)\end{aligned}

 Using equation (`[eq:general_laplacian] <#eq:general_laplacian>`__),
the Laplacian operator is given by

.. math::

   \begin{aligned}
   \nabla^2 = &\frac{\partial^2}{\partial x^2}\left|\nabla x\right|^2 +
       \frac{\partial^2}{\partial y^2}\left|\nabla y\right|^2 +
       \frac{\partial^2}{\partial z^2}\left|\nabla z\right|^2 \nonumber \\
       &-2\frac{\partial^2}{\partial x\partial z}I\left(R\ensuremath{B_{\text{pol}}}\right)^2 -
       2\frac{\partial^2}{\partial y\partial z}\frac{\nu}{h_\theta^2}\\
       &+\frac{\partial}{\partial x}\nabla^2x + \frac{\partial}{\partial
   y}\nabla^2y + \frac{\partial}{\partial z}\nabla^2z \nonumber\end{aligned}

 Using equation (`[eq:laplace_expand] <#eq:laplace_expand>`__) for
:math:`\nabla^2x = G^x` etc, the values are

.. math::

   \begin{aligned}
   \nabla^2x = \frac{\ensuremath{B_{\text{pol}}}}{h_\theta}\frac{\partial}{\partial x}\left(h_\theta
   R^2\ensuremath{B_{\text{pol}}}\right) \qquad \nabla^2y = \frac{\ensuremath{B_{\text{pol}}}}{h_\theta}\frac{\partial}{\partial
   y}\left(\frac{1}{\ensuremath{B_{\text{pol}}}h_\theta}\right)\end{aligned}

.. math::

   \begin{aligned}
   \nabla^2z = -\frac{\ensuremath{B_{\text{pol}}}}{h_\theta}\left[\frac{\partial}{\partial x}\left(IR^2\ensuremath{B_{\text{pol}}}
   h_\theta\right) + \frac{\partial}{\partial y}\left(\frac{\nu}{\ensuremath{B_{\text{pol}}}h_\theta}\right)\right]\end{aligned}

 Neglecting some parallel derivative terms, the perpendicular Laplacian
can be written:

.. math::

   \begin{aligned}
   \nabla_\perp^2= \ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\left[\ensuremath{\frac{\partial^2 }{\partial {x}^2}} - 2I\frac{\partial^2}{\partial z\partial x} +
   \left(I^2 + \frac{B^2}{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^4}\right)\ensuremath{\frac{\partial^2 }{\partial {z}^2}}\right] + \nabla^2 x \ensuremath{\frac{\partial }{\partial x}} +
   \nabla^2 z\ensuremath{\frac{\partial }{\partial z}}\end{aligned}

 The second derivative along the equilibrium field

.. math::

   \begin{aligned}
   \partial^2_{||}\phi = \partial^0_{||}\left(\partial^0_{||}\phi\right) =
   \frac{1}{\sqrt{g_{yy}}}\ensuremath{\frac{\partial }{\partial y}}\left(\frac{1}{\sqrt{g_{yy}}}\right)\ensuremath{\frac{\partial  \phi}{\partial y}}
   + \frac{1}{g_{yy}}\frac{\partial^2\phi}{\partial y^2}\end{aligned}

 A common expression (the Poisson bracket in reduced MHD) is (from
equation (`[eq:brackets] <#eq:brackets>`__)):

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{b}}_0\cdot\nabla\phi\times\nabla A =
   \frac{1}{J\sqrt{g_{yy}}}\left[\left(g_{yy}\ensuremath{\frac{\partial \phi}{\partial z}} -
   g_{yz}\ensuremath{\frac{\partial \phi}{\partial y}}\right)\ensuremath{\frac{\partial A}{\partial x}} + \left(g_{yz}\ensuremath{\frac{\partial \phi}{\partial x}} -
   g_{xy}\ensuremath{\frac{\partial \phi}{\partial z}}\right)\ensuremath{\frac{\partial A}{\partial y}} + \left(g_{xy}\ensuremath{\frac{\partial \phi}{\partial y}} -
   g_{yy}\ensuremath{\frac{\partial \phi}{\partial x}}\right)\ensuremath{\frac{\partial A}{\partial z}}\right]\end{aligned}

 The perpendicular nabla operator:

.. math::

   \begin{aligned}
   \nabla_\perp \equiv& \nabla - \ensuremath{\boldsymbol{b}}\left(\ensuremath{\boldsymbol{b}}\cdot\nabla\right) \\ =& \nabla
       x\left(\ensuremath{\frac{\partial }{\partial x}} - \frac{g_{xy}}{\left(JB\right)^2}\ensuremath{\frac{\partial }{\partial y}}\right) + \nabla
       z\left(\ensuremath{\frac{\partial }{\partial z}} - \frac{g_{yz}}{\left(JB\right)^2}\ensuremath{\frac{\partial }{\partial y}}\right)\end{aligned}

.. _sec:jxb_fac:

J x B in field-aligned coordinates
----------------------------------

Components of the magnetic field in field-aligned coordinates:

.. math::

   \begin{aligned}
   B^y = \frac{\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}} \qquad B^x = B^z = 0\end{aligned}

 and

.. math::

   \begin{aligned}
   B_x = \ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}I R \qquad B_y = \frac{B^2\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}} \qquad B_z = \ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}R\end{aligned}

 Calculate current
:math:`\ensuremath{\boldsymbol{J}}= \frac{1}{\mu}\ensuremath{\nabla\times \ensuremath{\boldsymbol{B}} }`

.. math::

   \begin{aligned}
   \left(\ensuremath{\nabla\times \ensuremath{\boldsymbol{B}} }\right)^x = \frac{1}{J}\left(\ensuremath{\frac{\partial B_z}{\partial y}} - \ensuremath{\frac{\partial B_y}{\partial z}}\right) = 0\end{aligned}

 since :math:`\ensuremath{B_{\text{tor}}}R` is a flux-surface quantity,
and :math:`\ensuremath{\boldsymbol{B}}` is axisymmetric.

.. math::

   \begin{aligned}
   \left(\ensuremath{\nabla\times \ensuremath{\boldsymbol{B}} }\right)^y =& -\ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial x}}\left(\ensuremath{B_{\text{tor}}}R\right) \\
       \left(\ensuremath{\nabla\times \ensuremath{\boldsymbol{B}} }\right)^z =&
       \frac{\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}}\left[\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B^2\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) -
       \ensuremath{\sigma_{B\theta}}\ensuremath{\frac{\partial }{\partial y}}\left(\ensuremath{B_{\text{tor}}}I R\right)\right]\end{aligned}

 The second term can be simplified, again using
:math:`\ensuremath{B_{\text{tor}}}R` constant on flux-surfaces:

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial }{\partial y}}\left(\ensuremath{B_{\text{tor}}}I R\right) = \ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}R\ensuremath{\frac{\partial \nu}{\partial x}} \qquad \nu =
   \frac{\ensuremath{h_\theta}\ensuremath{B_{\text{tor}}}}{R\ensuremath{B_{\text{pol}}}}\end{aligned}

 From these, calculate covariant components:

.. math::

   \begin{aligned}
   \left(\ensuremath{\nabla\times \ensuremath{\boldsymbol{B}} }\right)_x =& -\ensuremath{B_{\text{tor}}}I R \ensuremath{\frac{\partial }{\partial x}}\left(\ensuremath{B_{\text{tor}}}R\right) +
       \frac{IR^2\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}}\left[\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B^2\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) - \ensuremath{B_{\text{tor}}}
       R\ensuremath{\frac{\partial \nu}{\partial x}}\right] \nonumber\\
   %
   \left(\ensuremath{\nabla\times \ensuremath{\boldsymbol{B}} }\right)_y =& -\ensuremath{\sigma_{B\theta}}\frac{B^2\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\ensuremath{\frac{\partial }{\partial x}}\left(\ensuremath{B_{\text{tor}}}R\right) +
       \ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}R\left[\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B^2\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) - \ensuremath{B_{\text{tor}}}R\ensuremath{\frac{\partial \nu}{\partial x}}\right]
       \label{eq:curlb_y}\\
   %
   \left(\ensuremath{\nabla\times \ensuremath{\boldsymbol{B}} }\right)_z =& -\ensuremath{B_{\text{tor}}}R\ensuremath{\frac{\partial }{\partial x}}\left(\ensuremath{B_{\text{tor}}}R\right) +
       \frac{R^2\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}}\left[\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B^2\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) - \ensuremath{B_{\text{tor}}}
       R\ensuremath{\frac{\partial \nu}{\partial x}}\right] \nonumber\end{aligned}

 Calculate
:math:`\ensuremath{\boldsymbol{J}}\times\ensuremath{\boldsymbol{B}}`
using

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{e}}^i = \frac{1}{J}\left(\ensuremath{\boldsymbol{e}}_j \times \ensuremath{\boldsymbol{e}}_k\right) \qquad \ensuremath{\boldsymbol{e}}_i =
   J\left(\ensuremath{\boldsymbol{e}}^j \times \ensuremath{\boldsymbol{e}}^k\right) \qquad i,j,k \texttt{ cyc } 1,2,3\end{aligned}

 gives

.. math::

   \begin{aligned}
   \mu_0 \left(\ensuremath{\boldsymbol{J}}\times\ensuremath{\boldsymbol{B}}\right)^x =& \frac{1}{J}\left[\left(\ensuremath{\nabla\times \ensuremath{\boldsymbol{B}} }\right)_y B_z -
   \left(\ensuremath{\nabla\times \ensuremath{\boldsymbol{B}} }\right)_z B_y \right]\\ =& -\frac{\ensuremath{B_{\text{pol}}}^3
   R^2}{\ensuremath{h_\theta}}\left[\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B^2\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) - \ensuremath{B_{\text{tor}}}R\ensuremath{\frac{\partial \nu}{\partial x}}\right]\end{aligned}

 Covariant components of :math:`\nabla P`:

.. math::

   \begin{aligned}
   \left(\nabla P\right)_x = \ensuremath{\frac{\partial P}{\partial x}} \qquad \left(\nabla P\right)_y = \left(\nabla P\right)_z = 0\end{aligned}

 and contravariant:

.. math::

   \begin{aligned}
   \left(\nabla P\right)^x = \ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\ensuremath{\frac{\partial P}{\partial x}} \qquad \left(\nabla P\right)^y = 0 \qquad
   \left(\nabla P\right)^z = -I\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\ensuremath{\frac{\partial P}{\partial x}}\end{aligned}

 Hence equating contravariant x components of
:math:`\ensuremath{\boldsymbol{J}}\times\ensuremath{\boldsymbol{B}}= \nabla P`,

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial }{\partial x}}\left(\frac{B^2\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) - \ensuremath{B_{\text{tor}}}
   R\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{B_{\text{tor}}}\ensuremath{h_\theta}}{R\ensuremath{B_{\text{pol}}}}\right) + \frac{\mu_0\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\ensuremath{\frac{\partial P}{\partial x}} =
   0
   \label{eq:xbalance}\end{aligned}

 Use this to calculate :math:`\ensuremath{h_\theta}` profiles (need to
fix :math:`\ensuremath{h_\theta}` at one radial location).

Close to x-points, the above expression becomes singular, so a better
way to write it is:

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial }{\partial x}}\left(B^2\ensuremath{h_\theta}\right) - \ensuremath{h_\theta}\ensuremath{B_{\text{pol}}}\ensuremath{\frac{\partial \ensuremath{B_{\text{pol}}}}{\partial x}} - \ensuremath{B_{\text{tor}}}
   R\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{B_{\text{tor}}}\ensuremath{h_\theta}}{R}\right) + \mu_0\ensuremath{h_\theta}\ensuremath{\frac{\partial P}{\partial x}} = 0\end{aligned}

 For solving force-balance by adjusting :math:`P` and :math:`f`
profiles, the form used is

.. math::

   \begin{aligned}
   \ensuremath{B_{\text{tor}}}\ensuremath{h_\theta}\ensuremath{\frac{\partial \ensuremath{B_{\text{tor}}}}{\partial x}} + \frac{\ensuremath{B_{\text{tor}}}^2\ensuremath{h_\theta}}{R}\ensuremath{\frac{\partial R}{\partial x}} +
   \mu_0\ensuremath{h_\theta}\ensuremath{\frac{\partial P}{\partial x}} = -\ensuremath{B_{\text{pol}}}\ensuremath{\frac{\partial }{\partial x}}\left(\ensuremath{B_{\text{pol}}}\ensuremath{h_\theta}\right)\end{aligned}

 A quick way to calculate f is to rearrange this to:

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial \ensuremath{B_{\text{tor}}}}{\partial x}} = \ensuremath{B_{\text{tor}}}\left[-\frac{1}{R}\ensuremath{\frac{\partial R}{\partial x}}\right] +
   \frac{1}{\ensuremath{B_{\text{tor}}}}\left[-\mu_0\ensuremath{\frac{\partial P}{\partial x}} -
   \ensuremath{\frac{\partial \ensuremath{B_{\text{pol}}}}{\partial \ensuremath{h_\theta}}}\ensuremath{\frac{\partial }{\partial x}}\left(\ensuremath{B_{\text{pol}}}\ensuremath{h_\theta}\right)\right]\end{aligned}

 and then integrate this using LSODE.

Parallel current
----------------

.. math::

   \begin{aligned}
   J_{||} = \ensuremath{\boldsymbol{b}}\cdot\ensuremath{\boldsymbol{J}}\qquad b^y = \frac{\ensuremath{B_{\text{pol}}}}{B\ensuremath{h_\theta}}\end{aligned}

 and from equation `[eq:curlb_y] <#eq:curlb_y>`__:

.. math::

   \begin{aligned}
   J_y = \frac{\ensuremath{\sigma_{B\theta}}}{\mu_0}\left\{-\frac{B^2\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\ensuremath{\frac{\partial }{\partial x}}\left(\ensuremath{B_{\text{tor}}}R\right) + \ensuremath{B_{\text{tor}}}
   R\left[\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B^2\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) - \ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}R\ensuremath{\frac{\partial \nu}{\partial x}}\right]\right\}\end{aligned}

 since :math:`J_{||} = b^yJ_y`,

.. math::

   \begin{aligned}
   \mu_0 J_{||} =\ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}\ensuremath{B_{\text{tor}}}
   R}{B\ensuremath{h_\theta}}\left[\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B^2\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) - \ensuremath{B_{\text{tor}}}R\ensuremath{\frac{\partial \nu}{\partial x}}\right] -
   \ensuremath{\sigma_{B\theta}}B\ensuremath{\frac{\partial }{\partial x}}\left(\ensuremath{B_{\text{tor}}}R\right)\end{aligned}

Curvature
---------

For reduced MHD, need to calculate curvature term
:math:`\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}`,
where
:math:`\ensuremath{\boldsymbol{\kappa}} = \left(\ensuremath{\boldsymbol{b}}\cdot\nabla\right)\ensuremath{\boldsymbol{b}}=
-\ensuremath{\boldsymbol{b}}\times\left(\nabla\times\ensuremath{\boldsymbol{b}}\right)`.
Re-arranging, this becomes:

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}} = \nabla\times\ensuremath{\boldsymbol{b}}-
   \ensuremath{\boldsymbol{b}}\left(\ensuremath{\boldsymbol{b}}\cdot\left(\nabla\times\ensuremath{\boldsymbol{b}}\right)\right)\end{aligned}

 Components of :math:`\nabla\times\ensuremath{\boldsymbol{b}}` are:

.. math::

   \begin{aligned}
   \left(\nabla\times\ensuremath{\boldsymbol{b}}\right)^x =& \ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial y}}\left(\frac{\ensuremath{B_{\text{tor}}}
   R}{B}\right) \\ \left(\nabla\times\ensuremath{\boldsymbol{b}}\right)^y =&
       -\ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{B_{\text{tor}}}R}{B}\right) \\
       \left(\nabla\times\ensuremath{\boldsymbol{b}}\right)^z =&
       \frac{\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) - \ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}\ensuremath{B_{\text{tor}}}
       R}{\ensuremath{h_\theta}B}\ensuremath{\frac{\partial \nu}{\partial x}} - \ensuremath{\sigma_{B\theta}}I\frac{\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial y}}\left(\frac{\ensuremath{B_{\text{tor}}}
       R}{B}\right) \\\end{aligned}

 giving:

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{\kappa}} =& -\frac{\ensuremath{B_{\text{pol}}}}{B h_\theta}\left[\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B
   h_\theta}{\ensuremath{B_{\text{pol}}}}\right) - \ensuremath{\sigma_{B\theta}}\ensuremath{\frac{\partial }{\partial y}}\left(\frac{\ensuremath{B_{\text{tor}}}I R}{B}\right)\right]\nabla x \nonumber
   \\ &+ \ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}}{B h_\theta}\ensuremath{\frac{\partial }{\partial y}}\left(\frac{\ensuremath{B_{\text{tor}}}R}{B}\right)\nabla z
   \label{eq:curvature}\end{aligned}

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{b}}\cdot\left(\nabla\times\ensuremath{\boldsymbol{b}}\right) = -\ensuremath{\sigma_{B\theta}}B\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{B_{\text{tor}}}R}{B}\right) +
   \ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{tor}}}\ensuremath{B_{\text{pol}}}R}{B\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) -
   \frac{\ensuremath{B_{\text{pol}}}\ensuremath{B_{\text{tor}}}^2R^2}{\ensuremath{h_\theta}B^2}\ensuremath{\frac{\partial \nu}{\partial x}}\end{aligned}

 therefore,

.. math::

   \begin{aligned}
   \left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\right)^x =& \ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial y}}\left(\frac{\ensuremath{B_{\text{tor}}}
   R}{B}\right) = -\ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}\ensuremath{B_{\text{tor}}}R}{\ensuremath{h_\theta}B^2}\ensuremath{\frac{\partial B}{\partial y}} \\
   \left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\right)^y =& \frac{\ensuremath{B_{\text{pol}}}^2\ensuremath{B_{\text{tor}}}^2
   R^2}{B^3\ensuremath{h_\theta}^2}\ensuremath{\frac{\partial \nu}{\partial x}} - \ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}^2\ensuremath{B_{\text{tor}}}
   R}{B^2\ensuremath{h_\theta}^2}\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) \\
   \left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\right)^z =&
   \frac{\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) - \ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}\ensuremath{B_{\text{tor}}}
   R}{\ensuremath{h_\theta}B}\ensuremath{\frac{\partial \nu}{\partial x}} - I\left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\right)^x\end{aligned}

 Using equation \ `[eq:xbalance] <#eq:xbalance>`__:

.. math::

   \begin{aligned}
   B\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) + \frac{B\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\ensuremath{\frac{\partial B}{\partial x}} - \ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}
   R\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{B_{\text{tor}}}\ensuremath{h_\theta}}{R\ensuremath{B_{\text{pol}}}}\right) + \frac{\mu_0\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\ensuremath{\frac{\partial P}{\partial x}} =
   0\end{aligned}

 we can re-write the above components as:

.. math::

   \begin{aligned}
   \left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\right)^y =& \ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}\ensuremath{B_{\text{tor}}}
   R}{B^2\ensuremath{h_\theta}}\left[\frac{\mu_0}{B}\ensuremath{\frac{\partial P}{\partial x}} + \ensuremath{\frac{\partial B}{\partial x}}\right] \\
   \left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\right)^z =& -\frac{\mu_0}{B}\ensuremath{\frac{\partial P}{\partial x}} - \ensuremath{\frac{\partial B}{\partial x}} -
   I\left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\right)^x\end{aligned}

Curvature from div (b/B)
------------------------

The vector
:math:`\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}`
is an approximation of

.. math::

   \begin{aligned}
   \frac{B}{2}\nabla\times\left(\frac{\ensuremath{\boldsymbol{b}}}{B}\right) \simeq \ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\end{aligned}

 so can just derive from the original expression. Using the
contravariant components of :math:`\ensuremath{\boldsymbol{b}}`, and the
curl operator in curvilinear coordinates (see appendix):

.. math::

   \begin{aligned}
   \nabla\times\left(\frac{\ensuremath{\boldsymbol{b}}}{B}\right) =&
       \frac{\ensuremath{B_{\text{pol}}}}{\ensuremath{h_\theta}}\left[\left(\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) -
       \ensuremath{\frac{\partial }{\partial y}}\left(\frac{\ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}IR}{B^2}\right)\right)\ensuremath{\boldsymbol{e}}_z \right.  \\ &+
       \ensuremath{\frac{\partial }{\partial y}}\left(\frac{\ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}R}{B^2}\right)\ensuremath{\boldsymbol{e}}_x \\ &+
       \left.\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}R}{B^2}\right)\ensuremath{\boldsymbol{e}}_y\right]\end{aligned}

 This can be simplified using

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial }{\partial y}}\left(\frac{\ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}IR}{B^2}\right) = I\ensuremath{\sigma_{B\theta}}\ensuremath{B_{\text{tor}}}
   R\ensuremath{\frac{\partial }{\partial y}}\left(\frac{1}{B^2}\right) + \frac{\ensuremath{B_{\text{tor}}}R}{B^2}\ensuremath{\frac{\partial \nu}{\partial x}}\end{aligned}

 to give

.. math::

   \begin{aligned}
     \left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\right)^x =& -\ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}\ensuremath{B_{\text{tor}}}R}{\ensuremath{h_\theta}B^2}\ensuremath{\frac{\partial B}{\partial y}} \\
       \left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\right)^y =& -\ensuremath{\sigma_{B\theta}}\frac{B\ensuremath{B_{\text{pol}}}}{2\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{B_{\text{tor}}}
   R}{B^2}\right) \\ \left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\right)^z =&
       \frac{B\ensuremath{B_{\text{pol}}}}{2\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) - \frac{\ensuremath{B_{\text{pol}}}\ensuremath{B_{\text{tor}}}
       R}{2\ensuremath{h_\theta}B}\ensuremath{\frac{\partial \nu}{\partial x}} - I\left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\cdot\nabla\right)^x\end{aligned}

 The first and second terms in
:math:`\left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\cdot\nabla\right)^z`
almost cancel, so by expanding out :math:`\nu` a better expression is

.. math::

   \begin{aligned}
   \left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\right)^z = \frac{\ensuremath{B_{\text{pol}}}^3}{2\ensuremath{h_\theta}
   B}\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) - \frac{\ensuremath{B_{\text{tor}}}
   R}{2B}\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right)\end{aligned}

Curvature of a single line
--------------------------

The curvature vector can be calculated from the field-line toroidal
coordinates :math:`\left(R,Z,\phi\right)` as follows. The line element
is given by

.. math::

   \begin{aligned}
   d\ensuremath{\boldsymbol{r}} = dR\ensuremath{\hat{\ensuremath{\boldsymbol{R}}}}+ dZ\ensuremath{\hat{\ensuremath{\boldsymbol{Z}}}}+ Rd\phi\ensuremath{\hat{\ensuremath{\boldsymbol{\phi}}}}\end{aligned}

 Hence the tangent vector is

.. math::

   \begin{aligned}
   \hat{\ensuremath{\boldsymbol{T}}} \equiv \ensuremath{\frac{d \ensuremath{\boldsymbol{r}}}{d s}} = \ensuremath{\frac{d R}{d s}}\ensuremath{\hat{\ensuremath{\boldsymbol{R}}}}+ \ensuremath{\frac{d Z}{d s}}\ensuremath{\hat{\ensuremath{\boldsymbol{Z}}}}+
   R\ensuremath{\frac{d \phi}{d s}}\ensuremath{\hat{\ensuremath{\boldsymbol{\phi}}}}\end{aligned}

 where :math:`s` is the distance along the field-line. From this, the
curvature vector is given by

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{\kappa}}\equiv \ensuremath{\frac{d \ensuremath{\boldsymbol{T}}}{d s}} =& \ensuremath{\frac{d^2 R}{d s^2}}\ensuremath{\hat{\ensuremath{\boldsymbol{R}}}}+ \ensuremath{\frac{d R}{d s}}\ensuremath{\frac{d \phi}{d s}}\ensuremath{\hat{\ensuremath{\boldsymbol{\phi}}}}
       \\ &+ \ensuremath{\frac{d^2 Z}{d s^2}}\ensuremath{\hat{\ensuremath{\boldsymbol{Z}}}}\\ &+ \ensuremath{\frac{d R}{d s}}\ensuremath{\frac{d \phi}{d s}}\ensuremath{\hat{\ensuremath{\boldsymbol{\phi}}}}+
       R\ensuremath{\frac{d^2 \phi}{d s^2}}\ensuremath{\hat{\ensuremath{\boldsymbol{\phi}}}}- R\left(\ensuremath{\frac{d \phi}{d s}}\right)^2 \ensuremath{\hat{\ensuremath{\boldsymbol{R}}}}\end{aligned}

 i.e.

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{\kappa}}= \left[\ensuremath{\frac{d^2 R}{d s^2}} - R\left(\ensuremath{\frac{d \phi}{d s}}\right)^2\right]\ensuremath{\hat{\ensuremath{\boldsymbol{R}}}}+ \ensuremath{\frac{d^2 Z}{d s^2}}\ensuremath{\hat{\ensuremath{\boldsymbol{Z}}}}+
   \left[2\ensuremath{\frac{d R}{d s}}\ensuremath{\frac{d \phi}{d s}} + R\ensuremath{\frac{d^2 \phi}{d s^2}}\right]\ensuremath{\hat{\ensuremath{\boldsymbol{\phi}}}}
   \label{eq:kappaline}\end{aligned}

 Want the components of
:math:`\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}`,
and since the vector :math:`\ensuremath{\boldsymbol{b}}` is just the
tangent vector :math:`\ensuremath{\boldsymbol{T}}` above, this can be
written using the cross-products

.. math::

   \begin{aligned}
   \ensuremath{\hat{\ensuremath{\boldsymbol{R}}}}\times\ensuremath{\hat{\ensuremath{\boldsymbol{Z}}}}= -\ensuremath{\hat{\ensuremath{\boldsymbol{\phi}}}}\qquad \ensuremath{\hat{\ensuremath{\boldsymbol{\phi}}}}\times\ensuremath{\hat{\ensuremath{\boldsymbol{Z}}}}= \ensuremath{\hat{\ensuremath{\boldsymbol{R}}}}\qquad
   \ensuremath{\hat{\ensuremath{\boldsymbol{R}}}}\times\ensuremath{\hat{\ensuremath{\boldsymbol{\phi}}}}= \ensuremath{\hat{\ensuremath{\boldsymbol{Z}}}}\end{aligned}

 This vector must then be dotted with :math:`\nabla\psi`,
:math:`\nabla\theta`, and :math:`\nabla\phi`. This is done by writing
these vectors in cylindrical coordinates:

.. math::

   \begin{aligned}
   \nabla\psi =& \ensuremath{\frac{\partial \psi}{\partial R}}\hat{\ensuremath{\boldsymbol{R}}} + \ensuremath{\frac{\partial \psi}{\partial Z}}\hat{\ensuremath{\boldsymbol{Z}}} \\ \nabla\theta =&
       \frac{1}{\ensuremath{B_{\text{pol}}}\ensuremath{h_\theta}}\nabla\phi\times\nabla\psi =
       \frac{1}{R\ensuremath{B_{\text{pol}}}\ensuremath{h_\theta}}\left(\ensuremath{\frac{\partial \psi}{\partial Z}}\hat{\ensuremath{\boldsymbol{R}}} - \ensuremath{\frac{\partial \psi}{\partial R}}\hat{\ensuremath{\boldsymbol{Z}}}\right) \\\end{aligned}

 An alternative is to use

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{b}}\times \nabla\phi = \frac{\ensuremath{\sigma_{B\theta}}}{BR^2}\nabla\psi\end{aligned}

 and that the tangent vector
:math:`\ensuremath{\boldsymbol{T}} = \ensuremath{\boldsymbol{b}}`. This
gives

.. math::

   \begin{aligned}
   \nabla\psi = \ensuremath{\sigma_{B\theta}}BR\left[\frac{dR}{ds}\ensuremath{\boldsymbol{Z}} - \frac{dZ}{ds}\ensuremath{\boldsymbol{R}}\right]
   \label{eq:flinenablapsi}\end{aligned}

 and so because
:math:`d\phi / ds = \ensuremath{B_{\text{tor}}}/ \left(RB\right)`

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{\kappa}}\cdot\nabla\psi = \ensuremath{\sigma_{B\theta}}BR\left[ \left( \frac{\ensuremath{B_{\text{tor}}}^2}{RB^2} -
   \ensuremath{\frac{d^2 R}{d s^2}}\right)\ensuremath{\frac{d Z}{d s}} + \ensuremath{\frac{d^2 Z}{d s^2}}\frac{dR}{ds} \right]
   \label{eq:flinekappsi}\end{aligned}

 Taking the cross-product of the tangent vector with the curvature in
equation \ `[eq:kappaline] <#eq:kappaline>`__ above gives

.. math::

   \begin{aligned}
     \ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}=& \left[\frac{\ensuremath{B_{\text{tor}}}}{B}\ensuremath{\frac{d^2 Z}{d s^2}} -
   \ensuremath{\frac{d Z}{d s}}\left(2\ensuremath{\frac{d R}{d s}}\ensuremath{\frac{d \phi}{d s}} + R\ensuremath{\frac{d^2 \phi}{d s^2}}\right)\right]\ensuremath{\boldsymbol{R}} \\ &+
       \left[\ensuremath{\frac{d R}{d s}}\left(2\ensuremath{\frac{d R}{d s}}\ensuremath{\frac{d \phi}{d s}} + R\ensuremath{\frac{d^2 \phi}{d s^2}}\right) -
       \frac{\ensuremath{B_{\text{tor}}}}{B}\left(\ensuremath{\frac{d^2 R}{d s^2}} - R\left(\ensuremath{\frac{d \phi}{d s}}\right)^2\right)\right]\ensuremath{\boldsymbol{Z}} \\ &+
           \left[\ensuremath{\frac{d Z}{d s}}\left(\ensuremath{\frac{d^2 R}{d s^2}} - R\left(\ensuremath{\frac{d \phi}{d s}}\right)^2\right) -
           \ensuremath{\frac{d R}{d s}}\ensuremath{\frac{d^2 Z}{d s^2}}\right]\ensuremath{\hat{\ensuremath{\boldsymbol{\phi}}}}\end{aligned}

 The components in field-aligned coordinates can then be calculated:

.. math::

   \begin{aligned}
   \left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\right)^x =& \ensuremath{\sigma_{B\theta}}\left(\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\right)\cdot\nabla\psi \\ =&
       \frac{R\ensuremath{B_{\text{pol}}}^2}{B}\left(2\ensuremath{\frac{d R}{d s}}\ensuremath{\frac{d \phi}{d s}} + R\ensuremath{\frac{d^2 \phi}{d s^2}}\right) -
       R\ensuremath{B_{\text{tor}}}\left(\ensuremath{\frac{d R}{d s}}\ensuremath{\frac{d^2 R}{d s^2}} + \ensuremath{\frac{d Z}{d s}}\ensuremath{\frac{d^2 Z}{d s^2}}\right) +
       \frac{\ensuremath{B_{\text{tor}}}^3}{B^2}\ensuremath{\frac{d R}{d s}}\end{aligned}

Curvature in toroidal coordinates
---------------------------------

In toroidal coordinates :math:`\left(\psi,\theta,\phi\right)`, the
:math:`\ensuremath{\boldsymbol{b}}` vector is

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{b}}=& \frac{\ensuremath{B_{\text{pol}}}}{B}\ensuremath{\hat{\ensuremath{\boldsymbol{e}}}}_\theta + \frac{\ensuremath{B_{\text{tor}}}}{B}\ensuremath{\hat{\ensuremath{\boldsymbol{e}}}}_\phi \\ =&
       \frac{\ensuremath{B_{\text{pol}}}\ensuremath{h_\theta}}{B}\nabla\theta + \frac{R\ensuremath{B_{\text{tor}}}}{B}\nabla\phi\end{aligned}

 The curl of this vector is

.. math::

   \begin{aligned}
   \left(\nabla\times\ensuremath{\boldsymbol{b}}\right)^\psi =& \frac{1}{\sqrt{g}}\left(\ensuremath{\frac{\partial b_\phi}{\partial \theta}} -
       \ensuremath{\frac{\partial b_\theta}{\partial \phi}}\right) \\ \left(\nabla\times\ensuremath{\boldsymbol{b}}\right)^\theta =&
       \frac{1}{\sqrt{g}}\left(\ensuremath{\frac{\partial b_\psi}{\partial \phi}} - \ensuremath{\frac{\partial b_\phi}{\partial \psi}}\right) \\
       \left(\nabla\times\ensuremath{\boldsymbol{b}}\right)^\phi =& \frac{1}{\sqrt{g}}\left(\ensuremath{\frac{\partial b_\theta}{\partial \psi}}
       - \ensuremath{\frac{\partial b_\psi}{\partial \theta}}\right)\end{aligned}

 where
:math:`1/\sqrt{g} = \ensuremath{B_{\text{pol}}}/\ensuremath{h_\theta}`.
Therefore, in terms of unit vectors:

.. math::

   \begin{aligned}
   \nabla\times\ensuremath{\boldsymbol{b}}=
   \frac{1}{R\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial \theta}}\left(\frac{R\ensuremath{B_{\text{tor}}}}{B}\right)\ensuremath{\hat{\ensuremath{\boldsymbol{e}}}}_\psi -
   \ensuremath{B_{\text{pol}}}\ensuremath{\frac{\partial }{\partial \psi}}\left(\frac{R\ensuremath{B_{\text{tor}}}}{B}\right)\ensuremath{\hat{\ensuremath{\boldsymbol{e}}}}_\theta + \frac{\ensuremath{B_{\text{pol}}}
   R}{\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial \psi}}\left(\frac{\ensuremath{h_\theta}\ensuremath{B_{\text{pol}}}}{B}\right)\ensuremath{\hat{\ensuremath{\boldsymbol{e}}}}_\phi\end{aligned}

psi derivative of the B field
-----------------------------

Needed to calculate magnetic shear, and one way to get the curvature.
The simplest way is to use finite differencing, but there is another way
using local derivatives (implemented using DCT).

.. math::

   \begin{aligned}
   \ensuremath{B_{\text{pol}}}= \frac{\left|\nabla\psi\right|}{R} = \frac{1}{R}\sqrt{\left(\ensuremath{\frac{\partial \psi}{\partial R}}\right)^2 +
   \left(\ensuremath{\frac{\partial \psi}{\partial R}}\right)^2}\end{aligned}

 Using

.. math::

   \begin{aligned}
   \nabla\ensuremath{B_{\text{pol}}}= \ensuremath{\frac{\partial \ensuremath{B_{\text{pol}}}}{\partial \psi}}\nabla\psi + \ensuremath{\frac{\partial \ensuremath{B_{\text{pol}}}}{\partial \theta}}\nabla\theta +
   \ensuremath{\frac{\partial \ensuremath{B_{\text{pol}}}}{\partial \phi}}\nabla\phi\end{aligned}

 we get

.. math::

   \begin{aligned}
   \nabla\ensuremath{B_{\text{pol}}}\cdot\nabla\psi = \ensuremath{\frac{\partial \ensuremath{B_{\text{pol}}}}{\partial \psi}}\left|\nabla\psi\right|^2\end{aligned}

 and so

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial \ensuremath{B_{\text{pol}}}}{\partial \psi}} = \nabla\ensuremath{B_{\text{pol}}}\cdot\nabla\psi / \left(R\ensuremath{B_{\text{pol}}}\right)^2\end{aligned}

 The derivatives of :math:`\ensuremath{B_{\text{pol}}}` in :math:`R` and
:math:`Z` are:

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial \ensuremath{B_{\text{pol}}}}{\partial R}} =& -\frac{\ensuremath{B_{\text{pol}}}}{R} + \frac{1}{\ensuremath{B_{\text{pol}}}
   R^2}\left[\ensuremath{\frac{\partial \psi}{\partial R}}\ensuremath{\frac{\partial^2 \psi}{\partial {R}^2}} +
   \ensuremath{\frac{\partial \psi}{\partial Z}}\frac{\partial^2\psi}{\partial R\partial Z}\right] \\ \ensuremath{\frac{\partial \ensuremath{B_{\text{pol}}}}{\partial Z}}
   =& \frac{1}{\ensuremath{B_{\text{pol}}}R^2}\left[\ensuremath{\frac{\partial \psi}{\partial Z}}\ensuremath{\frac{\partial^2 \psi}{\partial {Z}^2}} +
   \ensuremath{\frac{\partial \psi}{\partial R}}\frac{\partial^2\psi}{\partial R\partial Z}\right]\end{aligned}

 For the toroidal field, :math:`\ensuremath{B_{\text{tor}}}= f/R`

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial \ensuremath{B_{\text{tor}}}}{\partial \psi}} = \frac{1}{R}\ensuremath{\frac{\partial f}{\partial \psi}} - \frac{f}{R^2}\ensuremath{\frac{\partial R}{\partial \psi}}\end{aligned}

 As above,
:math:`\ensuremath{\frac{\partial R}{\partial \psi}} = \nabla R \cdot\nabla\psi / \left(R\ensuremath{B_{\text{pol}}}\right)^2`,
and since :math:`\nabla R\cdot\nabla R = 1`,

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial R}{\partial \psi}} = \ensuremath{\frac{\partial \psi}{\partial R}} / \left(R\ensuremath{B_{\text{pol}}}\right)^2\end{aligned}

 similarly,

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial Z}{\partial \psi}} = \ensuremath{\frac{\partial \psi}{\partial Z}} / \left(R\ensuremath{B_{\text{pol}}}\right)^2\end{aligned}

 and so the variation of toroidal field with :math:`\psi` is

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial \ensuremath{B_{\text{tor}}}}{\partial \psi}} = \frac{1}{R}\ensuremath{\frac{\partial f}{\partial \psi}} -
   \frac{\ensuremath{B_{\text{tor}}}}{R^3\ensuremath{B_{\text{pol}}}^2}\ensuremath{\frac{\partial \psi}{\partial R}}\end{aligned}

 From the definition
:math:`B=\sqrt{\ensuremath{B_{\text{tor}}}^2 + \ensuremath{B_{\text{pol}}}^2}`,

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial B}{\partial \psi}} = \frac{1}{B}\left(\ensuremath{B_{\text{tor}}}\ensuremath{\frac{\partial \ensuremath{B_{\text{tor}}}}{\partial \psi}} + \ensuremath{B_{\text{pol}}}\ensuremath{\frac{\partial \ensuremath{B_{\text{pol}}}}{\partial \psi}}\right)\end{aligned}

Parallel derivative of the B field
----------------------------------

To get the parallel nablaients of the :math:`B` field components, start
with

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial }{\partial s}}\left(B^2\right) = \ensuremath{\frac{\partial }{\partial s}}\left(\ensuremath{B_{\text{tor}}}^2\right) + \ensuremath{\frac{\partial }{\partial s}}\left(\ensuremath{B_{\text{pol}}}^2\right)\end{aligned}

 Using the fact that :math:`R\ensuremath{B_{\text{tor}}}` is constant
along :math:`s`,

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial }{\partial s}}\left(R^2\ensuremath{B_{\text{tor}}}^2\right) = R^2\ensuremath{\frac{\partial }{\partial s}}\left(\ensuremath{B_{\text{tor}}}^2\right) +
   \ensuremath{B_{\text{tor}}}^2\ensuremath{\frac{\partial }{\partial s}}\left(R^2\right) = 0\end{aligned}

 which gives

.. math::

   \begin{aligned}
     \ensuremath{\frac{\partial }{\partial s}}\left(\ensuremath{B_{\text{tor}}}^2\right) = -\frac{\ensuremath{B_{\text{tor}}}^2}{R^2}\ensuremath{\frac{\partial }{\partial s}}\left(R^2\right)\end{aligned}

 The poloidal field can be calculated from

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial }{\partial s}}\left(\nabla\psi \cdot \nabla\psi\right) = \ensuremath{\frac{\partial }{\partial s}}\left(R^2\ensuremath{B_{\text{pol}}}^2\right) =
   R^2\ensuremath{\frac{\partial }{\partial s}}\left(\ensuremath{B_{\text{pol}}}^2\right) + \ensuremath{B_{\text{pol}}}^2\ensuremath{\frac{\partial }{\partial s}}\left(R^2\right)\end{aligned}

 Using equation \ `[eq:flinenablapsi] <#eq:flinenablapsi>`__,
:math:`\nabla\psi \cdot \nabla\psi` can also be written as

.. math::

   \begin{aligned}
   \nabla\psi \cdot \nabla\psi = B^2R^2\left[\left(\ensuremath{\frac{\partial R}{\partial s}}\right)^2 +
   \left(\ensuremath{\frac{\partial Z}{\partial s}}\right)^2\right]\end{aligned}

 and so (unsurprisingly)

.. math::

   \begin{aligned}
   \frac{\ensuremath{B_{\text{pol}}}^2}{B^2} = \left[\left(\ensuremath{\frac{\partial R}{\partial s}}\right)^2 + \left(\ensuremath{\frac{\partial Z}{\partial s}}\right)^2\right]\end{aligned}

 Hence

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial }{\partial s}}\left(\ensuremath{B_{\text{pol}}}^2\right) = B^2\ensuremath{\frac{\partial }{\partial s}}\left[\left(\ensuremath{\frac{\partial R}{\partial s}}\right)^2 +
   \left(\ensuremath{\frac{\partial Z}{\partial s}}\right)^2\right] + \frac{\ensuremath{B_{\text{pol}}}^2}{B^2}\ensuremath{\frac{\partial }{\partial s}}\left(B^2\right)\end{aligned}

 Which gives

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial }{\partial s}}\left(B^2\right) = -\frac{B^2}{R^2}\ensuremath{\frac{\partial }{\partial s}}\left(R^2\right) +
   \frac{B^4}{\ensuremath{B_{\text{tor}}}^2}\ensuremath{\frac{\partial }{\partial s}}\left[\left(\ensuremath{\frac{\partial R}{\partial s}}\right)^2 + \left(\ensuremath{\frac{\partial Z}{\partial s}}\right)^2\right]\end{aligned}

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial }{\partial s}}\left(\ensuremath{B_{\text{pol}}}^2\right) = \left(1 +
   \frac{\ensuremath{B_{\text{pol}}}^2}{\ensuremath{B_{\text{tor}}}^2}\right)B^2\ensuremath{\frac{\partial }{\partial s}}\left[\left(\ensuremath{\frac{\partial R}{\partial s}}\right)^2 +
   \left(\ensuremath{\frac{\partial Z}{\partial s}}\right)^2\right] - \frac{\ensuremath{B_{\text{pol}}}^2}{R^2}\ensuremath{\frac{\partial }{\partial s}}\left(R^2\right)\end{aligned}

Magnetic shear from J x B
-------------------------

Re-arranging the radial force balance
equation \ `[eq:xbalance] <#eq:xbalance>`__ gives

.. math::

   \begin{aligned}
   \frac{\ensuremath{B_{\text{pol}}}^2R}{\ensuremath{B_{\text{tor}}}}\ensuremath{\frac{\partial \nu}{\partial \psi}} + \nu\left(\frac{2RB}{\ensuremath{B_{\text{tor}}}}\ensuremath{\frac{\partial B}{\partial \psi}} +
   \frac{B^2}{\ensuremath{B_{\text{tor}}}}\ensuremath{\frac{\partial R}{\partial \psi}} - \frac{B^2R}{\ensuremath{B_{\text{tor}}}^2}\ensuremath{\frac{\partial \ensuremath{B_{\text{tor}}}}{\partial \psi}}\right) +
   \frac{\mu_0\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\ensuremath{\frac{\partial P}{\partial \psi}} = 0\end{aligned}

Magnetic shear
--------------

The field-line pitch is given by

.. math::

   \begin{aligned}
   \nu = \frac{\ensuremath{h_\theta}\ensuremath{B_{\text{tor}}}}{\ensuremath{B_{\text{pol}}}R}\end{aligned}

 and so

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial \nu}{\partial \psi}} = \frac{\nu}{\ensuremath{h_\theta}}\ensuremath{\frac{\partial \ensuremath{h_\theta}}{\partial \psi}} +
   \frac{\nu}{\ensuremath{B_{\text{tor}}}}\ensuremath{\frac{\partial \ensuremath{B_{\text{tor}}}}{\partial \psi}} - \frac{\nu}{\ensuremath{B_{\text{pol}}}}\ensuremath{\frac{\partial \ensuremath{B_{\text{pol}}}}{\partial \psi}} -
   \frac{\nu}{R}\ensuremath{\frac{\partial R}{\partial \psi}}\end{aligned}

 The last three terms are given in the previous section, but
:math:`\partial\ensuremath{h_\theta}/\partial\psi` needs to be evaluated

psi derivative of h
-------------------

From the expression for curvature `[eq:curvature] <#eq:curvature>`__,
and using
:math:`\nabla x \cdot \nabla \psi = \ensuremath{\sigma_{B\theta}}\left(R\ensuremath{B_{\text{pol}}}\right)^2`
and
:math:`\nabla z\cdot\nabla \psi = -\ensuremath{\sigma_{B\theta}}I \left(R\ensuremath{B_{\text{pol}}}\right)^2`

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{\kappa}}\cdot\nabla\psi =& -\ensuremath{\sigma_{B\theta}}
       \frac{\ensuremath{B_{\text{pol}}}}{B\ensuremath{h_\theta}}\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\left[\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) -
       \ensuremath{\sigma_{B\theta}}\ensuremath{\frac{\partial }{\partial y}}\left(\frac{\ensuremath{B_{\text{tor}}}IR}{B}\right)\right] \\ &- I\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}
           \frac{\ensuremath{B_{\text{pol}}}}{B\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial y}}\left(\frac{\ensuremath{B_{\text{tor}}}R}{B}\right)\end{aligned}

 The second and third terms partly cancel, and using
:math:`\ensuremath{\frac{\partial I}{\partial y}} = \ensuremath{\sigma_{B\theta}}
\ensuremath{\frac{\partial \nu}{\partial x}}`

.. math::

   \begin{aligned}
     \frac{\ensuremath{\boldsymbol{\kappa}}\cdot\nabla\psi}{\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}} =&
       -\ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}}{B\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) +
       \ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}}{B\ensuremath{h_\theta}}\frac{\ensuremath{B_{\text{tor}}}R}{B}\ensuremath{\frac{\partial \nu}{\partial x}} \\ =&
       -\ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}}{B\ensuremath{h_\theta}}\left[\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}}\right) - \frac{\ensuremath{B_{\text{tor}}}
       R}{B}\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{B_{\text{tor}}}\ensuremath{h_\theta}}{\ensuremath{B_{\text{pol}}}R}\right)\right] \\ =&
               -\ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}}{B\ensuremath{h_\theta}}\left[\ensuremath{h_\theta}\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B}{\ensuremath{B_{\text{pol}}}}\right) -
               \ensuremath{h_\theta}\frac{\ensuremath{B_{\text{tor}}}R}{B}\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{B_{\text{tor}}}}{\ensuremath{B_{\text{pol}}}R}\right) +
           \frac{B^2}{B\ensuremath{B_{\text{pol}}}}\ensuremath{\frac{\partial \ensuremath{h_\theta}}{\partial x}} -
       \frac{\ensuremath{B_{\text{tor}}}^2}{B\ensuremath{B_{\text{pol}}}}\ensuremath{\frac{\partial \ensuremath{h_\theta}}{\partial x}}\right] \\ =& -\ensuremath{\sigma_{B\theta}}
           \frac{\ensuremath{B_{\text{pol}}}}{B^2\ensuremath{h_\theta}}\ensuremath{\frac{\partial \ensuremath{h_\theta}}{\partial x}} -
           \ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}}{B^2}\left[B\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B}{\ensuremath{B_{\text{pol}}}}\right) - \ensuremath{B_{\text{tor}}}
           R\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{B_{\text{tor}}}}{\ensuremath{B_{\text{pol}}}R}\right)\right]\end{aligned}

 Writing

.. math::

   \begin{aligned}
   B\ensuremath{\frac{\partial }{\partial x}}\left(\frac{B}{\ensuremath{B_{\text{pol}}}}\right) =& \ensuremath{\frac{\partial }{\partial x}}\left(\frac{B^2}{\ensuremath{B_{\text{pol}}}}\right) -
       \frac{B}{\ensuremath{B_{\text{pol}}}}\ensuremath{\frac{\partial B}{\partial x}} \\ \ensuremath{B_{\text{tor}}}R\ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{B_{\text{tor}}}}{\ensuremath{B_{\text{pol}}}R}\right) =&
       \ensuremath{\frac{\partial }{\partial x}}\left(\frac{\ensuremath{B_{\text{tor}}}^2}{\ensuremath{B_{\text{pol}}}}\right) - \frac{\ensuremath{B_{\text{tor}}}}{\ensuremath{B_{\text{pol}}}R}\ensuremath{\frac{\partial }{\partial x}}\left(\ensuremath{B_{\text{tor}}}
       R\right)\end{aligned}

 and using
:math:`B\ensuremath{\frac{\partial B}{\partial x}} = \ensuremath{B_{\text{tor}}}\ensuremath{\frac{\partial \ensuremath{B_{\text{tor}}}}{\partial x}} + \ensuremath{B_{\text{pol}}}\ensuremath{\frac{\partial \ensuremath{B_{\text{pol}}}}{\partial x}}`,
this simplifies to give

.. math::

   \begin{aligned}
   \frac{\ensuremath{\boldsymbol{\kappa}}\cdot\nabla\psi}{\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}} =
   -\ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{pol}}}^2}{B^2\ensuremath{h_\theta}}\ensuremath{\frac{\partial \ensuremath{h_\theta}}{\partial x}} - \ensuremath{\sigma_{B\theta}}\frac{\ensuremath{B_{\text{tor}}}^2}{B^2
   R}\ensuremath{\frac{\partial R}{\partial x}}
   \label{eq:dhdpsi}\end{aligned}

 This can be transformed into an expression for
:math:`\ensuremath{\frac{\partial \ensuremath{h_\theta}}{\partial x}}`
involving only derivatives along field-lines. Writing :math:`\nabla R =
\ensuremath{\frac{\partial R}{\partial \psi}}\nabla\psi + \ensuremath{\frac{\partial R}{\partial \theta}}\nabla\theta`,

.. math::

   \begin{aligned}
   \nabla R \cdot \nabla\psi = \ensuremath{\frac{\partial R}{\partial \psi}}\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\end{aligned}

 Using `[eq:flinenablapsi] <#eq:flinenablapsi>`__,

.. math::

   \begin{aligned}
   \nabla\psi \cdot \nabla R = -\ensuremath{\sigma_{B\theta}}B R\frac{dZ}{ds}\end{aligned}

 and so

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial R}{\partial x}} = -\frac{BR}{\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}}\frac{dZ}{ds}\end{aligned}

 Substituting this and equation `[eq:flinekappsi] <#eq:flinekappsi>`__
for :math:`\ensuremath{\boldsymbol{\kappa}}\cdot\nabla\psi` into
equation \ `[eq:dhdpsi] <#eq:dhdpsi>`__ the
:math:`\ensuremath{\frac{\partial R}{\partial x}}` term cancels with
part of the :math:`\ensuremath{\boldsymbol{\kappa}}\cdot\nabla\psi`
term, simplifying to

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial \ensuremath{h_\theta}}{\partial x}} =
   -\ensuremath{h_\theta}\frac{B^3R}{\ensuremath{B_{\text{pol}}}^2\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}}\left[\frac{d^2Z}{ds^2}\frac{dR}{ds} -
   \frac{d^2R}{ds^2}\frac{dZ}{ds}\right]\end{aligned}

.. _sec:shiftcoords:

Shifted radial derivatives
==========================

The coordinate system given by
equation \ `[eq:coordtransform] <#eq:coordtransform>`__ and used in the
above sections has a problem: There is a special poloidal location
:math:`\theta_0` where the radial basis vector
:math:`\ensuremath{\boldsymbol{e}}_x` is purely in the
:math:`\nabla\psi` direction. Moving away from this location, the
coordinate system becomes sheared in the toroidal direction.

Making the substitution

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial }{\partial x}} = \ensuremath{\frac{\partial }{\partial \psi}} + I\ensuremath{\frac{\partial }{\partial z}}\end{aligned}

 we also get the mixed derivative

.. math::

   \begin{aligned}
   \frac{\partial}{\partial z\partial x} =& \ensuremath{\frac{\partial }{\partial z}}\ensuremath{\frac{\partial }{\partial \psi}} +
       \ensuremath{\frac{\partial I}{\partial z}}\ensuremath{\frac{\partial }{\partial z}} + I\frac{\partial^2}{\partial z^2} \nonumber \\ =&
       \frac{\partial^2}{\partial z\partial \psi} + I\frac{\partial^2}{\partial
       z^2}\end{aligned}

 and second-order :math:`x` derivative

.. math::

   \begin{aligned}
   \frac{\partial^2}{\partial x^2} =& \frac{\partial^2}{\partial \psi^2} +
       \ensuremath{\frac{\partial }{\partial \psi}}\left(I\ensuremath{\frac{\partial }{\partial z}}\right) + I\ensuremath{\frac{\partial }{\partial z}}\left(\ensuremath{\frac{\partial }{\partial \psi}} +
       I\ensuremath{\frac{\partial }{\partial z}}\right) \nonumber \\ =& \frac{\partial^2}{\partial \psi^2} +
       I^2\frac{\partial^2}{\partial z^2} + 2I\frac{\partial^2}{\partial z\partial
       \psi} + \ensuremath{\frac{\partial I}{\partial \psi}}\ensuremath{\frac{\partial }{\partial z}}\end{aligned}

Perpendicular Laplacian
-----------------------

.. math::

   \begin{aligned}
   \nabla_\perp^2= \ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\left[\ensuremath{\frac{\partial^2 }{\partial {x}^2}} - 2I\frac{\partial^2}{\partial z\partial x} +
   \left(I^2 + \frac{B^2}{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^4}\right)\ensuremath{\frac{\partial^2 }{\partial {z}^2}}\right]\end{aligned}

 transforms to

.. math::

   \begin{aligned}
   \nabla_\perp^2= \ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\left[\ensuremath{\frac{\partial^2 }{\partial {\psi}^2}} + \ensuremath{\frac{\partial I}{\partial \psi}}\ensuremath{\frac{\partial }{\partial z}} +
   \frac{B^2}{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^4}\ensuremath{\frac{\partial^2 }{\partial {z}^2}}\right]
   \label{eq:delp}\end{aligned}

 The extra term involving :math:`I` disappears, but only if both the
:math:`x` and :math:`z` first derivatives are taken into account:

.. math::

   \begin{aligned}
   \nabla_\perp^2= \ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\left[\ensuremath{\frac{\partial^2 }{\partial {x}^2}} - 2I\frac{\partial^2}{\partial z\partial x} +
   \left(I^2 + \frac{B^2}{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^4}\right)\ensuremath{\frac{\partial^2 }{\partial {z}^2}}\right] + \nabla^2 x \ensuremath{\frac{\partial }{\partial x}} +
   \nabla^2 z\ensuremath{\frac{\partial }{\partial z}}\end{aligned}

 with

.. math::

   \begin{aligned}
   \nabla^2 x = \frac{1}{J}\ensuremath{\frac{\partial }{\partial x}}\left[J\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\right]\end{aligned}

.. math::

   \begin{aligned}
   \nabla^2 z =& \frac{1}{J}\left[-\ensuremath{\frac{\partial }{\partial x}}\left(JI\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\right) -
   \ensuremath{\frac{\partial }{\partial y}}\left(\frac{\ensuremath{B_{\text{tor}}}}{\ensuremath{B_{\text{pol}}}^2R}\right)\right] \nonumber \\ =&
       \frac{1}{J}\left[-I\ensuremath{\frac{\partial }{\partial x}}\left(J\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\right) - \ensuremath{\frac{\partial I}{\partial x}}J\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}-
       \ensuremath{\frac{\partial }{\partial y}}\left(\frac{\ensuremath{B_{\text{tor}}}}{\ensuremath{B_{\text{pol}}}^2R}\right)\right] \label{eq:delpz}\end{aligned}

 where :math:`J=\ensuremath{h_\theta}/ \ensuremath{B_{\text{pol}}}` is
the Jacobian. Transforming into :math:`\psi` derivatives, the middle
term of equation \ `[eq:delpz] <#eq:delpz>`__ cancels the :math:`I` term
in equation \ `[eq:delp] <#eq:delp>`__, but introduces another :math:`I`
term (first term in equation \ `[eq:delpz] <#eq:delpz>`__). This term
cancels with the :math:`\nabla^2 x` term when
:math:`\ensuremath{\frac{\partial }{\partial x}}` is expanded, so the
full expression for :math:`\nabla_\perp^2` using :math:`\psi`
derivatives is:

.. math::

   \begin{aligned}
   \nabla_\perp^2=& \ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\left[\ensuremath{\frac{\partial^2 }{\partial {\psi}^2}} + \frac{B^2}{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^4}\ensuremath{\frac{\partial^2 }{\partial {z}^2}}\right]
       \nonumber \\ &+ \frac{1}{J}\ensuremath{\frac{\partial }{\partial \psi}}\left[J\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\right]\ensuremath{\frac{\partial }{\partial \psi}} -
       \frac{1}{J}\ensuremath{\frac{\partial }{\partial y}}\left(\frac{\ensuremath{B_{\text{tor}}}}{\ensuremath{B_{\text{pol}}}^2R}\right)\ensuremath{\frac{\partial }{\partial z}}
   \label{eq:delp_shift}\end{aligned}

In orthogonal (psi, theta, zeta) flux coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For comparison, the perpendicular Laplacian can be derived in orthogonal
“flux” coordinates

.. math::

   \begin{aligned}
   \left|\nabla\psi\right| = \ensuremath{R\ensuremath{B_{\text{pol}}}}\qquad \left|\nabla\theta\right| = 1/\ensuremath{h_\theta}\qquad
   \left|\nabla\zeta\right| = 1/R\end{aligned}

 The Laplacian operator is given by

.. math::

   \begin{aligned}
   \nabla^2 A =& \ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\ensuremath{\frac{\partial^2 A}{\partial {\psi}^2}} + \frac{1}{\ensuremath{h_\theta}^2}\ensuremath{\frac{\partial^2 A}{\partial {\theta}^2}} +
       \frac{1}{R^2}\ensuremath{\frac{\partial^2 A}{\partial {\zeta}^2}} \nonumber \\ &+
       \frac{1}{J}\ensuremath{\frac{\partial }{\partial \psi}}\left[J\ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\right]\ensuremath{\frac{\partial A}{\partial \psi}} +
       \frac{1}{J}\ensuremath{\frac{\partial }{\partial \theta}}\left(J/\ensuremath{h_\theta}^2\right)\ensuremath{\frac{\partial A}{\partial \theta}}\end{aligned}

 parallel derivative by

.. math::

   \begin{aligned}
   \partial_{||} \equiv \ensuremath{\boldsymbol{b}}\cdot\nabla = \frac{\ensuremath{B_{\text{pol}}}}{B\ensuremath{h_\theta}}\ensuremath{\frac{\partial }{\partial \theta}} +
   \frac{\ensuremath{B_{\text{tor}}}}{RB}\ensuremath{\frac{\partial }{\partial \zeta}}\end{aligned}

 and so

.. math::

   \begin{aligned}
   \partial^2_{||}A \equiv \partial_{||}\left(\partial_{||}A\right) =&
       \left(\frac{\ensuremath{B_{\text{pol}}}}{B\ensuremath{h_\theta}}\right)^2\ensuremath{\frac{\partial^2 A}{\partial {\theta}^2}} +
       \left(\frac{\ensuremath{B_{\text{tor}}}}{RB}\right)^2\ensuremath{\frac{\partial^2 A}{\partial {\zeta}^2}} \nonumber \\ &+
       2\frac{\ensuremath{B_{\text{pol}}}\ensuremath{B_{\text{tor}}}}{B^2\ensuremath{h_\theta}R}\frac{\partial^2 A}{\partial\theta\partial\zeta}
       \nonumber \\ &+ \ensuremath{\frac{\partial }{\partial \theta}}\left(\frac{\ensuremath{B_{\text{pol}}}}{B\ensuremath{h_\theta}}\right)\ensuremath{\frac{\partial A}{\partial \theta}} +
       \ensuremath{\frac{\partial }{\partial \theta}}\left(\frac{\ensuremath{B_{\text{tor}}}}{RB}\right)\ensuremath{\frac{\partial A}{\partial \zeta}}\end{aligned}

 Hence in orthogonal flux coordinates, the perpendicular Laplacian is:

.. math::

   \begin{aligned}
   \nabla_\perp^2\equiv \nabla^2 - \partial_{||}^2 = \ensuremath{\left(\ensuremath{R\ensuremath{B_{\text{pol}}}}\right)^2}\left[\ensuremath{\frac{\partial^2 }{\partial {\psi}^2}} +
   \frac{1}{R^4B^2}\ensuremath{\frac{\partial^2 }{\partial {\zeta^2}^2}}\right] +
   \frac{\ensuremath{B_{\text{tor}}}^2}{\ensuremath{h_\theta}^2B^2}\ensuremath{\frac{\partial^2 }{\partial {\theta}^2}} + \cdots
   \label{eq:delp_flux}\end{aligned}

 where the neglected terms are first-order derivatives. The coefficient
for the second-order :math:`z` derivative differs from
equation \ `[eq:delp_shift] <#eq:delp_shift>`__, and
equation \ `[eq:delp_flux] <#eq:delp_flux>`__ still contains a
derivative in :math:`\theta`. This shows that the transformation made to
get equation \ `[eq:delp_shift] <#eq:delp_shift>`__ doesn’t result in
the same answer as orthogonal flux coordinates:
equation \ `[eq:delp_shift] <#eq:delp_shift>`__ is in field-aligned
coordinates.

Note that in the limit of :math:`\ensuremath{B_{\text{pol}}}= B`, both
equations \ `[eq:delp_shift] <#eq:delp_shift>`__ and
`[eq:delp_flux] <#eq:delp_flux>`__ are the same, as they should be.

Operator B x Nabla Phi Dot Nabla A
----------------------------------

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{B}}\times\nabla\phi\cdot\nabla A =& \left(\ensuremath{\frac{\partial \phi}{\partial x}}\ensuremath{\frac{\partial A}{\partial y}} -
       \ensuremath{\frac{\partial \phi}{\partial y}}\ensuremath{\frac{\partial A}{\partial x}}\right)\left(-\ensuremath{B_{\text{tor}}}\frac{\ensuremath{R\ensuremath{B_{\text{pol}}}}}{\ensuremath{h_\theta}}\right) \\ &+
       \left(\ensuremath{\frac{\partial \phi}{\partial x}}\ensuremath{\frac{\partial A}{\partial z}} - \ensuremath{\frac{\partial \phi}{\partial z}}\ensuremath{\frac{\partial A}{\partial x}}\right)\left(-B^2\right)
       \\ &- \left(\ensuremath{\frac{\partial \phi}{\partial y}}\ensuremath{\frac{\partial A}{\partial z}} -
       \ensuremath{\frac{\partial \phi}{\partial z}}\ensuremath{\frac{\partial A}{\partial y}}\right)\left(I\ensuremath{B_{\text{tor}}}\frac{\ensuremath{R\ensuremath{B_{\text{pol}}}}}{\ensuremath{h_\theta}}\right)\end{aligned}

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{B}}\times\nabla\phi\cdot\nabla A =& \left(\ensuremath{\frac{\partial \phi}{\partial \psi}}\ensuremath{\frac{\partial A}{\partial y}} + I
       \ensuremath{\frac{\partial \phi}{\partial z}}\ensuremath{\frac{\partial A}{\partial y}} - \ensuremath{\frac{\partial \phi}{\partial y}}\ensuremath{\frac{\partial A}{\partial \psi}} -
       I\ensuremath{\frac{\partial \phi}{\partial y}}\ensuremath{\frac{\partial A}{\partial z}}\right)\left(-\ensuremath{B_{\text{tor}}}\frac{\ensuremath{R\ensuremath{B_{\text{pol}}}}}{\ensuremath{h_\theta}}\right) \\ &+
       \left(\ensuremath{\frac{\partial \phi}{\partial \psi}}\ensuremath{\frac{\partial A}{\partial z}} + I\ensuremath{\frac{\partial \phi}{\partial z}}\ensuremath{\frac{\partial A}{\partial z}} -
       \ensuremath{\frac{\partial \phi}{\partial z}}\ensuremath{\frac{\partial A}{\partial \psi}} - I\ensuremath{\frac{\partial \phi}{\partial z}}\ensuremath{\frac{\partial A}{\partial z}}\right)\left(-B^2\right)
       \\ &- \left(\ensuremath{\frac{\partial \phi}{\partial y}}\ensuremath{\frac{\partial A}{\partial z}} -
       \ensuremath{\frac{\partial \phi}{\partial z}}\ensuremath{\frac{\partial A}{\partial y}}\right)\left(I\ensuremath{B_{\text{tor}}}\frac{\ensuremath{R\ensuremath{B_{\text{pol}}}}}{\ensuremath{h_\theta}}\right)\end{aligned}

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{B}}\times\nabla\phi\cdot\nabla A =& \left(\ensuremath{\frac{\partial \phi}{\partial \psi}}\ensuremath{\frac{\partial A}{\partial y}} -
       \ensuremath{\frac{\partial \phi}{\partial y}}\ensuremath{\frac{\partial A}{\partial \psi}}\right)\left(-\ensuremath{B_{\text{tor}}}\frac{\ensuremath{R\ensuremath{B_{\text{pol}}}}}{\ensuremath{h_\theta}}\right) \nonumber \\
       &+ \left(\ensuremath{\frac{\partial \phi}{\partial \psi}}\ensuremath{\frac{\partial A}{\partial z}} - \ensuremath{\frac{\partial \phi}{\partial z}}\ensuremath{\frac{\partial A}{\partial \psi}}
       \right)\left(-B^2\right)\end{aligned}

Useful identities
=================

:math:`\mathbf{b}\times\mathbf{\kappa}\cdot\nabla\psi \simeq -RB_\zeta\partial_{||}\ln B`
-----------------------------------------------------------------------------------------

Using
:math:`\mathbf{b}\times\mathbf{\kappa} \simeq \frac{B}{2}\nabla\times\frac{\mathbf{b}}{B}`,
and working in orthogonal :math:`\left(\psi, \theta, \zeta\right)`
coordinates. The magnetic field unit vector is:

.. math:: \mathbf{b} = \frac{B_\theta h_\theta}{B}\nabla\theta + \frac{B_\zeta R}{B}\nabla\zeta

 and using the definition of curl
(equation `[eq:curlcurvilinear] <#eq:curlcurvilinear>`__) we can write

.. math:: \mathbf{b}\times\mathbf{\kappa} \simeq \frac{B}{2}\nabla\times\frac{\mathbf{b}}{B} = \frac{B}{2}\frac{B_\theta}{h_\theta}\left[\frac{\partial}{\partial\theta}\left(\frac{B_\zeta R}{B^2}\right) - \frac{\partial}{\partial\zeta}\left(\frac{B_\theta h_\theta}{B^2}\right)\right]\mathbf{e}_\psi + \left[\cdot\right]\mathbf{e}_\theta + \left[\cdot\right]\mathbf{e}_\zeta

 so that when dotted with :math:`\nabla\psi`, only the first bracket
survives. The parallel gradient is

.. math:: \partial_{||} = \mathbf{b}\cdot\nabla = \frac{B_\theta}{Bh_\theta}\frac{\partial}{\partial\theta} + \frac{B_\theta}{BR}\frac{\partial}{\partial\zeta}

 Neglecting derivatives for axisymmetric equilibrium

.. math:: \frac{B}{2}\nabla\times\frac{\mathbf{b}}{B}\cdot\nabla\psi = \frac{B}{2}B\partial_{||}\left(\frac{B_\zeta R}{B^2}\right)

 Since :math:`B_\zeta R` is a flux function, this can be written as

.. math:: \frac{B}{2}\nabla\times\frac{\mathbf{b}}{B}\cdot\nabla\psi = -B_\zeta R\frac{1}{B}\partial_{||} B

 and so

.. math:: \mathbf{b}\times\mathbf{\kappa}\cdot\nabla\psi \simeq -RB_\zeta\partial_{||}\ln B

.. raw:: latex

   \bibliographystyle{unsrt}

.. raw:: latex

   \appendix

Differential geometry
=====================

| WARNING: Several mistakes have been found (and is now corrected) in
  this section, so it should be proof read before removing this warning!
| The following is notes from :raw-latex:`\cite{haeseler-1}`.

Sets of vectors :math:`\left\{\mathbf{A, B, C}\right\}` and
:math:`\left\{\mathbf{a, b, c}\right\}` are reciprocal if

.. math::

   \begin{aligned}
   \mathbf{A\cdot a} = \mathbf{B\cdot b} = \mathbf{C\cdot c} = 1\\ \mathbf{A\cdot
   b} = \mathbf{A\cdot c} = \mathbf{B\cdot a} = \mathbf{B\cdot c} = \mathbf{C\cdot
   a} = \mathbf{C\cdot b} = 0 \\\end{aligned}

 which implies that :math:`\left\{\mathbf{A, B, C}\right\}` and
:math:`\left\{\mathbf{a, b, c}\right\}` are each linearly independent.
Equivalently,

.. math::

   \begin{aligned}
   \mathbf{a} = \frac{\mathbf{B\times C}}{\mathbf{A\cdot\left(B\times C\right)}}\qquad
   \ensuremath{\boldsymbol{b}}= \frac{\mathbf{C\times A}}{\mathbf{B\cdot\left(C\times A\right)}}\qquad
   \mathbf{c} = \frac{\mathbf{A\times B}}{\mathbf{C\cdot\left(A\times B\right)}}\end{aligned}

 Either of these sets can be used as a basis, and any vector
:math:`\mathbf{w}` can be represented as
:math:`\mathbf{w} = \left(\mathbf{w\cdot a}\right)\mathbf{A} +
\left(\mathbf{w\cdot b}\right)\ensuremath{\boldsymbol{B}}+ \left(\mathbf{w\cdot c}\right)\mathbf{C}`
or as
:math:`\mathbf{w} = \left(\mathbf{w\cdot A}\right)\mathbf{a} + \left(\mathbf{w\cdot B}\right)\ensuremath{\boldsymbol{b}}
+ \left(\mathbf{w\cdot C}\right)\mathbf{c}`. In the Cartesian coordinate
system, the basis vectors are reciprocal to themselves so this
distinction is not needed. For a general set of coordinates
:math:`\left\{u^1, u^2, u^3\right\}`, tangent basis vectors can be
defined. If the Cartesian coordinates of a point are given by
:math:`\left(x, y, z\right) = \mathbf{R}\left(u^1, u^2, u^3\right)` then
the tangent basis vectors are:

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{e}}_i = \frac{\partial\mathbf{R}}{\partial u^i}\end{aligned}

 and in general these will vary from point to point. The scale factor or
metric coefficient
:math:`h_i =\left|\ensuremath{\boldsymbol{e}}_i\right|` is the distance
moved for a unit change in :math:`u^i`. The unit vector
:math:`\hat{\ensuremath{\boldsymbol{e}}}_i = \ensuremath{\boldsymbol{e}}_i/h_i`.
Definition of nabla operator:

.. raw:: latex

   \framebox{$\nabla\Phi$ of a function $\Phi$ is defined so that $d\Phi =
   \nabla\Phi\cdot d{\mathbf{R}}$}

From the chain rule,
:math:`d\mathbf{R} = \frac{\partial\mathbf{R}}{\partial u^i}du^i
= \ensuremath{\boldsymbol{e}}_idu^i` and substituting :math:`\Phi = u^i`

.. math::

   \begin{aligned}
   du^i = \nabla u^i\cdot\ensuremath{\boldsymbol{e}}_jdu^j\end{aligned}

 which can only be true if
:math:`\nabla u^i\cdot\ensuremath{\boldsymbol{e}}_j = \delta^i_j` i.e.
if

.. raw:: latex

   \framebox{Sets of vectors $\ve{e}^i\equiv\nabla u^i$ and $\ve{e}_j$ are
   reciprocal}

Since the sets of vectors
:math:`\left\{\ensuremath{\boldsymbol{e}}^i\right\}` and
:math:`\left\{\ensuremath{\boldsymbol{e}}_i\right\}` are reciprocal, any
vector :math:`\mathbf{D}` can be written as
:math:`\mathbf{D} = D_i\ensuremath{\boldsymbol{e}}^i
= D^i\ensuremath{\boldsymbol{e}}_i` where
:math:`D_i = \mathbf{D\cdot e}_i` are the covariant components and
:math:`D^i = \mathbf{D\cdot e}^i` are the contravariant components. To
convert between covariant and contravariant components, define the
metric coefficients :math:`g_{ij} = \mathbf{e_i\cdot e_j}` and
:math:`g^{ij} =
\mathbf{e^i\cdot e^j}` so that
:math:`\ensuremath{\boldsymbol{e}}_i = g_{ij}\ensuremath{\boldsymbol{e}}^j`.
:math:`g_{ij}` and :math:`g^{ij}` are symmetric and if the basis is
orthogonal then :math:`g_{ij}=g^{ij} = 0` for :math:`i\neq j` i.e. the
metric is diagonal.

.. raw:: latex

   \framebox{$g_{ij} = h_ih_j\hv{e}_i\cdot\hv{e}_j$ and so $g_{ii} = h_i^2$}

For a general set of coordinates, the nabla operator can be expressed as

.. math::

   \begin{aligned}
   \nabla = \nabla u^i\frac{\partial}{\partial u^i} =
   \ensuremath{\boldsymbol{e}}^i\frac{\partial}{\partial u^i}\end{aligned}

 and for a general set of (differentiable) coordinates
:math:`\left\{u^i\right\}`, the Laplacian is given by

.. math::

   \begin{aligned}
   \nabla^2\phi = \frac{1}{J}\frac{\partial}{\partial
   u^i}\left(Jg^{ij}\frac{\partial\phi}{\partial u^j}\right)
   \label{eq:laplacegen}\end{aligned}

 which can be expanded as

.. math::

   \begin{aligned}
   \nabla^2\phi = g^{ij}\frac{\partial^2\phi}{\partial u^i\partial u^j} +
   \underbrace{\frac{1}{J}\frac{\partial}{\partial
   u^i}\left(Jg^{ij}\right)}_{G^j}\frac{\partial\phi}{\partial u^j}
   \label{eq:laplace_expand}\end{aligned}

 where :math:`G^j` must **not** be mistaken as the so called connection
coefficients (i.e. the Christoffel symbols of second kind). Setting
:math:`\phi =
u^k` in equation (`[eq:laplacegen] <#eq:laplacegen>`__) gives
:math:`\nabla^2u^k = G^k`. Expanding
(`[eq:laplacegen] <#eq:laplacegen>`__) and setting
:math:`\left\{u^i\right\} = \left\{x, y, z\right\}` gives

.. math::

   \begin{aligned}
   \nabla^2f = \nabla\cdot\nabla f = \nabla\cdot\left(\frac{\partial}{\partial
   x}\nabla x + \frac{\partial}{\partial y}\nabla y + \frac{\partial}{\partial
   z}\nabla z\right) \nonumber \\
   \label{eq:general_laplacian}
   = \frac{\partial^2 f}{\partial x^2}\left|\nabla x\right|^2 + \frac{\partial^2
   f}{\partial y^2}\left|\nabla y\right|^2 + \frac{\partial^2 f}{\partial z^2}\left|\nabla
   z\right|^2 \\ +2\frac{\partial^2 f}{\partial x\partial y}\left(\nabla x\cdot\nabla
   y\right) +2\frac{\partial^2 f}{\partial x\partial z}\left(\nabla x\cdot\nabla z\right)
   +2\frac{\partial^2 f}{\partial y\partial z}\left(\nabla y\cdot\nabla z\right)
   \nonumber \\ +\nabla^2x\frac{\partial f}{\partial x} +\nabla^2y\frac{\partial
   f}{\partial y} + \nabla^2z\frac{\partial f}{\partial z} \nonumber\end{aligned}

 Curl defined as:

.. math::

   \begin{aligned}
   \nabla\times\mathbf{A} = \frac{1}{\sqrt{g}}\sum_k\left(\frac{\partial
   A_j}{\partial u_i} - \frac{\partial A_i}{\partial u_j}\right)\ensuremath{\boldsymbol{e}}_k \qquad i,j,k
   \texttt{ cyc } 1,2,3 \label{eq:curlcurvilinear}\end{aligned}

 Cross-product relation between contravariant and covariant vectors:

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{e}}^i = \frac{1}{J}\left(\ensuremath{\boldsymbol{e}}_j \times \ensuremath{\boldsymbol{e}}_k\right) \qquad \ensuremath{\boldsymbol{e}}_i =
   J\left(\ensuremath{\boldsymbol{e}}^j \times \ensuremath{\boldsymbol{e}}^k\right) \qquad i,j,k \texttt{ cyc } 1,2,3\end{aligned}

Derivation of operators in the BOUT++ Clebsch system
====================================================

The Clebsch system in BOUT++ goes like this

.. math::

   \begin{aligned}
       \ensuremath{\boldsymbol{B}}=&\nabla z \times \nabla x\\ =&\ensuremath{\boldsymbol{e}}^z \times \ensuremath{\boldsymbol{e}}^x\\
       J^{-1}\ensuremath{\boldsymbol{e}}_y=&\ensuremath{\boldsymbol{e}}^z \times \ensuremath{\boldsymbol{e}}^x\end{aligned}

 We have

.. math::

   \begin{aligned}
       B\ensuremath{\overset{\text{def}}{=}}& \sqrt{\ensuremath{\boldsymbol{B}}\cdot\ensuremath{\boldsymbol{B}}} = \sqrt{J^{-1}\ensuremath{\boldsymbol{e}}_y\cdot
   J^{-1}\ensuremath{\boldsymbol{e}}_y} = \sqrt{J^{-2}g_{yy}} = J^{-1}\sqrt{g_{yy}}\end{aligned}

 Further on

.. math::

   \begin{aligned}
       \ensuremath{\boldsymbol{B}}=&B\ensuremath{\boldsymbol{b}}\\ \ensuremath{\boldsymbol{b}}=&\frac{\ensuremath{\boldsymbol{B}}}{B}
       =\frac{J^{-1}\ensuremath{\boldsymbol{e}}_y}{J^{-1}\sqrt{g_{yy}}} =\frac{\ensuremath{\boldsymbol{e}}_y}{\sqrt{g_{yy}}}\end{aligned}

The parallel and perpendicular gradients
----------------------------------------

We have that

.. math::

   \begin{aligned}
       \ensuremath{\nabla}=& \ensuremath{\boldsymbol{e}}^i \partial_i = \ensuremath{\boldsymbol{e}}^x \partial_x + \ensuremath{\boldsymbol{e}}^y \partial_y +
       \ensuremath{\boldsymbol{e}}^z \partial_z\end{aligned}

 and that

.. math::

   \begin{aligned}
       \ensuremath{\nabla}_\| =& \left(\ensuremath{\boldsymbol{b}} \cdot \ensuremath{\nabla}\right) \ensuremath{\boldsymbol{b}} = \ensuremath{\boldsymbol{b}} \ensuremath{\boldsymbol{b}} \cdot \ensuremath{\nabla}=
       \frac{\ensuremath{\boldsymbol{e}}_y \ensuremath{\boldsymbol{e}}_y}{g_{yy}} \cdot \ensuremath{\nabla}= \frac{\ensuremath{\boldsymbol{e}}_y
       \ensuremath{\boldsymbol{e}}_y}{g_{yy}} \cdot \ensuremath{\boldsymbol{e}}^i \partial_i = \frac{\ensuremath{\boldsymbol{e}}_y}{g_{yy}}
       \partial_y\end{aligned}

 so that

.. math::

   \begin{aligned}
       \ensuremath{\nabla}_\perp =& \ensuremath{\nabla}- \ensuremath{\nabla}_\|\\
   %
                   =& \ensuremath{\boldsymbol{e}}^x \partial_x + \ensuremath{\boldsymbol{e}}^y \partial_y + \ensuremath{\boldsymbol{e}}^z
       \partial_z - \frac{\ensuremath{\boldsymbol{e}}_y}{g_{yy}} \partial_y\\
   %
                   =& \ensuremath{\boldsymbol{e}}^x \partial_x + \ensuremath{\boldsymbol{e}}^y \partial_y + \ensuremath{\boldsymbol{e}}^z
       \partial_z - \frac{g_{yi}\ensuremath{\boldsymbol{e}}^i}{g_{yy}} \partial_y\\
   %
                   =& \ensuremath{\boldsymbol{e}}^x \partial_x + \ensuremath{\boldsymbol{e}}^y \partial_y + \ensuremath{\boldsymbol{e}}^z
       \partial_z - \frac{g_{yx}\ensuremath{\boldsymbol{e}}^x +g_{yy}\ensuremath{\boldsymbol{e}}^y +g_{yz}\ensuremath{\boldsymbol{e}}^z
       }{g_{yy}}\partial_y\\
   %
                   =& \ensuremath{\boldsymbol{e}}^x \left(\partial_x - \frac{g_{yx}}{g_{yy}}\partial_y\right)
       +  \ensuremath{\boldsymbol{e}}^z \left(\partial_z - \frac{g_{yz}}{g_{yy}}\partial_y\right)\end{aligned}

The perpendicular gradients in Laplacian inversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the Laplacian inversion BOUT++ currently neglects the parallel
:math:`y` derivatives if :math:`g_{xy}` and :math:`g_{yz}` are non-zero,
thus

.. math::

   \begin{aligned}
       \ensuremath{\nabla}_\perp \simeq& \ensuremath{\boldsymbol{e}}^x \partial_x +  \ensuremath{\boldsymbol{e}}^z \partial_z
       \label{eq:reduced_grad_perp}\end{aligned}

The Laplacian
-------------

We would here like to find an expression for the Laplacian

.. math::

   \begin{aligned}
       \ensuremath{\nabla}^2 = \ensuremath{\nabla\cdot}\ensuremath{\nabla}\end{aligned}

 In general we have (using equation (2.6.39) in D’Haeseleer
:raw-latex:`\cite{haeseler-1}`)

.. math::

   \begin{aligned}
       \ensuremath{\nabla\cdot}\ensuremath{\boldsymbol{A}} = \frac{1}{J} \partial_i \left(JA^i\right)
       \label{eq:divA}\end{aligned}

 and that

.. math::

   \begin{aligned}
       A^i = \ensuremath{\boldsymbol{A}}\cdot \ensuremath{\boldsymbol{e}}^i\end{aligned}

 In our case :math:`A \to \ensuremath{\nabla}`, so that

.. math::

   \begin{aligned}
       \ensuremath{\nabla}^i = \left(\ensuremath{\nabla}\right)\cdot \ensuremath{\boldsymbol{e}}^i = \ensuremath{\boldsymbol{e}}^i \cdot \left(\ensuremath{\nabla}\right) = \ensuremath{\boldsymbol{e}}^i
       \cdot \left(\ensuremath{\boldsymbol{e}}^j \partial_j\right) = g^{ij} \partial_j\end{aligned}

 Thus

.. math::

   \begin{aligned}
       \ensuremath{\nabla}^2 =& \frac{1}{J} \partial_i \left(J g^{ij} \partial_j\right)\\ =&
       \frac{1}{J} g^{ij} J \partial_i \partial_j + \frac{1}{J} \partial_i \left(J
       g^{ij} \right) \partial_j\\ =& g^{ij} \partial_i \partial_j + G^j \partial_j\\\end{aligned}

 where we have defined  [1]_

.. math::

   \begin{aligned}
       G^j =& \frac{1}{J} \partial_i \left(J g^{ij} \right)\\ =& \frac{1}{J} \left(
       \partial_x \left[J g^{xj} \right] + \partial_y \left[J g^{yj} \right] + \partial_z \left[J
       g^{zj} \right] \right)\end{aligned}

 By writing the terms out, we get

.. math::

   \begin{aligned}
       \ensuremath{\nabla}^2 =& g^{ij} \partial_i \partial_j + G^j \partial_j\\
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
:math:`g^{ij}=g^{ji}`, and :math:`g_{ij}=g_{ji}`, and that the partial
derivatives commutes for smooth functions
:math:`\partial_i\partial_j=\partial_j\partial_i`. This gives

.. math::

   \begin{aligned}
       \ensuremath{\nabla}^2 =&\quad \, \left(g^{xx} \partial_x^2 \right) + \left(G^x \partial_x\right)\\ &+
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

 Notice that :math:`G^i` does not operate on :math:`\partial_i`, but
rather that the two are multiplied together.

The parallel Laplacian
----------------------

We have that

.. math::

   \begin{aligned}
       \ensuremath{\nabla}_\| =& \left(\ensuremath{\boldsymbol{b}} \cdot \ensuremath{\nabla}\right) \ensuremath{\boldsymbol{b}}\ = \ensuremath{\boldsymbol{b}} \ensuremath{\boldsymbol{b}} \cdot \ensuremath{\nabla}=
       \frac{\ensuremath{\boldsymbol{e}}_y \ensuremath{\boldsymbol{e}}_y}{g_{yy}} \cdot \ensuremath{\nabla}= \frac{\ensuremath{\boldsymbol{e}}_y
       \ensuremath{\boldsymbol{e}}_y}{g_{yy}} \cdot \ensuremath{\boldsymbol{e}}^i \partial_i = \frac{\ensuremath{\boldsymbol{e}}_y}{g_{yy}}
       \partial_y\end{aligned}

 we have that

.. math::

   \begin{aligned}
       \ensuremath{\nabla}_\|^i =& \left(\frac{\ensuremath{\boldsymbol{e}}_y}{g_{yy}} \partial_y\right)\cdot \ensuremath{\boldsymbol{e}}^i =
       \ensuremath{\boldsymbol{e}}^i \cdot \left(\frac{\ensuremath{\boldsymbol{e}}_y}{g_{yy}} \partial_y\right)\end{aligned}

 so that by equation (`[eq:divA] <#eq:divA>`__),

.. math::

   \begin{aligned}
       \ensuremath{\nabla}_\|^2 =& \ensuremath{\nabla\cdot}\left(\ensuremath{\boldsymbol{b}} \ensuremath{\boldsymbol{b}} \cdot \ensuremath{\nabla}\right)\\ =&
       \ensuremath{\nabla\cdot}\left(\frac{\ensuremath{\boldsymbol{e}}_y}{g_{yy}} \cdot \partial_y\right)\\ =& \frac{1}{J}
       \partial_i \left( J\ensuremath{\boldsymbol{e}}^i \cdot \left[\frac{\ensuremath{\boldsymbol{e}}_y}{g_{yy}} \partial_y\right]
       \right)\\ =& \frac{1}{J} \partial_y \left(\frac{J}{g_{yy}} \partial_y\right)\end{aligned}

The perpendicular Laplacian
---------------------------

For the perpendicular Laplacian, we have that

.. math::

   \begin{aligned}
       \ensuremath{\nabla}_\perp^2 =& \ensuremath{\nabla}^2 - \ensuremath{\nabla}_\|^2\\ =& g^{ij} \partial_i \partial_j +
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
dependent variable in Laplacian inversion if :math:`g_{xy}` and
:math:`g_{yz}` are non-zero (if these are zero, the derivation can be
done directly from equation
(`[eq:reduced_grad_perp] <#eq:reduced_grad_perp>`__) instead), so that

.. math::

   \begin{aligned}
       \ensuremath{\nabla}_\perp^2 \simeq& \quad \, \left(g^{xx} \partial_x^2\right) + \left( \frac{1}{J}
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
:math:`\ensuremath{\boldsymbol{v}}_E =
-\frac{\nabla\phi\times\ensuremath{\boldsymbol{b}}}{B}`, which is
similar to
:math:`\ensuremath{\boldsymbol{v}}=\ensuremath{\boldsymbol{k}}\times\nabla\psi`
found in incompressible fluid flow

.. math::

   \begin{aligned}
       \ensuremath{\boldsymbol{v}}_E =& -\frac{\nabla\phi\times\ensuremath{\boldsymbol{b}}}{B}\\
                %
                =&-\frac{\nabla\phi\times\ensuremath{\boldsymbol{e}}_y}{
   \sqrt{g_{yy}}J^{-1}\sqrt{g_{yy}}}\\
                %
                =&-\frac{J}{g_{yy}}\nabla\phi\times\ensuremath{\boldsymbol{e}}_y\\
                %
                =&\frac{J}{g_{yy}}\ensuremath{\boldsymbol{e}}_y\times\nabla\phi\\
                %
                =&\frac{J}{g_{yy}}\ensuremath{\boldsymbol{e}}_y\times \left(\ensuremath{\boldsymbol{e}}^x\partial_x + \ensuremath{\boldsymbol{e}}^y\partial_y +
   \ensuremath{\boldsymbol{e}}^z\partial_z\right)\phi\\
                %
                =&\frac{J}{g_{yy}} \left(g_{yx}\ensuremath{\boldsymbol{e}}^x + g_{yy}\ensuremath{\boldsymbol{e}}^y +
   g_{yz}\ensuremath{\boldsymbol{e}}^z\right) \times \left(\ensuremath{\boldsymbol{e}}^x\partial_x + \ensuremath{\boldsymbol{e}}^y\partial_y +
   \ensuremath{\boldsymbol{e}}^z\partial_z\right)\phi\\
                %
                =&\frac{J}{g_{yy}} \left( g_{yx}\ensuremath{\boldsymbol{e}}^x\times\ensuremath{\boldsymbol{e}}^x\partial_x +
   g_{yy}\ensuremath{\boldsymbol{e}}^y\times\ensuremath{\boldsymbol{e}}^x\partial_x + g_{yz}\ensuremath{\boldsymbol{e}}^z\times\ensuremath{\boldsymbol{e}}^x\partial_x
   \right.  \\ &\quad\; + g_{yx}\ensuremath{\boldsymbol{e}}^x\times\ensuremath{\boldsymbol{e}}^y\partial_y +
   g_{yy}\ensuremath{\boldsymbol{e}}^y\times\ensuremath{\boldsymbol{e}}^y\partial_y + g_{yz}\ensuremath{\boldsymbol{e}}^z\times\ensuremath{\boldsymbol{e}}^y\partial_y
   \\ &\quad\; \left.  + g_{yx}\ensuremath{\boldsymbol{e}}^x\times\ensuremath{\boldsymbol{e}}^z\partial_z +
   g_{yy}\ensuremath{\boldsymbol{e}}^y\times\ensuremath{\boldsymbol{e}}^z\partial_z + g_{yz}\ensuremath{\boldsymbol{e}}^z\times\ensuremath{\boldsymbol{e}}^z\partial_z
   \right) \phi\\
                %
                =&\frac{J}{g_{yy}} \left( - g_{yy}\ensuremath{\boldsymbol{e}}^y\times\ensuremath{\boldsymbol{e}}^x\partial_x +
   g_{yz}\ensuremath{\boldsymbol{e}}^z\times\ensuremath{\boldsymbol{e}}^x\partial_x \right.  \\ &\quad +
   g_{yx}\ensuremath{\boldsymbol{e}}^x\times\ensuremath{\boldsymbol{e}}^y\partial_y - g_{yz}\ensuremath{\boldsymbol{e}}^z\times\ensuremath{\boldsymbol{e}}^y\partial_y
   \\ &\quad \left.  - g_{yx}\ensuremath{\boldsymbol{e}}^x\times\ensuremath{\boldsymbol{e}}^z\partial_z +
   g_{yy}\ensuremath{\boldsymbol{e}}^y\times\ensuremath{\boldsymbol{e}}^z\partial_z \right) \phi\\
                %
                =&\frac{1}{g_{yy}} \left( - g_{yy}\ensuremath{\boldsymbol{e}}_z\partial_x +
   g_{yz}\ensuremath{\boldsymbol{e}}_y\partial_x + g_{yx}\ensuremath{\boldsymbol{e}}_z\partial_y - g_{yz}\ensuremath{\boldsymbol{e}}_x\partial_y
   - g_{yx}\ensuremath{\boldsymbol{e}}_y\partial_z + g_{yy}\ensuremath{\boldsymbol{e}}_x\partial_z \right) \phi\end{aligned}

The electrostatic ExB advection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The electrostatic :math:`E\times B` advection operator thus becomes

.. math::

   \begin{aligned}
       \ensuremath{\boldsymbol{v}}_E\cdot\nabla =& -\frac{\nabla\phi\times\ensuremath{\boldsymbol{b}}}{B}\cdot\nabla\\
       %
       =&\frac{1}{g_{yy}} \left( - g_{yy}\ensuremath{\boldsymbol{e}}_z\partial_x +
       g_{yz}\ensuremath{\boldsymbol{e}}_y\partial_x + g_{yx}\ensuremath{\boldsymbol{e}}_z\partial_y -
       g_{yz}\ensuremath{\boldsymbol{e}}_x\partial_y - g_{yx}\ensuremath{\boldsymbol{e}}_y\partial_z +
       g_{yy}\ensuremath{\boldsymbol{e}}_x\partial_z \right) \phi \cdot\left(\ensuremath{\boldsymbol{e}}^x\partial_x +
       \ensuremath{\boldsymbol{e}}^y\partial_y + \ensuremath{\boldsymbol{e}}^z\partial_z\right)\\
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

   \begin{aligned}
       B\ensuremath{\boldsymbol{v}}_E\cdot\nabla =& -\nabla\phi\times\ensuremath{\boldsymbol{b}}\cdot\nabla\\
       %
       =& \ensuremath{\boldsymbol{b}} \times \nabla\phi\cdot\nabla\\
       %
       =& \frac{\sqrt{g_{yy}}}{Jg_{yy}} \left( g_{yx}\{\phi, \cdot\}_{y,z} +
       g_{yy}\{\phi, \cdot\}_{z,x} + g_{yz}\{\phi, \cdot\}_{x,y} \right) \\
       %
       =& \frac{1}{J\sqrt{g_{yy}}} \left( g_{yx}\{\phi, \cdot\}_{y,z} + g_{yy}\{\phi,
   \cdot\}_{z,x} + g_{yz}\{\phi, \cdot\}_{x,y} \right) \addtocounter{equation}{1}\tag{\theequation}
                  \label{eq:brackets}\end{aligned}

The brackets operator in BOUT++
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Notice that the (phi,f)@ operators in BOUT++ returns
:math:`-\frac{\nabla\phi\times\ensuremath{\boldsymbol{b}}}{B}\cdot\nabla f`
rather than
:math:`-\nabla\phi\times\ensuremath{\boldsymbol{b}}\cdot\nabla f`.

Notice also that the Arakawa brackets neglects the :math:`\partial_y`
derivative terms (the :math:`y`-derivative terms) if :math:`g_{xy}` and
:math:`g_{yz}` are non-zero, so for the Arakawa brackets, BOUT++ returns

.. math::

   \begin{aligned}
       \ensuremath{\boldsymbol{v}}_E\cdot\nabla =& -\frac{\nabla\phi\times\ensuremath{\boldsymbol{b}}}{B}\cdot\nabla\\
       %
       \simeq& \frac{1}{g_{yy}} \left( g_{yy}\{\phi, \cdot\}_{z,x} \right)\\
       %
       =& \partial_z\phi\partial_x - \partial_x\phi\partial_z\end{aligned}

Divergence of ExB velocity
==========================

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{v}}_{ExB} = \frac{\ensuremath{\boldsymbol{b}}\times\nabla\phi}{B}\end{aligned}

 Using

.. math::

   \begin{aligned}
   \nabla\cdot\left(\ensuremath{\boldsymbol{F}}\times\ensuremath{\boldsymbol{G}}\right) = \left(\nabla\times\ensuremath{\boldsymbol{F}}\right)\cdot\ensuremath{\boldsymbol{G}} -
   \ensuremath{\boldsymbol{F}}\cdot\left(\nabla\times\ensuremath{\boldsymbol{G}}\right)\end{aligned}

 the divergence of the
:math:`\ensuremath{\boldsymbol{E}}\times\ensuremath{\boldsymbol{B}}`
velocity can be written as

.. math::

   \begin{aligned}
   \nabla\cdot\left(\frac{1}{B}\ensuremath{\boldsymbol{b}}\times\nabla\phi\right) =
   \left[\nabla\times\left(\frac{1}{B}\ensuremath{\boldsymbol{b}}\right)\right]\cdot\nabla\phi -
   \frac{1}{B}\ensuremath{\boldsymbol{b}}\cdot\nabla\times\nabla\phi
   \label{eq:exb1}\end{aligned}

 The second term on the right is identically zero (curl of a nablaient).
The first term on the right can be expanded as

.. math::

   \begin{aligned}
   \left[\nabla\times\left(\frac{1}{B}\ensuremath{\boldsymbol{b}}\right)\right]\cdot\nabla\phi =
   \left[\nabla\left(\frac{1}{B}\right)\times\ensuremath{\boldsymbol{b}} +
   \frac{1}{B}\nabla\times\ensuremath{\boldsymbol{b}}\right]\cdot\nabla\phi\end{aligned}

 Using

.. math::

   \begin{aligned}
   \ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}} = \nabla\times\ensuremath{\boldsymbol{b}} -
   \ensuremath{\boldsymbol{b}}\left[\ensuremath{\boldsymbol{b}}\cdot\left(\nabla\times\ensuremath{\boldsymbol{b}}\right)\right]\end{aligned}

 this becomes:

.. math::

   \begin{aligned}
     \nabla\cdot\left(\frac{1}{B}\ensuremath{\boldsymbol{b}}\times\nabla\phi\right) =
     &-\ensuremath{\boldsymbol{b}}\times\nabla\left(\frac{1}{B}\right)\cdot\nabla\phi \\ &+
     \frac{1}{B}\ensuremath{\boldsymbol{b}}\times\ensuremath{\boldsymbol{\kappa}}\cdot\nabla\phi \\ &+
     \left[\ensuremath{\boldsymbol{b}}\cdot\left(\nabla\times\ensuremath{\boldsymbol{b}}\right)\right]\ensuremath{\boldsymbol{b}}\cdot\nabla\phi\end{aligned}

 Alternatively, equation \ `[eq:exb1] <#eq:exb1>`__ can be expanded as

.. math::

   \begin{aligned}
     \nabla\cdot\left(\frac{1}{B}\ensuremath{\boldsymbol{b}}\times\nabla\phi\right) =&
       -B\ensuremath{\boldsymbol{b}}\times\nabla\left(\frac{1}{B^2}\right)\cdot\nabla\phi +
       \frac{1}{B^2}\nabla\times\ensuremath{\boldsymbol{B}}\cdot\nabla\phi \\ =&
       -B\ensuremath{\boldsymbol{b}}\times\nabla\left(\frac{1}{B^2}\right)\cdot\nabla\phi +
       \frac{1}{B^2}\ensuremath{\boldsymbol{J}}\cdot\nabla\phi\end{aligned}

.. math::

   \begin{aligned}
   \nabla\cdot\left(n\frac{\mathbf{b}\times\nabla\phi}{B}\right) &=& \frac{1}{J}\frac{\partial}{\partial\psi}\left(Jn\frac{\partial\phi}{\partial z} \right) - \frac{1}{J}\frac{\partial}{\partial z}\left(Jn\frac{\partial\phi}{\partial\psi}\right)  \\
                                                                 &+& \frac{1}{J}\frac{\partial}{\partial\psi}\left(Jn\frac{g^{\psi\psi}g^{yz}}{B^2}\frac{\partial\phi}{\partial y}\right) - \frac{1}{J}\frac{\partial}{\partial y}\left(Jn\frac{g^{\psi\psi}g^{yz}}{B^2}\frac{\partial\phi}{\partial\psi}\right)\end{aligned}

.. [1]
   | Notice that :math:`G^i` is **not** the same as the *Christoffel
     symbols of second kind* (also known as the *connection
     coefficients* or
     :math:`\Gamma^i_{jk}=\ensuremath{\boldsymbol{e}}^i\cdot\partial_k \ensuremath{\boldsymbol{e}}_j`),
     although the derivation of the two are quite similar.
   | We find that
     :math:`\Gamma^i_{ji}=\ensuremath{\boldsymbol{e}}^i\cdot\partial_i \ensuremath{\boldsymbol{e}}_j = \ensuremath{\nabla\cdot}\ensuremath{\boldsymbol{e}}_j`,
     whereas using equation `[eq:divA] <#eq:divA>`__ leads to
     :math:`G^i=\ensuremath{\boldsymbol{e}}^i\cdot\partial_i \ensuremath{\boldsymbol{e}}^j = \ensuremath{\nabla\cdot}
     \ensuremath{\boldsymbol{e}}^j`, since :math:`g^{ji}=g^{ij}` due to
     symmetry.
