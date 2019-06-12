.. default-role:: math

.. _sec-nonlocal-heatflux:


Nonlocal heat flux models
=========================

Spitzer-Harm heat flux
----------------------

The Spitzer-Harm heat flux `q_{SH}` is calculated using

.. math::

   q_{SH} = - \frac{n_e e T_e}{m_e}\frac{3\sqrt{\pi}}{4}\tau_{ei,T}\kappa_0\frac{Z+0.24}{Z+4.2} \partial_{||} T_e

where `n_e` is the electron density in `m^{-3}`, `T_e` is the electron temperature in eV, `kappa_0 = 13.58`,
`Z` is the average ion charge. The resulting expression is in units of `eV/m^2/s`. 

The thermal collision time `tau_{ei,T} = \lambda_{ei,T} / v_{T}` is calculated using the thermal mean free path
and thermal velocity:

.. math::

   \lambda_{ee,T} = \frac{v_T^3}{Yn_e \ln\Lambda}
   
   \lambda_{ei,T} = \frac{v_T^3}{YZ^2n_i \ln\Lambda}
   
   v_T = \sqrt{\frac{2eT_e}{m_e}}

where it is assumed that `n_i = n_e`, and the following are used:

.. math::

   Y = 4\pi\left(\frac{e^2}{4\pi \epsilon_0 m_e}\right)^2

   \ln\Lambda = 6.6 - 0.5\log\left(\frac{n_e}{10^{20}}\right) + 1.5 \log\left(T_e\right);

SNB model
---------
   
The SNB model calculates a correction to the Spitzer-Harm heat flux, solving a
diffusion equation for each of a set of energy groups with normalised
energy `\beta = E_g / eT_e` where `E_g` is the energy of the group.
   
.. math::

   \left[\frac{1}{\lambda'_{g,ee}} - \nabla_{||}\left(\frac{\lambda'_{g,ei}}{3}\partial_{||}\right)\right]H_g = -\nabla_{||} U_g


where `\nabla_{||}` is the divergence of a parallel flux, and `\partial_{||}` is a parallel gradient.
`U_g = W_g q_{SH}` is the contribution to the Spitzer-Harm heat flux from a group:

.. math::

   W_g = \frac{1}{24}\int_{\beta_{g-1}}^{\beta^{g+1}} \beta^4 e^{-\beta} d\beta

The modified mean free paths for each group are:

.. math::

   \lambda'_{g,ee} = \beta^2 \lambda_{ee,T} / r

   \lambda'_{g,ei} = \beta^2 \lambda_{ei,T} \frac{Z + 0.24}{Z + 4.2}

From the quantities `H_g` for each group, the SNB heat flux is:

.. math::

   q_{SNB} = q_{SH} - \sum_g\frac{\lambda_g,ei}{3}\nabla H_g

In flud models we actually want the divergence of the heat flux, rather than the heat flux itself.
We therefore rearrange to get:

.. math::

   \nabla_{||}\left(\frac{\lambda'_{g,ei}}{3}\partial_{||}\right)H_g = \nabla_{||} U_g + H_g / \lambda'_{g,ee}

and so calculate the divergence of the heat flux as:

.. math::

   \nabla_{||} q_{SNB} = \nabla_{||} q_{SH} - \sum_g\left(\nabla_{||} U_g + H_g / \lambda'_{g,ee}\right)


The Helmholtz type equation along the magnetic field is solved using a tridiagonal solver.
The parallel divergence term is currently split into a second derivative term, and a first derivative correction:

.. math::

   \nabla_{||}\left(k\partial_{||} T\right) = \frac{1}{J}\frac{\partial}{\partial y}\left(\frac{k J}{g_{22}}\frac{\partial T}{\partial y}\right)
   = k\frac{1}{g_22}\frac{\partial^2 T}{\partial y^2} + \frac{1}{J}\frac{\partial}{\partial y}\left(\frac{k J}{g_{22}}\right)\frac{\partial T}{\partial y}

   
