Compressible gas dynamics
=========================

This example solves the Euler equations for compressible gas dynamics. The evolving
variables are the density "N", pressure "P", and three components of the velocity "V".
The equations are the continuity equation for density:

    dN/dt + Div(V * N) = 0

the momentum equation

    dV/dt + V*Grad(V) = -Grad(P)/N + g + nu * Laplace(V)

where g is the acceleration due to gravity, and nu is the viscosity.
Note that only a simple viscosity form is used here.
The pressure equation is

    dP/dt + Div(V * P) = -(gamma-1)*P*Div(V)

where gamma is the ratio of specific heats.

Test cases
----------

* advect1d is a simple advection in 1D periodic domain. This is just to test advection schemes

    $ ./gas_compress -d advect1d

* sod-shock solves the Sod shock tube problem, a standard 1D test case for fluid codes

    $ ./gas_compress -d sod-shock

* rayleigh-taylor is a 2D simulation, which starts with a dense fluid above a lighter fluid with
  gravity pointing downwards. This situation is unstable, so the initial perturbation grows.

    $ ./gas_compress -d rayleigh-taylor

