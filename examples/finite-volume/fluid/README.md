Fluid model in 1D
=================

This example solves a 1D set of adiabatic fluid equations,
evolving the density (n), pressure (p) and parallel momentum (nv).

Continuity (density) equation:

    dn/dt + Div(n v) = 0

Pressure equation:

    dp/dt + Div(p v) = - (gamma-1) * p * Div(v)

Momentum equation:

    d(nv)/dt + Div(nv v) = -Grad(p)

The advection terms (divergences on the left) are solved using
the FV::Div_par function in `bout/fv_ops.hxx`. This uses the MC
slope limiter, together with a Lax flux at the local sound speed
to provide dissipation and minimise unphysical oscillations. 

