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
the `FV::Div_par` function in `bout/fv_ops.hxx`. This uses the MC
slope limiter, together with a Lax flux at the local sound speed
to provide dissipation and minimise unphysical oscillations.

See also:
https://bout-dev.readthedocs.io/en/latest/user_docs/differential_operators.html#parallel-divergence-div-par

MMS test
--------

To run an MMS convergence test, use the Python script:

    $ ./runtest

This will use the analytic solution and sources calculated in `mms.py`,
and given in the `mms/BOUT.inp` inputs. The return code indicates success (0)
or failure (1). If matplotlib is installed, then this should also output a
figure "fluid_norm.pdf" and "fluid_norm.png".

The default slope limiter used by `FV::Div_par` is Monotonised Central (MC).
To test other limiters, replace

    FV::Div_par

with

    FV::Div_par<FV::Upwind>

for the first order upwinding method, or

    FV::Div_par<FV::MinMod>

for the second order MinMod method.



