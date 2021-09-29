test-beuler
===========

Integrate a stiff system:

    ddt(f) = 998 * f + 1998 * (g - 1.0);
    ddt(g) = -999 * f - 1999 * (g - 1.0);

starting with f=1, g=0. The solution has an exp(-t) term, and
stiff exp(-1000t) term which can be challenging to integrate.
The solution should converge to f=0, g=1.

This is an example of a problem where many time integration
solvers will fail (including CVODE with standard settings),
but are quite easily solved with Backward Euler.
