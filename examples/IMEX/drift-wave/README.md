Drift wave test case
====================

This example evolves electron density and vorticity, and includes
a resistive Ohm's law with parallel electron pressure.

    dn/dt + [phi, n] + Div_par(ve) = 0

    dU/dt + [phi, U] + Div_par(Ve) = 0

    nu * Ve = Grad_par(phi) - Grad_par(Ne)

These equations are solved in a slab geometry, starting
with a sinusoidal perturbation.

The value of resistivity nu is set in the BOUT.inp file
(or on the command line). By decreasing nu the parallel dynamics
can be made faster, and the problem more "stiff".

Running with CVODE, treating both convective and diffusive parts together:

    $ ./test-drift

    1.000e+02        574              574       1.87e+00    57.2   29.8    4.7    0.3    8.0
    2.000e+02        233              233       7.91e-01    57.2   29.5    4.9    0.6    7.7
    3.000e+02        237              237       7.73e-01    57.1   29.6    4.7    0.7    7.8

and with IMEXBDF2:

    $ ./test-drift solver:type=imexbdf2

    1.000e+02          2               67       1.71e-01    44.6   33.4    6.6    3.0   12.4
    2.000e+02          2               66       1.81e-01    46.3   32.0    6.7    2.8   12.3
    3.000e+02          2               63       9.37e-02    43.4   31.9    5.4    6.0   13.2



