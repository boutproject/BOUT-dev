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

    1.000e+02        556              556       1.81e+00    62.0   29.4    3.3    0.3    5.0
    2.000e+02        253              253       8.30e-01    61.8   29.3    3.2    0.8    4.9
    3.000e+02        366              366       1.17e+00    62.1   29.3    3.3    0.5    4.9

and with IMEXBDF2 (adaptive timestep):

    $ ./test-drift solver:type=imexbdf2 solver:maxl=50


    1.000e+02          5              310       1.28e+00    39.9   14.9    1.9    0.5   42.8
    2.000e+02         11             1690       4.43e+00    42.3   15.7    2.1    0.1   39.7
    3.000e+02          8              272       6.92e-01    42.4   16.0    2.1    0.9   38.6

Increasing the maximum number of linear iterations to 100:

    $ ./test-drift solver:type=imexbdf2 solver:maxl=100

    1.000e+02          2               64       2.54e-01    41.8   15.6    2.1    2.5   38.1
    2.000e+02          2               57       2.24e-01    42.3   15.9    2.1    2.8   36.9
    3.000e+02          2               53       2.08e-01    42.3   15.9    2.1    3.1   36.6
   
Without adaptive timestep:

    $ ./test-drift solver:type=imexbdf2 solver:maxl=100 solver:adaptive=false

    1.000e+02          2               64       2.53e-01    42.0   15.7    2.0    2.7   37.6
    2.000e+02          2               57       2.23e-01    42.4   15.9    2.0    3.1   36.5
    3.000e+02          2               53       2.08e-01    42.4   15.9    2.1    3.3   36.3



