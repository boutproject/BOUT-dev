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
    
    1.000e+02        467              467       5.40e-01    50.9   38.5    2.1    0.5    8.0
    2.000e+02        159              159       1.91e-01    49.2   37.3    2.0    1.4   10.2
    3.000e+02         99               99       1.22e-01    47.8   36.1    2.0    2.1   12.0

and with IMEXBDF2 (adaptive timestep):

    $ ./test-drift solver:type=imexbdf2 solver:maxl=50
    
    1.000e+02          2                7       1.55e-02    20.0   13.2    0.9   16.1   49.8
    2.000e+02          2               18       2.36e-02    31.7   19.1    1.2   10.7   37.3
    3.000e+02          2               25       2.88e-02    35.5   21.1    1.3    8.6   33.4


Things to try:

1. Turning off adaptive timestepping: solver:adaptive=false
2. Modifying maximum linear iterations (solver:maxl) and nonlinear iterations (max_nonlinear_it)




