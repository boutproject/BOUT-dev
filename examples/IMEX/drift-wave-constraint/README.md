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

This version of the code solves for the potential as a constraint.
This makes the Jacobian of the RHS function sparse, since it no
longer contains inverse Laplacian operations.

Some care must be taken to ensure that the time derivatives
only depend on nearby grid points, since these are all that are
currently included in the coloring.

Using the IMEX-BDF2 solver:

    $ ./test-drift solver:type=imexbdf2 solver:matrix_free=false

    1.000e+02          2               35       1.25e-01    64.3    0.0    3.7    0.9   31.1
    2.000e+02          2               35       1.21e-01    66.3    0.0    3.8    0.9   29.0
    3.000e+02          2               35       1.22e-01    65.8    0.0    3.7    0.9   29.6

    $ ./test-drift solver:type=imexbdf2 solver:matrix_free=false -pc_type sor

    1.000e+02          2               35       1.66e-01    48.6    0.0    2.8    0.6   47.9
    2.000e+02          2               35       1.49e-01    54.2    0.0    3.1    0.7   42.0
    3.000e+02          2               35       1.62e-01    50.1    0.0    2.9    0.7   46.4

This is comparable to the time taken for the drift-wave test with matrix-free IMEX-BDF2,
though takes a smaller number of RHS evaluations



