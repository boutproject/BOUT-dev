Nonlinear diffusion test
========================

Solves a diffusion equation in which the coefficient
depends on the evolving variable. Simplified example
of heat conduction in a plasma.

Running the code without any options:

    $ ./diffusion-nl

will use the default time integration solver: CVODE if available, otherwise PVODE
and produce an output similar to:

    1.000e+00        145              145       1.84e-02    48.9    0.0    3.9   12.3   34.9
    2.000e+00         64               64       1.18e-02    40.0    0.0    2.3   18.9   38.8
    3.000e+00         49               49       9.31e-03    33.1    0.0    2.5   23.2   41.2

Use the IMEX-BDF2 multistep scheme by running:

    $ ./diffusion-nl solver:type=imexbdf2

which should produce something like

    1.000e+00          2              868       2.33e-02    35.8    0.0    9.6    4.6   49.9
    2.000e+00          2              339       1.36e-02    30.3    0.0    9.2    7.6   52.9
    3.000e+00          2              279       1.08e-02    30.0    0.0    7.7   10.1   52.2


Preconditioning
---------------

The preconditioner can be enabled by adding another flag (or setting in BOUT.inp)

    $ ./diffusion-nl solver:use_precon=true

which then greatly reduces the number of iterations needed by CVODE:

    1.000e+00         88               88       9.34e-03    31.3    0.0    1.9   16.5   50.3
    2.000e+00         22               22       4.83e-03    16.5    0.0    1.6   31.4   50.5
    3.000e+00         16               16       4.02e-03    14.5    0.0    1.2   34.6   49.6

and for IMEX-BDF2:

    $ ./diffusion-nl solver:type=imexbdf2 solver:use_precon=true

    1.000e+00          2               99       7.21e-03    14.3    0.0    2.9   14.1   68.7
    2.000e+00          2               71       5.43e-03    13.4    0.0    3.8   21.9   60.8
    3.000e+00          2               58       5.30e-03    13.2    0.0    3.1   21.5   62.2

Jacobian coloring
-----------------

By default the IMEX-BDF2 solver uses a matrix-free method, but the Jacobian
can be calculated using finite differences. This uses coloring to improve efficiency

    $ ./diffusion-nl solver:type=imexbdf2 solver:matrix_free=false

    1.000e+00          2               26       3.10e-03    31.9    0.0    5.9    7.8   54.4
    2.000e+00          2               23       2.24e-03    35.7    0.0    7.0    7.9   49.4
    3.000e+00          2               18       2.04e-03    35.0    0.0    6.2    8.5   50.4

