Nonlinear diffusion test
========================

Solves a diffusion equation in which the coefficient
depends on the evolving variable. Simplified example
of heat conduction in a plasma.

Running the code without any options:

    $ ./diffusion-nl

will use the default time integration solver: CVODE if available, otherwise PVODE
and produce an output similar to:

    1.000e+00        138              138       1.04e-02    21.7    0.0    1.0   19.2   58.1
    2.000e+00         67               67       8.78e-03    12.4    0.0    0.6   23.3   63.8
    3.000e+00         56               56       8.50e-03    11.1    0.0    0.5   23.4   65.1

Use the IMEX-BDF2 multistep scheme by running:

    $ ./diffusion-nl solver:type=imexbdf2

which should produce something like

    1.000e+00         22             1718       3.06e-01     5.5    0.0    0.7    0.7   93.1
    2.000e+00         24             1190       9.72e-02    11.9    0.0    1.5    2.1   84.5
    3.000e+00         17             1007       1.19e-01     8.3    0.0    1.1    1.7   88.9

The settings can be adjusted, so to see what the solver is doing internally add "solver:verbose=true".
This shows that the KSP (linear) and SNES solvers are not converging. This is probably because they
are reaching the set maximum number of iterations. To increase the number of allowed iterations:

    $ ./diffusion-nl solver:type=imexbdf2 solver:maxl=100 solver:max_nonlinear_it=20

should produce:

    1.000e+00          2              579       5.28e-02    10.1    0.0    1.2    3.8   84.8
    2.000e+00          2              215       2.39e-02     8.5    0.0    1.0    8.4   82.1
    3.000e+00          2              170       2.05e-02     7.8    0.0    0.9    9.9   81.4

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

    1.000e+00          9              452       1.63e-01     3.0    0.0    0.4    1.3   95.4
    2.000e+00          8              287       3.57e-02     8.2    0.0    1.0    5.6   85.2
    3.000e+00          6              169       2.39e-02     7.2    0.0    0.9    8.4   83.4

Jacobian coloring
-----------------

By default the IMEX-BDF2 solver uses a matrix-free method, but the Jacobian
can be calculated using finite differences. This uses coloring to improve efficiency

    $ ./diffusion-nl solver:type=imexbdf2 solver:matrix_free=false

    1.000e+00          2              137       1.85e-02    10.9    0.0    1.0   10.7   77.4
    2.000e+00          2               72       1.45e-02     7.4    0.0    0.7   14.6   77.3
    3.000e+00          2               71       1.40e-02     8.3    0.0    0.8   14.1   76.8


