IMEX-BDF2
=========

This is a second-order multistep IMplicit-EXplicit (IMEX) scheme
taken from this paper:

W.Hundsdorfer, S.J.Ruuth "IMEX extensions of linear multistep methods with general monotonicity and boundedness properties" JCP 225 (2007) 2016-2042

which is available at: http://homepages.cwi.nl/~willem/DOCART/JCP07.pdf

Each time step a single nonlinear solve is required, which is done
using PETSc's SNES routines.

Options
-------

The following options are used:

* timestep     = The internal time step
* matrix_free  = True/false (default True). Determines whether the Jacobian in SNES is matrix free
* use_coloring = True/false (default True). If not matrix free, use coloring to calculate Jacobian?
* lag_jacobian = Integer number of times to (re-)use Jacobian. Default is 4
* atol         = Absolute tolerance (1e-16)
* rtol         = Relative tolerance (1e-10)
* predictor    = Predictor method (default 1)
* use_precon   = True/false (default False). Use a user-defined preconditioner?

Predictor
---------

The implicit solve is an iterative method (typically GMRES) through SNES,
and works best if given a good starting guess for the solution. This is
done by extrapolating from previous times. The predictor option
chooses between:

0. Constant: Next step is the same as the last step
1. Linear: Extrapolate from the last two steps
2. Quadratic: Extrapolate from the last three steps
3. Explicit only: Start with the result of the explicit step

Linear extrapolation (default, 1) has been found to work well generally.

