.. _sec-time-integration:

Time integration
================

.. _sec-timeoptions:

Options
-------

BOUT++ can be compiled with several different time-integration solvers ,
and at minimum should have Runge-Kutta (RK4) and PVODE (BDF/Adams)
solvers available.

The solver library used is set using the ``solver:type`` option, so
either in BOUT.inp:

.. code-block:: cfg

    [solver]
    type = rk4  # Set the solver to use

or on the command line by adding ``solver:type=pvode`` for example:

.. code-block:: bash

    mpirun -np 4 ./2fluid solver:type=rk4

**NB**: Make sure there are no spaces around the “=” sign:
``solver:type =pvode`` won’t work (probably). Table :numref:`tab-solvers` gives
a list of time integration solvers, along with any compile-time options
needed to make the solver available.

.. _tab-solvers:
.. table:: Available time integration solvers
	   
   +---------------+-----------------------------------------+------------------------+
   | Name          | Description                             | Compile options        |
   +===============+=========================================+========================+
   | euler         | Euler explicit method (example only)    | Always available       |
   +---------------+-----------------------------------------+------------------------+
   | rk4           | Runge-Kutta 4th-order explicit method   | Always available       |
   +---------------+-----------------------------------------+------------------------+
   | rkgeneric     | Generic Runge Kutta explicit methods    | Always available       |
   +---------------+-----------------------------------------+------------------------+
   | rk3ssp        | 3rd-order Strong Stability Preserving   | Always available       |
   +---------------+-----------------------------------------+------------------------+
   | splitrk       | Split RK3-SSP and RK-Legendre           | Always available       |
   +---------------+-----------------------------------------+------------------------+
   | pvode         | 1998 PVODE with BDF method              | Always available       |
   +---------------+-----------------------------------------+------------------------+
   | cvode         | SUNDIALS CVODE. BDF and Adams methods   | -DBOUT_USE_SUNDIALS=ON |
   +---------------+-----------------------------------------+------------------------+
   | ida           | SUNDIALS IDA. DAE solver                | -DBOUT_USE_SUNDIALS=ON |
   +---------------+-----------------------------------------+------------------------+
   | arkode        | SUNDIALS ARKODE IMEX solver             | -DBOUT_USE_SUNDIALS=ON |
   +---------------+-----------------------------------------+------------------------+
   | petsc         | PETSc TS methods                        | -DBOUT_USE_PETSC=ON    |
   +---------------+-----------------------------------------+------------------------+
   | imexbdf2      | IMEX-BDF2 scheme                        | -DBOUT_USE_PETSC=ON    |
   +---------------+-----------------------------------------+------------------------+
   | beuler / snes | Backward Euler with SNES solvers        | -DBOUT_USE_PETSC=ON    |
   +---------------+-----------------------------------------+------------------------+

Each solver can have its own settings which work in slightly different
ways, but some common settings and which solvers they are used in are
given in table :numref:`tab-solveropts`.

.. _tab-solveropts:
.. table:: Time integration solver options
	   
   +--------------------------+--------------------------------------------+-------------------------------------+
   | Option                   | Description                                | Solvers used                        |
   +==========================+============================================+=====================================+
   | atol                     | Absolute tolerance                         | rk4, pvode, cvode, ida, imexbdf2,   |
   |                          |                                            | beuler                              |
   +--------------------------+--------------------------------------------+-------------------------------------+
   | rtol                     | Relative tolerance                         | rk4, pvode, cvode, ida, imexbdf2,   |
   |                          |                                            | beuler                              |
   +--------------------------+--------------------------------------------+-------------------------------------+
   | mxstep                   | Maximum internal steps                     | rk4, imexbdf2                       |
   |                          | per output step                            |                                     |
   +--------------------------+--------------------------------------------+-------------------------------------+
   | max\_timestep            | Maximum timestep                           | rk4, cvode                          |
   +--------------------------+--------------------------------------------+-------------------------------------+
   | timestep                 | Starting timestep                          | rk4, euler, imexbdf2, beuler        |
   +--------------------------+--------------------------------------------+-------------------------------------+
   | adaptive                 | Adapt timestep? (Y/N)                      | rk4, imexbdf2                       |
   +--------------------------+--------------------------------------------+-------------------------------------+
   | use\_precon              | Use a preconditioner? (Y/N)                | pvode, cvode, ida, imexbdf2         |
   +--------------------------+--------------------------------------------+-------------------------------------+
   | mudq, mldq               | BBD preconditioner settings                | pvode, cvode, ida                   |
   +--------------------------+--------------------------------------------+-------------------------------------+
   | mukeep, mlkeep           |                                            |                                     |
   +--------------------------+--------------------------------------------+-------------------------------------+
   | maxl                     | Maximum number of linear iterations        | cvode, imexbdf2                     |
   +--------------------------+--------------------------------------------+-------------------------------------+
   | max_nonlinear_iterations | Maximum number of nonlinear iterations     | cvode, imexbdf2, beuler             |
   +--------------------------+--------------------------------------------+-------------------------------------+
   | use\_jacobian            | Use user-supplied Jacobian? (Y/N)          | cvode                               |
   +--------------------------+--------------------------------------------+-------------------------------------+
   | adams\_moulton           | Use Adams-Moulton method                   | cvode                               |
   |                          | rather than BDF                            |                                     |
   +--------------------------+--------------------------------------------+-------------------------------------+
   | diagnose                 | Collect and print additional diagnostics   | cvode, imexbdf2, beuler             |
   +--------------------------+--------------------------------------------+-------------------------------------+

|

The most commonly changed options are the absolute and relative solver
tolerances, ``atol`` and ``rtol`` which should be varied to check
convergence.

CVODE
-----

The most commonly used time integration solver is CVODE, or its older
version PVODE. CVODE has several advantages over PVODE, including better
support for preconditioning and diagnostics.

Enabling diagnostics output using ``solver:diagnose=true`` will print a
set of outputs for each timestep similar to:

.. code-block:: bash

    CVODE: nsteps 51, nfevals 69, nniters 65, npevals 126, nliters 79
        -> Newton iterations per step: 1.274510e+00
        -> Linear iterations per Newton iteration: 1.215385e+00
        -> Preconditioner evaluations per Newton: 1.938462e+00
        -> Last step size: 1.026792e+00, order: 5
        -> Local error fails: 0, nonlinear convergence fails: 0
        -> Stability limit order reductions: 0
    1.000e+01        149       2.07e+01    78.3    0.0   10.0    0.9   10.8

When diagnosing slow performance, key quantities to look for are
nonlinear convergence failures, and the number of linear iterations per
Newton iteration. A large number of failures, and close to 5 linear
iterations per Newton iteration are a sign that the linear solver is not
converging quickly enough, and hitting the default limit of 5
iterations. This limit can be modified using the ``solver:maxl``
setting. Giving it a large value e.g. ``solver:maxl=1000`` will show how
many iterations are needed to solve the linear system. If the number of
iterations becomes large, this may be an indication that the system is
poorly conditioned, and a preconditioner might help improve performance.
See :ref:`sec-preconditioning`.

CVODE can set constraints to keep some quantities positive, non-negative,
negative or non-positive. These constraints can be activated by setting the
option ``solver:apply_positivity_constraints=true``, and then in the section
for a certain variable (e.g. ``[n]``), setting the option
``positivity_constraint`` to one of ``positive``, ``non_negative``,
``negative``, or ``non_positive``.

Additional options can be used to modify the behaviour of the linear and
nonlinear solvers:

- ``cvode_nonlinear_convergence_coef`` specifies the safety factor
  used in the nonlinear convergence test. Passed as a parameter to
  `CVodeSetNonlinConvCoef
  <https://sundials.readthedocs.io/en/latest/cvodes/Usage/SIM.html#c.CVodeSetNonlinConvCoef>`_.

- ``cvode_linear_convergence_coef`` specifies the factor by which the
  Krylov linear solver’s convergence test constant is reduced from the
  nonlinear solver test constant. Passed as a parameter to
  `CVodeSetEpsLin
  <https://sundials.readthedocs.io/en/latest/cvodes/Usage/SIM.html#c.CVodeSetEpsLin>`_.

The linear solver type can be set using the ``linear_solver`` option.
Valid choices include ``gmres`` (the default), ``fgmres``, ``tfqmr``, ``bcgs``.

IMEX-BDF2
---------

This is an IMplicit-EXplicit time integration solver, which allows the
evolving function to be split into two parts: one which has relatively
long timescales and can be integrated using explicit methods, and a
part which has short timescales and must be integrated implicitly. The
order of accuracy is variable (up to 4th-order currently), and an
adaptive timestep can be used.

To use the IMEX-BDF2 solver, set the solver type to ``imexbdf2``,
e.g. on the command-line add ``solver:type=imexbdf2`` or in the
options file:

.. code-block:: cfg

    [solver]
    type = imexbdf2


The order of the method is set to 2 by default, but can be increased up to a maximum of 4:

.. code-block:: cfg

    [solver]
    type = imexbdf2
    maxOrder = 3

This is a multistep method, so the state from previous steps are used
to construct the next one. This means that at the start, when there
are no previous steps, the order is limited to 1 (backwards Euler
method). Similarly, the second step is limited to order 2, and so
on. At the moment the order is not adapted, so just increases until
reaching `maxOrder`.

At each step the explicit (non-stiff) part of the function is called,
and combined with previous timestep values. The implicit part of the
function is then solved using PETSc's SNES, which consists of a
nonlinear solver (usually modified Newton iteration), each iteration
of which requires a linear solve (usually GMRES). Settings which
affect this implicit part of the solve are:

+------------------+-----------+----------------------------------------------------+
| Option           | Default   |Description                                         |
+==================+===========+====================================================+
| atol             | 1e-16     | Absolute tolerance on SNES solver                  |
+------------------+-----------+----------------------------------------------------+
| rtol             | 1e-10     | Relative tolerance on SNES solver                  |
+------------------+-----------+----------------------------------------------------+
| max_nonlinear_it | 5         | Maximum number of nonlinear iterations             |
|                  |           | If adaptive timestepping is used then              |
|                  |           | failure will cause timestep reduction              |
+------------------+-----------+----------------------------------------------------+
| maxl             | 20        | Maximum number of linear iterations                |
|                  |           | If adaptive, failure will cause timestep reduction |
+------------------+-----------+----------------------------------------------------+
| predictor        | 1         | Starting guess for the nonlinear solve             |
|                  |           | Specifies order of extrapolating polynomial        |
+------------------+-----------+----------------------------------------------------+
| use_precon       | false     | Use user-supplied preconditioner?                  |
+------------------+-----------+----------------------------------------------------+
| matrix_free      | true      | Use Jacobian-free methods? If false, calculates    |
|                  |           | the Jacobian matrix using finite difference        |
+------------------+-----------+----------------------------------------------------+
| use_coloring     | true      | If not matrix free, use coloring to speed up       |
|                  |           | calculation of the Jacobian                        |
+------------------+-----------+----------------------------------------------------+


Note that the SNES tolerances `atol` and `rtol` are set very conservatively by default. More reasonable
values might be 1e-10 and 1e-5, but this must be explicitly asked for in the input options.

The predictor extrapolates from previous timesteps to get a starting estimate for the value
at the next timestep. This estimate is then used to initialise the SNES nonlinear solve.
The value is the order of the extrapolating polynomial, so 1 (the default) is a linear extrapolation
from the last two steps, 0 is the same as the last step. A value of -1 uses the explicit
update to the state as the starting guess, i.e. assuming that the implicit part of the problem is small.
This is usually not a good guess.

To diagnose what is happening in the time integration, for example to see why it is
failing to converge or why timesteps are small, there are two settings which can be
set to ``true`` to enable:

- `diagnose` outputs a summary at each output time, similar to CVODE. This
  contains information like the last timestep, average number of iterations
  and number of convergence failures.
- `verbose` prints information at every internal step, with more information
  on the values used to modify timesteps, and the reasons for solver failures.

By default adaptive timestepping is turned on, using several factors to
modify the timestep:

#. If the nonlinear solver (SNES) fails to converge, either because it diverges or exceeds the iteration limits
   `max_nonlinear_its` or `maxl`. Reduces the timestep by 2 and tries again, giving up after 10 failures.

#. Every `nadapt` internal timesteps (default 4), the error is checked by taking the timestep twice:
   Once with the current order of accuracy, and once with one order of accuracy lower. The difference
   between the solutions is then used to estimate the timestep required to achieve the required
   tolerances. If this is much larger or smaller than the current timestep, then the timestep is modified.

#. The timestep is kept within user-specified maximum and minimum ranges.


The options which control this behaviour are:

+------------------+-----------+----------------------------------------------------+
| Option           | Default   |Description                                         |
+==================+===========+====================================================+
| adaptive         | true      | Turns on adaptive timestepping                     |
+------------------+-----------+----------------------------------------------------+
| timestep         | output    | If adaptive sets the starting timestep.            |
|                  | timestep  | If not adaptive, timestep fixed at this value      |
+------------------+-----------+----------------------------------------------------+
| dtMin            | 1e-10     | Minimum timestep                                   |
+------------------+-----------+----------------------------------------------------+
| dtMax            | output    | Maximum timestep                                   |
|                  | timestep  |                                                    |
+------------------+-----------+----------------------------------------------------+
| mxstep           | 1e5       | Maximum number of internal steps between outputs   |
+------------------+-----------+----------------------------------------------------+
| nadapt           | 4         | How often is error checked and timestep adjusted?  |
+------------------+-----------+----------------------------------------------------+
| adaptRtol        | 1e-3      | Target relative tolerance for adaptive timestep    |
+------------------+-----------+----------------------------------------------------+
| scaleCushDown    | 1.0       | Timestep scale factor below which the timestep is  |
|                  |           | modified. By default the timestep is always reduced|
+------------------+-----------+----------------------------------------------------+
| scaleCushUp      | 1.5       | Minimum timestep scale factor based on adaptRtol   |
|                  |           | above which the timestep will be modified.         |
|                  |           | Currently the timestep increase is limited to 25%  |
+------------------+-----------+----------------------------------------------------+


Split-RK
--------

The `splitrk` solver type uses Strang splitting to combine two
explicit Runge Kutta schemes:

#. `2nd order Runge-Kutta-Legendre method <https://doi.org/10.1016/j.jcp.2013.08.021>`_
   for the diffusion (parabolic) part. These schemes use
   multiple stages to increase stability, rather than accuracy; this
   is always 2nd order, but the stable timestep for diffusion
   problems increases as the square of the number of stages. The
   number of stages is an input option, and can be arbitrarily large.

#. 3rd order SSP-RK3 scheme for the advection (hyperbolic) part
   http://www.cscamm.umd.edu/tadmor/pub/linear-stability/Gottlieb-Shu-Tadmor.SIREV-01.pdf

Each timestep consists of

#. A half timestep of the diffusion part
#. A full timestep of the advection part
#. A half timestep of the diffusion part

Options to control the behaviour of the solver are:

+------------------+-----------+----------------------------------------------------+
| Option           | Default   |Description                                         |
+==================+===========+====================================================+
| timestep         | output    | If adaptive sets the starting timestep.            |
|                  | timestep  | If not adaptive, timestep fixed at this value      |
+------------------+-----------+----------------------------------------------------+
| nstages          | 10        | Number of stages in RKL step. Must be > 1          |
+------------------+-----------+----------------------------------------------------+
| diagnose         | false     |  Print diagnostic information                      |
+------------------+-----------+----------------------------------------------------+

And the adaptive timestepping options:

+---------------------+-----------+----------------------------------------------------+
| Option              | Default   |Description                                         |
+=====================+===========+====================================================+
| adaptive            | true      | Turn on adaptive timestepping                      |
+---------------------+-----------+----------------------------------------------------+
| atol                | 1e-10     | Absolute tolerance                                 |
+---------------------+-----------+----------------------------------------------------+
| rtol                | 1e-5      | Relative tolerance                                 |
+---------------------+-----------+----------------------------------------------------+
| max_timestep        | output    | Maximum internal timestep                          |
|                     | timestep  |                                                    |
+---------------------+-----------+----------------------------------------------------+
| max_timestep_change | 2         | Maximum factor by which the timestep by which the  |
|                     |           | time step can be changed at each step              |
+---------------------+-----------+----------------------------------------------------+
| mxstep              | 1000      | Maximum number of internal steps before output     |
+---------------------+-----------+----------------------------------------------------+
| adapt_period        | 1         | Number of internal steps between tolerance checks  |
+---------------------+-----------+----------------------------------------------------+

Backward Euler - SNES
---------------------

The `beuler` or `snes` solver type (either name can be used) is a PETSc-based implicit
solver for finding steady-state solutions to systems of partial differential equations.
It supports multiple solution strategies including backward Euler timestepping,
direct Newton iteration, and Pseudo-Transient Continuation (PTC) with Switched
Evolution Relaxation (SER).

Basic Configuration
~~~~~~~~~~~~~~~~~~~

The SNES solver is configured through the ``[solver]`` section of the input file:

.. code-block:: ini

   [solver]
   type = snes

   # Nonlinear solver settings
   snes_type = newtonls          # anderson, newtonls, newtontr, nrichardson
   atol = 1e-7                   # Absolute tolerance
   rtol = 1e-6                   # Relative tolerance
   stol = 1e-12                  # Solution change tolerance
   max_nonlinear_iterations = 20 # Maximum SNES iterations per solve

   # Linear solver settings
   ksp_type = fgmres             # Linear solver: gmres, bicgstab, etc.
   maxl = 20                     # Maximum linear iterations
   pc_type = ilu                 # Preconditioner: ilu, bjacobi, hypre, etc.

Timestepping Modes
~~~~~~~~~~~~~~~~~~

The solver supports several timestepping strategies controlled by ``equation_form``:

**Backward Euler (default)**
   Standard implicit backward Euler method. Good for general timestepping.

   .. code-block:: ini

      equation_form = rearranged_backward_euler  # Default

   This method has low accuracy in time but its dissipative properties
   are helpful when evolving to steady state solutions.

**Direct Newton**
   Solves the steady-state problem F(u) = 0 directly without timestepping.

   .. code-block:: ini

      equation_form = direct_newton

   This method is unlikely to converge unless the system is very close
   to steady state.

**Pseudo-Transient Continuation**
   Uses pseudo-time to guide the solution to steady state. Recommended for
   highly nonlinear problems where Newton's method fails.

   .. code-block:: ini

      equation_form = pseudo_transient

   This uses the same form as rearranged_backward_euler, but the time step
   can be different for each cell.

Adaptive Timestepping
~~~~~~~~~~~~~~~~~~~~~

When ``equation_form = rearranged_backward_euler`` (default), the
solver uses global timestepping with adaptive timestep control based
on nonlinear iteration count.

.. code-block:: ini

   [solver]
   type = snes
   equation_form = rearranged_backward_euler

   # Initial and maximum timesteps
   timestep = 1.0                         # Initial timestep
   max_timestep = 1e10                    # Upper limit on timestep
   dt_min_reset = 1e-6                    # Reset the solver when timestep < this

   # Timestep adaptation
   lower_its = 3                          # Increase dt if iterations < this
   upper_its = 10                         # Decrease dt if iterations > this
   timestep_factor_on_lower_its = 1.4     # Growth factor
   timestep_factor_on_upper_its = 0.9     # Reduction factor
   timestep_factor_on_failure = 0.5       # Reduction on convergence failure

PID Controller
^^^^^^^^^^^^^^

An alternative adaptive strategy using a PID controller:

.. code-block:: ini

   [solver]
   pid_controller = true
   target_its = 7        # Target number of nonlinear iterations
   kP = 0.7              # Proportional gain
   kI = 0.3              # Integral gain
   kD = 0.2              # Derivative gain

The PID controller adjusts the timestep to maintain approximately ``target_its``
nonlinear iterations per solve, providing smoother adaptation than threshold-based
methods.

Pseudo-Transient Continuation and Switched Evolution Relaxation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``equation_form = pseudo_transient`` the solver uses
Pseudo-Transient Continuation (PTC). This is a robust numerical
technique for solving steady-state problems that are too nonlinear for
direct Newton iteration. Instead of solving the steady-state system
**F(u) = 0** directly, PTC solves a modified time-dependent problem:

.. math::

   M(u) \frac{\partial u}{\partial \tau} + F(u) = 0

where :math:`\tau` is a pseudo-time variable (not physical time) and :math:`M(u)`
is a preconditioning matrix. As :math:`\tau \to \infty`, the solution converges
to the steady state **F(u) = 0**.

The key advantage of PTC is that it transforms a difficult root-finding problem
into a sequence of easier initial value problems. Poor initial guesses that would
cause Newton's method to diverge can still reach the solution via a stable
pseudo-transient path.

The Switched Evolution Relaxation (SER) method is a spatially adaptive
variant of PTC that allows each cell to use a different
pseudo-timestep :math:`\Delta\tau_i`. The timestep in each cell adapts
based on the local residual, allowing the algorithm to take large
timesteps in well-behaved regions (fast convergence), while taking
small timesteps in difficult regions (stable advancement).  The the
same :math:`\Delta\tau_i` is used for all equations (density,
momentum, energy etc.) within each cell. This maintains coupling
between temperature, pressure, and composition through the equation of
state.

**Key parameters:**

``pseudo_max_ratio`` (default: 2.0)
   Maximum allowed ratio of timesteps between neighboring cells. This prevents
   sharp spatial gradients in convergence rate.

**Example PTC configuration:**

.. code-block:: ini

   [solver]
   type = snes
   equation_form = pseudo_transient

   timestep = 1.0                # Initial timestep

   # SER parameters
   pid_controller = true          # Scale timesteps based on iterations
   pseudo_max_ratio = 2.0         # Limit neighbor timestep ratio

   # Tolerances
   atol = 1e-7
   rtol = 1e-6
   stol = 1e-12

SER timestep strategy
^^^^^^^^^^^^^^^^^^^^^

After each nonlinear solve the timesteps in each cell are adjusted.
The strategy used depends on the ``pseudo_strategy`` option:

**inverse_residual** (default)

If ``pseudo_strategy = inverse_residual`` then the timestep is inversely
proportional to the RMS residual in each cell.
``pseudo_alpha`` (default: 100 × atol × timestep)
Controls the relationship between residual and timestep. The local timestep
is computed as:

.. math::

   \Delta\tau_i = \frac{\alpha}{||R_i||}

Larger values allow more aggressive timestepping. The default is to use
a fixed ``pseudo_alpha`` but a better strategy is to enable the PID controller
that adjusts this parameter based on the nonlinear solver convergence.

The timestep is limited to be between ``dt_min_reset`` and
``max_timestep``.  In addition the timestep is limited between 0.67 ×
previous timestep and 1.5 × previous timestep, to limit sudden changes
in timestep.

In practice this strategy seems to work well, though problems could
arise when residuals become very small.

**history_based**

When ``pseudo_strategy = history_based`` the history of residuals
within each cell is used to adjust the timestep. The key parameters
are:

``pseudo_growth_factor`` (default: 1.1)
   Factor by which timestep increases when residual decreases successfully.

``pseudo_reduction_factor`` (default: 0.5)
   Factor by which timestep decreases when residual increases (step rejected).

This method may be less susceptible to fluctuations when residuals
become small, but tends to be slower to converge when residuals are
large.

**hybrid**

When ``pseudo_strategy = hybrid`` the ``inverse_residual`` and
``history_based`` strategies are combined: When the residuals are
large the ``inverse_residual`` method is used, and when residuals
become small the method switches to ``history_based``.

PID Controller
^^^^^^^^^^^^^^

When using the PTC method the PID controller can be used to dynamically
adjust ``pseudo_alpha`` depending on the nonlinearity of the system:

.. code-block:: ini

   [solver]
   pid_controller = true
   target_its = 7        # Target number of nonlinear iterations
   kP = 0.7              # Proportional gain
   kI = 0.3              # Integral gain
   kD = 0.2              # Derivative gain

The PID controller adjusts ``pseudo_alpha``, scaling all cell
timesteps together, to maintain approximately ``target_its`` nonlinear
iterations per solve.

With this enabled the solver uses the number of nonlinear iterations
to scale timesteps globally, and residuals to scale timesteps locally.
Note that the PID controller has no effect on the ``history_based``
strategy because that strategy does not use ``pseudo_alpha``.

Jacobian Finite Difference with Coloring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The default and recommended approach for most problems:

.. code-block:: ini

   [solver]
   use_coloring = true               # Enable (default)
   lag_jacobian = 5                  # Reuse Jacobian for this many iterations

   # Stencil shape (determines Jacobian sparsity pattern)
   stencil:taxi = 2                  # Taxi-cab distance (default)
   stencil:square = 0                # Square stencil extent
   stencil:cross = 0                 # Cross stencil extent

The coloring algorithm exploits the sparse structure of the Jacobian to reduce
the number of function evaluations needed for finite differencing.

Jacobian coloring stencil
^^^^^^^^^^^^^^^^^^^^^^^^^

The stencil used to create the Jacobian colouring can be varied,
depending on which numerical operators are in use. It is important to
note that the coloring won't work for every problem: It assumes that
each evolving quantity is coupled to all other evolving quantities on
the same grid cell, and on all the neighbouring grid cells. If the RHS
function includes Fourier transforms, or matrix inversions
(e.g. potential solves) then these will introduce longer-range
coupling and the Jacobian calculation will give spurious
results. Generally the method will then fail to converge. Two
solutions are to a) switch to matrix-free (``matrix_free=true``),
or b) solve the matrix inversion as a constraint.


``solver:stencil:cross = N``
e.g. for N == 2

.. code-block:: bash

        *
        *
    * * x * *
        *
        *


``solver:stencil:square = N``
e.g. for N == 2

.. code-block:: bash

    * * * * *
    * * * * *
    * * x * *
    * * * * *
    * * * * *

``solver:stencil:taxi = N``
e.g. for N == 2

.. code-block:: bash

        *
      * * *
    * * x * *
      * * *
        *

Setting ``solver:force_symmetric_coloring = true``, will make sure
that the jacobian colouring matrix is symmetric.  This will often
include a few extra non-zeros that the stencil will miss otherwise

Diagnostics and Monitoring
---------------------------

.. code-block:: ini

   [solver]
   diagnose = true                # Print iteration info to screen
   diagnose_failures = true       # Detailed diagnostics on failures

When ``equation_form = pseudo_transient``, the solver saves additional diagnostic fields:

- ``snes_pseudo_residual``: Local residual in each cell
- ``snes_pseudo_timestep``: Local pseudo-timestep in each cell
- ``snes_pseudo_alpha``: Global timestep scaling

These can be visualized to understand convergence behavior and identify
problematic regions.

Summary of solver options
~~~~~~~~~~~~~~~~~~~~~~~~~

+---------------------------+---------------+----------------------------------------------------+
| Option                    | Default       |Description                                         |
+===========================+===============+====================================================+
| pseudo_time               | false         | Pseudo-Transient Continuation (PTC) method, using  |
|                           |               | a different timestep for each cell.                |
+---------------------------+---------------+----------------------------------------------------+
| pseudo_max_ratio          | 2.            | Maximum timestep ratio between neighboring cells   |
+---------------------------+---------------+----------------------------------------------------+
| snes_type                 | newtonls      | PETSc SNES nonlinear solver (try anderson, qn)     |
+---------------------------+---------------+----------------------------------------------------+
| ksp_type                  | gmres         | PETSc KSP linear solver                            |
+---------------------------+---------------+----------------------------------------------------+
| pc_type                   | ilu / bjacobi | PETSc PC preconditioner (try hypre in parallel)    |
+---------------------------+---------------+----------------------------------------------------+
| pc_hypre_type             | pilut         | If ``pc_type = hypre``.                            |
|                           |               | Hypre preconditioner type: euclid, boomeramg       |
+---------------------------+---------------+----------------------------------------------------+
| max_nonlinear_iterations  | 20            | If exceeded, solve restarts with timestep / 2      |
+---------------------------+---------------+----------------------------------------------------+
| maxl                      | 20            | Maximum number of linear iterations                |
+---------------------------+---------------+----------------------------------------------------+
| atol                      | 1e-12         | Absolute tolerance of SNES solve                   |
+---------------------------+---------------+----------------------------------------------------+
| rtol                      | 1e-5          | Relative tolerance of SNES solve                   |
+---------------------------+---------------+----------------------------------------------------+
| upper_its                 | 80% max       | If exceeded, next timestep reduced by 10%          |
+---------------------------+---------------+----------------------------------------------------+
| lower_its                 | 50% max       | If under this, next timestep increased by 10%      |
+---------------------------+---------------+----------------------------------------------------+
| timestep                  | 1             | Initial timestep                                   |
+---------------------------+---------------+----------------------------------------------------+
| predictor                 | true          | Use linear predictor?                              |
+---------------------------+---------------+----------------------------------------------------+
| matrix_free               | false         | Matrix-free preconditioning?                       |
+---------------------------+---------------+----------------------------------------------------+
| matrix_free_operator      | false         | Use matrix free Jacobian-vector product?           |
+---------------------------+---------------+----------------------------------------------------+
| use_coloring              | true          | If ``matrix_free=false``, use coloring to speed up |
|                           |               | calculation of the Jacobian elements.              |
+---------------------------+---------------+----------------------------------------------------+
| lag_jacobian              | 50            | Re-use the Jacobian for successive inner solves    |
+---------------------------+---------------+----------------------------------------------------+
| kspsetinitialguessnonzero | false         | If true, Use previous solution as KSP initial      |
+---------------------------+---------------+----------------------------------------------------+
| use_precon                | false         | If ``matrix_free=true``, use user-supplied         |
|                           |               | preconditioner?                                    |
|                           |               | If false, the default PETSc preconditioner is used |
+---------------------------+---------------+----------------------------------------------------+
| diagnose                  | false         | Print diagnostic information every iteration       |
+---------------------------+---------------+----------------------------------------------------+
| stencil:cross             | 0             | If ``matrix_free=false`` and ``use_coloring=true`` |
| stencil:square            | 0             | Set the size and shape of the Jacobian coloring    |
| stencil:taxi              | 2             | stencil.                                           |
+---------------------------+---------------+----------------------------------------------------+
| force_symmetric_coloring  | false         | Ensure that the Jacobian coloring is symmetric     |
+---------------------------+---------------+----------------------------------------------------+

The predictor is linear extrapolation from the last two timesteps. It seems to be
effective, but can be disabled by setting ``predictor = false``.

The default `newtonls` SNES type can be very effective if combined
with Jacobian coloring: The coloring enables the Jacobian to be
calculated relatively efficiently; once a Jacobian matrix has been
calculated, effective preconditioners can be used to speed up
convergence.

The `SNES type
<https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESType.html>`_
can be set through PETSc command-line options, or in the BOUT++
options as setting `snes_type`. Good choices for unpreconditioned
problems where the Jacobian is not available (``matrix_free=true``) seem to be `anderson
<https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESANDERSON.html#SNESANDERSON>`_
and `qn
<https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESQN.html#SNESQN>`_
(quasinewton).

Preconditioner types:

#. On one processor the ILU solver is typically very effective, and is usually the default
#. The Hypre package can be installed with PETSc and used as a preconditioner. One of the
   options available in Hypre is the Euler parallel ILU solver.
   Enable with command-line args ``-pc_type hypre -pc_hypre_type euclid -pc_hypre_euclid_levels k``
   where ``k`` is the level (1-8 typically).



ODE integration
---------------

The `Solver` class can be used to solve systems of ODEs inside a physics
model: Multiple Solver objects can exist besides the main one used for
time integration. Example code is in ``examples/test-integrate``.

To use this feature, systems of ODEs must be represented by a class
derived from `PhysicsModel`.

::

    class MyFunction : public PhysicsModel {
     public:
      int init(bool restarting) {
        // Initialise ODE
        // Add variables to solver as usual
        solver->add(result, "result");
        ...
      }

      int rhs(BoutReal time) {
        // Specify derivatives of fields as usual
        ddt(result) = ...
      }
     private:
      Field3D result;
    };

To solve this ODE, create a new `Solver` object::

    Solver* ode = Solver::create(Options::getRoot()->getSection("ode"));

This will look in the section ``[ode]`` in the options file.
**Important:** To prevent this solver overwriting the main restart files
with its own restart files, either disable restart files:

.. code-block:: cfg

    [ode]
    enablerestart = false

or specify a different directory to put the restart files:

.. code-block:: cfg

    [ode]
    restartdir = ode  # Restart files ode/BOUT.restart.0.nc, ...

Create a model object, and pass it to the solver::

    MyFunction* model = new MyFunction();
    ode->setModel(model);

Finally tell the solver to perform the integration::

    ode->solve(5, 0.1);

The first argument is the number of steps to take, and the second is the
size of each step. These can also be specified in the options, so
calling

::

    ode->solve();

will cause ode to look in the input for ``nout`` and ``timestep``
options:

.. code-block:: cfg

    [ode]
    nout = 5
    timestep = 0.1

Finally, delete the model and solver when finished::

    delete model;
    delete solver;

**Note:** If an ODE needs to be solved multiple times, at the moment it
is recommended to delete the solver, and create a new one each time.

.. _sec-preconditioning:

Preconditioning
---------------

At every time step, an implicit scheme such as BDF has to solve a
non-linear problem to find the next solution. This is usually done using
Newton’s method, each step of which involves solving a linear (matrix)
problem. For :math:`N` evolving variables is an :math:`N\times N` matrix
and so can be very large. By default matrix-free methods are used, in
which the Jacobian :math:`\mathcal{J}` is approximated by finite
differences (see next subsection), and so this matrix never needs to be
explicitly calculated. Finding a solution to this matrix can still be
difficult, particularly as :math:`\delta t` gets large compared with
some time-scales in the system (i.e. a stiff problem).

A preconditioner is a function which quickly finds an approximate
solution to this matrix, speeding up convergence to a solution. A
preconditioner does not need to include all the terms in the problem
being solved, as the preconditioner only affects the convergence rate
and not the final solution. A good preconditioner can therefore
concentrate on solving the parts of the problem with the fastest
time-scales.

A simple example  [1]_ is a coupled wave equation, solved in the
``test-precon`` example code:

.. math::

   \frac{\partial u}{\partial t} = \partial_{||}v \qquad \frac{\partial
   v}{\partial t} = \partial_{||} u

First, calculate the Jacobian of this set of equations by taking
partial derivatives of the time-derivatives with respect to each of the
evolving variables

.. math::

   \mathcal{J} = (\begin{array}{cc}
   \frac{\partial}{\partial u}\frac{\partial u}{\partial t} &
   \frac{\partial}{\partial v}\frac{\partial u}{\partial t}\\
   \frac{\partial}{\partial u}\frac{\partial v}{\partial t} &
   \frac{\partial}{\partial v}\frac{\partial v}{\partial t}
   \end{array}
   ) = (\begin{array}{cc}
   0 & \partial_{||} \\
   \partial_{||} & 0
   \end{array}
   )

In this case :math:`\frac{\partial u}{\partial t}` doesn’t depend on
:math:`u` nor :math:`\frac{\partial v}{\partial t}` on :math:`v`, so the
diagonal is empty. Since the equations are linear, the Jacobian doesn’t
depend on :math:`u` or :math:`v` and so

.. math::

   \frac{\partial}{\partial t}(\begin{array}{c} u \\
   v \end{array}) = \mathcal{J} (\begin{array}{c} u \\
   v \end{array} )

In general for non-linear functions :math:`\mathcal{J}` gives the
change in time-derivatives in response to changes in the state variables
:math:`u` and :math:`v`.

In implicit time stepping, the preconditioner needs to solve an equation

.. math::

   \mathcal{I} - \gamma \mathcal{J}

where :math:`\mathcal{I}` is the identity matrix, and :math:`\gamma`
depends on the time step and method (e.g. :math:`\gamma = \delta t` for
backwards Euler method). For the simple wave equation problem, this is

.. math::

   \mathcal{I} - \gamma \mathcal{J} = (\begin{array}{cc}
   1 & -\gamma\partial_{||} \\
   -\gamma\partial_{||} & 1
   \end{array}
   )

This matrix can be block inverted using Schur factorisation  [2]_

.. math::

   (\begin{array}{cc}
     {\mathbf{E}} & {\mathbf{U}} \\
     {\mathbf{L}} & {\mathbf{D}}
   \end{array})^{-1}
    = (\begin{array}{cc}
     {\mathbf{I}} & -{\mathbf{E}}^{-1}{\mathbf{U}} \\
     0 & {\mathbf{I}}
   \end{array}
   )(\begin{array}{cc}
     {\mathbf{E}}^{-1} & 0 \\
     0 & {\mathbf{P}}_{Schur}^{-1}
   \end{array}
   )(\begin{array}{cc}
     {\mathbf{I}} & 0 \\
     -{\mathbf{L}}{\mathbf{E}}^{-1} & {\mathbf{I}}
   \end{array}
   )

where
:math:`{\mathbf{P}}_{Schur} = {\mathbf{D}} - {\mathbf{L}}{\mathbf{E}}^{-1}{\mathbf{U}}`
Using this, the wave problem becomes:

.. math::
   :label: precon

   (\begin{array}{cc} 1 & -\gamma\partial_{||} \\
   -\gamma\partial_{||} & 1 \end{array})^{-1} = (\begin{array}{cc} 1 & \gamma\partial_{||}\\
   0 & 1 \end{array} )(\begin{array}{cc} 1 & 0 \\
   0 & (1 -\gamma^2\partial^2_{||})^{-1} \end{array} )(\begin{array}{cc} 1 & 0\\
   \gamma\partial_{||} & 1 \end{array} )

The preconditioner is implemented by defining a function of the form

::

    int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
      ...
    }

which takes as input the current time, the :math:`\gamma` factor
appearing above, and :math:`\delta` which is only important for
constrained problems (not discussed here... yet). The current state of
the system is stored in the state variables (here ``u`` and ``v`` ),
whilst the vector to be preconditioned is stored in the time derivatives
(here ``ddt(u)`` and ``ddt(v)`` ). At the end of the preconditioner the
result should be in the time derivatives. A preconditioner which is just
the identity matrix and so does nothing is therefore::

    int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
    }

To implement the preconditioner in equation :eq:`precon`, first apply the
rightmost matrix to the given vector:

.. math::

   (\begin{array}{c}
   \texttt{ddt(u)} \\
   \texttt{ddt(v)}
   \end{array}
   ) = (\begin{array}{cc}
   1 & 0 \\
   \gamma\partial_{||} & 1
   \end{array}
   )(\begin{array}{c}
   \texttt{ddt(u)} \\
   \texttt{ddt(v)}
   \end{array}
   )

::

    int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
      mesh->communicate(ddt(u));
      //ddt(u) = ddt(u);
      ddt(v) = gamma*Grad_par(ddt(u)) + ddt(v);

note that since the preconditioner is linear, it doesn’t depend on
:math:`u` or :math:`v`. As in the RHS function, since we are taking a
differential of ``ddt(u)``, it first needs to be communicated to
exchange guard cell values.

The second matrix

.. math::

   (\begin{array}{c}
   \texttt{ddt(u)} \\
   \texttt{ddt(v)}
   \end{array}
   ) \rightarrow (\begin{array}{cc}
   1 & 0 \\
   0 & (1 - \gamma^2\partial^2_{||})^{-1}
   \end{array}
   )(\begin{array}{c}
   \texttt{ddt(u)} \\
   \texttt{ddt(v)}
   \end{array}
   )

doesn’t alter :math:`u`, but solves a parabolic equation in the
parallel direction. There is a solver class to do this called
`InvertPar` which solves the equation :math:`(A + B\partial_{||}^2)x =
b` where :math:`A` and :math:`B` are `Field2D` or constants [3]_. In
`PhysicsModel::init` we create one of these solvers::

    InvertPar *inv; // Parallel inversion class
    int init(bool restarting) {
       ...
       inv = InvertPar::Create();
       inv->setCoefA(1.0);
       ...
    }

In the preconditioner we then use this solver to update :math:`v`::

      inv->setCoefB(-SQ(gamma));
      ddt(v) = inv->solve(ddt(v));

which solves
:math:`ddt(v) \rightarrow (1 - \gamma^2\partial_{||}^2)^{-1} ddt(v)`.
The final matrix just updates :math:`u` using this new solution for
:math:`v`

.. math::

   (\begin{array}{c}
   \texttt{ddt(u)} \\
   \texttt{ddt(v)}
   \end{array}
   ) \rightarrow (\begin{array}{cc}
   1 & \gamma\partial_{||} \\
   0 & 1
   \end{array}
   )(\begin{array}{c}
   \texttt{ddt(u)} \\
   \texttt{ddt(v)}
   \end{array}
   )

::

      mesh->communicate(ddt(v));
      ddt(u) = ddt(u) + gamma*Grad_par(ddt(v));

Finally, boundary conditions need to be imposed, which should be
consistent with the conditions used in the RHS::

      ddt(u).applyBoundary("dirichlet");
      ddt(v).applyBoundary("dirichlet");

To use the preconditioner, pass the function to the solver in
`PhysicsModel::init`::

    int init(bool restarting) {
      solver->setPrecon(precon);
      ...
    }

then in the ``BOUT.inp`` settings file switch on the preconditioner

.. code-block:: bash

    [solver]
    type = cvode          # Need CVODE or PETSc
    use_precon = true     # Use preconditioner
    rightprec = false     # Use Right preconditioner (default left)

Jacobian function
-----------------

DAE constraint equations
------------------------

Using the IDA or IMEX-BDF2 solvers, BOUT++ can solve Differential
Algebraic Equations (DAEs), in which algebraic constraints are used for
some variables. Examples of how this is used are in the
``examples/constraints`` subdirectory.

First the variable to be constrained is added to the solver, in a
similar way to time integrated variables. For example

::

    Field3D phi;
    ...
    solver->constraint(phi, ddt(phi), "phi");

The first argument is the variable to be solved for (constrained). The
second argument is the field to contain the residual (error). In this
example the time derivative field ``ddt(phi)`` is used, but it could
be another `Field3D` variable. The solver will attempt to
find a solution to the first argument (``phi`` here) such that the
second argument (``ddt(phi)``) is zero to within tolerances.

In the RHS function the residual should be calculated. In this example
(``examples/constraints/drift-wave-constraint``) we have::

    ddt(phi) = Delp2(phi) - Vort;

so the time integration solver includes the algebraic constraint
``Delp2(phi) = Vort`` i.e. (:math:`\nabla_\perp^2\phi = \omega`).

IMEX-BDF2
---------

This is an implicit-explicit multistep method, which uses the PETSc
library for the SNES nonlinear solver. To use this solver, BOUT++ must
have been configured with PETSc support, and the solver type set to
``imexbdf2``

::

    [solver]
    type = imexbdf2

For examples of using IMEX-BDF2, see the ``examples/IMEX/``
subdirectory, in particular the ``diffusion-nl``, ``drift-wave`` and
``drift-wave-constrain`` examples.

The time step is currently fixed (not adaptive), and defaults to the
output timestep. To set a smaller internal timestep, the
``solver:timestep`` option can be set. If the timestep is too large,
then the explicit part of the problem may become unstable, or the
implicit part may fail to converge.

The implicit part of the problem can be solved matrix-free, in which
case the Jacobian-vector product is approximated using finite
differences. This is currently the default, and can be set on the
command-line using the options::

     solver:matrix_free=true  -snes_mf

Note the ``-snes_mf`` flag which is passed to PETSc. When using a matrix
free solver, the Jacobian is not calculated and so the amount of memory
used is minimal. However, since the Jacobian is not known, many standard
preconditioning methods cannot be used, and so in many cases a custom
preconditioner is needed to obtain good convergence.

An experimental feature uses PETSc’s ability to calculate the Jacobian
using finite differences. This can then speed up the linear solve, and
allows more options for preconditioning. To enable this option::

     solver:matrix_free=false

There are two ways to calculate the Jacobian: A brute force method which
is set up by this call to PETSc which is generally very slow, and a
“coloring” scheme which can be quite fast and is the default. Coloring
uses knowledge of where the non-zero values are in the Jacobian, to work
out which rows can be calculated simultaneously. The coloring code in
IMEX-BDF2 currently assumes that every field is coupled to every other
field in a star pattern: one cell on each side, a 7 point stencil for 3D
fields. If this is not the case for your problem, then the solver may
not converge.

The brute force method can be useful for comparing the Jacobian
structure, so to turn off coloring::

     solver:use_coloring=false

Using MatView calls, or the ``-mat_view`` PETSc options, the non-zero
structure of the Jacobian can be plotted or printed.

Monitoring the simulation output
--------------------------------

Monitoring of the solution can be done at two levels: output monitoring,
and timestep monitoring. Output monitoring occurs only when data is
written to file, whereas timestep monitoring is every timestep and so
(usually) much more frequent. Examples of both are in
``examples/monitor`` and ``examples/monitor-newapi``.

**Output monitoring**: At every output timestep the solver calls a
monitor method of the BoutMonitor class, which writes the output dump file,
calculates and prints timing information and estimated time remaining. If you
want to run additional code or write data to a different file, you can
implement the outputMonitor method of PhysicsModel::

    int outputMonitor(BoutReal simtime, int iter, int nout)

The first input is the current simulation time, the second is the output
number, and the last is the total number of outputs requested.
This method is called by a monitor object PhysicsModel::modelMonitor, which
writes the restart files at the same time. You can change the frequency at which
the monitor is called by calling, in PhysicsModel::init::

    modelMonitor.setTimestep(new_timestep)

where ``new_timestep`` is a BoutReal which is either ``timestep*n`` or
``timestep/n`` for an integer ``n``. Note that this will change the frequency
of writing restarts as well as of calling ``outputMonitor()``.

You can also add custom monitor object(s) for more flexibility.

You can call your output monitor class whatever you like, but it must be a
subclass of Monitor and provide the method ``call`` which takes 4 inputs and
returns an int::

    class MyOutputMonitor : public Monitor {
      int call(Solver *solver, BoutReal simtime, int iter, int NOUT) {
        ...
      }
    };

The first input is the solver object, the second is the current
simulation time, the third is the output number, and the last is the
total number of outputs requested. To get the solver to call this
function every output time, define a `MyOutputMonitor` object as a member of your
PhysicsModel::

      MyOutputMonitor my_output_monitor;

and put in your `PhysicsModel::init` code::

      solver->addMonitor(&my_output_monitor);

Note that the solver only stores a pointer to the `Monitor`, so you must make sure
the object is persistent, e.g. a member of a `PhysicsModel` class, not a local
variable in a constructor. If you want to later remove a monitor, you can do so with::

      solver->removeMonitor(&my_output_monitor);

A simple example using this monitor is::

    class MyOutputMonitor: public Monitor{
    public:
      MyOutputMonitor(BoutReal timestep=-1):Monitor(timestep){};
      int call(Solver *solver, BoutReal simtime, int iter, int NOUT) override;
    };

    int MyOutputMonitor::call(Solver *solver, BoutReal simtime, int iter, int NOUT) {
      output.write("Output monitor, time = %e, step %d of %d\n",
                   simtime, iter, NOUT);
      return 0;
    }

    MyOutputMonitor my_monitor;

    int init(bool restarting) {
      solver->addMonitor(&my_monitor);
    }

See the monitor example (``examples/monitor``) for full code.

**Timestep monitoring**: This uses functions instead of objects. First define a
monitor function::

    int my_timestep_monitor(Solver *solver, BoutReal simtime, BoutReal lastdt) {
      ...
    }

where ``simtime`` will again contain the current simulation time, and
``lastdt`` the last timestep taken. Add this function to the solver::

      solver->addTimestepMonitor(my_timestep_monitor);

Timestep monitoring is disabled by default, unlike output monitoring. To
enable timestep monitoring, set in the options file (BOUT.inp)::

    [solver]
    monitor_timestep = true

or put on the command line ``solver:monitor_timestep=true`` . When this
is enabled, it will change how solvers like CVODE and PVODE (the default
solvers) are used. Rather than being run in NORMAL mode, they will
instead be run in SINGLE\_STEP mode (see the SUNDIALS notes
here:\ https://computation.llnl.gov/casc/sundials/support/notes.html).
This may in some cases be less efficient.


Implementation internals
------------------------

The solver is the interface between BOUT++ and the time-integration
code such as SUNDIALS. All solvers implement the `Solver`
class interface (see ``src/solver/generic_solver.hxx``).

First all the fields which are to be evolved need to be added to the
solver. These are always done in pairs, the first specifying the field,
and the second the time-derivative::

    void add(Field2D &v, Field2D &F_v, const char* name);

This is normally called in the `PhysicsModel::init` initialisation routine.
Some solvers (e.g. IDA) can support constraints, which need to be added
in the same way as evolving fields::

    bool constraints();
    void constraint(Field2D &v, Field2D &C_v, const char* name);

The ``constraints()`` function tests whether or not the current solver
supports constraints. The format of ``constraint(...)`` is the same as
``add``, except that now the solver will attempt to make ``C_v`` zero.
If ``constraint`` is called when the solver doesn’t support them then an
error should occur.

If the physics model implements a preconditioner or Jacobian-vector
multiplication routine, these can be passed to the solver during
initialisation::

    typedef int (*PhysicsPrecon)(BoutReal t, BoutReal gamma, BoutReal delta);
    void setPrecon(PhysicsPrecon f); // Specify a preconditioner
    typedef int (*Jacobian)(BoutReal t);
    void setJacobian(Jacobian j); // Specify a Jacobian

If the solver doesn’t support these functions then the calls will just
be ignored.

Once the problem to be solved has been specified, the solver can be
initialised using::

    int init();

which returns an error code (0 on success). This is currently called in
:doc:`bout++.cxx<../_breathe_autogen/file/bout_09_09_8cxx>`::

    if (solver.init()) {
      output.write("Failed to initialise solver. Aborting\n");
      return(1);
    }

which passes the (physics module) RHS function `PhysicsModel::rhs` to the
solver along with the number and size of the output steps.

::

    typedef int (*MonitorFunc)(BoutReal simtime, int iter, int NOUT);
    int run(MonitorFunc f);

.. [1]
   Taken from a talk by L.Chacon available here
   https://bout2011.llnl.gov/pdf/talks/Chacon_bout2011.pdf

.. [2]
   See paper https://arxiv.org/abs/1209.2054 for an application to
   2-fluid equations

.. [3] This `InvertPar` class can handle cases with closed
   field-lines and twist-shift boundary conditions for tokamak
   simulations
