Solver
======

The solver is the interface between BOUT++ and the time-integration code
such as SUNDIALS. All solvers implement the ``Solver`` class interface
(see ``src/solver/generic_solver.hxx``).

First all the fields which are to be evolved need to be added to the
solver. These are always done in pairs, the first specifying the field,
and the second the time-derivative:

::

    void add(Field2D &v, Field2D &F_v, const char* name);

This is normally called in the ``physics_init`` initialisation routine.
Some solvers (e.g. IDA) can support constraints, which need to be added
in the same way as evolving fields:

::

    bool constraints();
    void constraint(Field2D &v, Field2D &C_v, const char* name);

The ``constraints()`` function tests whether or not the current solver
supports constraints. The format of ``constraint(...)`` is the same as
``add``, except that now the solver will attempt to make ``C_v`` zero.
If ``constraint`` is called when the solver doesn’t support them then an
error should occur.

If the physics model implements a preconditioner or Jacobian-vector
multiplication routine, these can be passed to the solver during
initialisation:

::

    typedef int (*PhysicsPrecon)(BoutReal t, BoutReal gamma, BoutReal delta);
    void setPrecon(PhysicsPrecon f); // Specify a preconditioner
    typedef int (*Jacobian)(BoutReal t);
    void setJacobian(Jacobian j); // Specify a Jacobian

If the solver doesn’t support these functions then the calls will just
be ignored.

Once the problem to be solved has been specified, the solver can be
initialised using:

::

    int init(rhsfunc f, int argc, char **argv, bool restarting, int nout, BoutReal tstep);

which returns an error code (0 on success). This is currently called in
:doc:`bout++.cxx<../_breathe_autogen/file/bout_09_09_8cxx>`:

::

    if(solver.init(physics_run, argc, argv, restart, NOUT, TIMESTEP)) {
      output.write("Failed to initialise solver. Aborting\n");
      return(1);
    }

which passes the (physics module) RHS function ``physics_run`` to the
solver along with the number and size of the output steps.

To run the solver using the (already supplied) settings, there is the
function:

::

    typedef int (*MonitorFunc)(BoutReal simtime, int iter, int NOUT);
    int run(MonitorFunc f);

Implementation: PVODE
---------------------

Implementation: IDA
-------------------

Implementation: PETSc
---------------------

