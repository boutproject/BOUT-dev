Constraint equations
====================

These examples show how to include constraints as part of the time advance.
This can be used as a way to solve potential from vorticity, but involves 
adding an auxilliary algebraic equation for the time integrator to solve
alongside the time evolving variables.

Solvers which can include constraints are

* IDA, from the SUNDIALS suite
* IMEX-BDF2, which uses PETSc SNES
