Overview of BOUT++
==================

BOUT++ is a multiblock structured finite difference (/volume) code in
curvilinear coordinates, with some features to support unusual
coordinate systems used in fusion plasma physics. An important point
is that it treats time integration and spatial operators separately,
an approach called the Method of Lines (MOL). This means that BOUT++
consists of two main parts:

#. A set of Ordinary Differential Equation (ODE) integrators, including implicit, explicit and IMEX schemes,
   such as Runge-Kutta and the CVODE solver from SUNDIALS. These don't "know" anything about
   the equations being solved, only requiring the time derivative of the system state.
   For example they make no distinction between the different evolving fields, or the
   number of dimensions in the simulation. This kind of problem-specific information
   can be used to improve efficiency, and is usually supplied in the form of user-supplied
   preconditioners. See section :ref:`sec-timeoptions` for more details.

#. A set of operators and data types for calculating time derivatives, given the system state.
   These calculate things like algebraic operations (+,-,*,/ etc), spatial derivatives,
   and some integral operators.  


Each of these two parts treats the other as a black box (mostly), and they communicate by
exchanging arrays of data: The ODE integrator finds the system state at a given time and passes it
to the problem-dependent code, which uses a combination of operators to calculate the time derivative.
This time derivative is passed back to the ODE integrator, which updates the state and the cycle continues.

This scheme has some advantages in terms of flexibility: Each part of the code doesn't depend on the
details of the other, so can be changed without requiring modifications to the other. Unfortunately
for many problems the details can make a big difference, so ways to provide problem-specific
information to time integrators, such as preconditioners, are also provided.

Next: :ref:`Getting Started <sec-getting-started>`

