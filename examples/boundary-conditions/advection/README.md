1D test of boundary conditions for the advection equation
=========================================================

This solves the advection equation

   df/dt + df/dx = 0

on a 1D domain with boundaries. An initial Gaussian pulse starts
in the middle of the domain and moves to the right.

The different input files test combinations
of numerical method and boundary condition on the downstream boundary.

Example adapted from this question on StackExchange:

http://scicomp.stackexchange.com/questions/5425/strange-oscillation-when-solving-the-advection-equation-by-finite-difference-wit

upwind
------

Here 1st-order upwinding is used, so the downstream boundary condition
is never used. The pulse can be seen to correctly propagate out of the domain.

central-dirichlet
-----------------

Changing to 2nd-order central differencing, with the default Dirichlet
boundary conditions (2nd order) on both boundaries. 

When the pulse reaches the right boundary a large zigzag mode is formed
which propagates to the left. This is a parasitic (unphysical) mode.

central-free
------------

Changing to a 2nd order "free" boundary condition on the downstream boundary
improves the solution, but a parasitic mode can still be seen propagating to
the left.

central-free-o3
---------------

Using a 3rd-order "free" boundary, which uses three points in the domain to
extrapolate into the guard cell. This results in the correct propagation
out of the domain, without reflecting mode.



