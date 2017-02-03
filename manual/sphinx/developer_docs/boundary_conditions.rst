Boundary conditions
===================

The boundary condition system needs to be very flexible in order to
handle:

-  Meshes which can divide up the boundary into an arbitrary number of
   regions, giving each one a label. For example in BoutMesh (specific
   to tokamaks), the boundary regions are labelled “core”, “sol”, “pf”
   and “target”.

-  Each variable can have a different boundary condition in each region.
   It should be possible to have a global setting “all variables have
   dirichlet conditions on all boundaries”, which is over-ridden by more
   specific settings such as “All variables have neumann conditions on
   the inner x boundaries”, and finally to “variable ’Ni’ has laplacian
   boundary conditions in the ’sol’ regions”

-  Boundary conditions can be modified to be “relaxing”. This means that
   rather than enforcing a strict boundary condition, it’s a mixture of
   zero-gradient in the time-derivative combined with a damping
   (relaxation) towards the desired boundary condition. This can help
   improve the numerics of turbulence simulations.

-  Users should be able to implement their own boundary conditions, and
   add them to the system at run-time without modifying the core code.

-  After ``physics_init``, a boundary condition must be applied to the
   variables. During a simulation (at the end of ``physics_run``),
   boundary conditions need to be applied to the time-derivatives. The
   boundary system should ensure that these conditions are consistent.

