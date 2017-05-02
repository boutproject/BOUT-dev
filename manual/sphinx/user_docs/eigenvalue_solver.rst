Eigenvalue solver
=================

By using the SLEPc library, BOUT++ can be used as an eigenvalue solver
to find the eigenvectors and eigenvalues of sets of equations.

Configuring with SLEPc
----------------------

The BOUT++ interface has been tested with SLEPc version 3.4.3, itself
compiled with PETSc 3.4.2. SLEPc version 3.4 should work, but other
versions will not yet.

SLEPc options
-------------

Time derivatives can be taken directly from the RHS function, or by
advancing the simulation in time by a relatively large increment. This
second method acts to damp high frequency components

Examples
--------

Wave in a box
~~~~~~~~~~~~~

``examples/eigen-box``

