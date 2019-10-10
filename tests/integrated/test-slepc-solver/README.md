Test SLEPc Eigen Solver
=======================

A simple test for the SLEPc eigen solver. Checks that it can get the
eigenvalue of

    ddt(f) = f * exp(t)

which is `0.0 + 1.0i`.

Note that this is just a sanity check and doesn't fully test the
solver.
