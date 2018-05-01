test-solver
===========

Integrate `sin^2(t)` between `t = 0` and `t = pi / 2`. The result
should be `pi / 4`.

This test tries to integrate the above problem with all available
solvers. This is not a complete test of the solvers -- it doesn't
check expected error scaling or convergence rates for example -- but
is more of a basic sanity check: does this solver work at all?
