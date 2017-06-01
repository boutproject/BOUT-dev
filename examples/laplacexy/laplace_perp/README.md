Comparison of LaplaceXY and Laplace_perp
========================================

These operators are discretised in different ways,
but should be inverses of each other for axisymmetric modes
(no z dependence). This test case takes an analytic input,
integrates using LaplaceXY, differentiates with Laplace_perp,
and compares the result against the input.


To run the test case
    $ ./runtest

By default this will run the "torus" test. Modify the path in runtest
to "square" to run the test case in a square domain.

Example output
--------------

    Magnitude of orig:  0.999699
    Maximum difference:  0.0015555

Turning off y derivatives in LaplaceXY gives

    Maximum difference:  0.625351

Turning off non_uniform corrections in Laplace_perp gives

    Maximum difference:  0.197387
