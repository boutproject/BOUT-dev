Simple LaplaceXY example
========================

This generates a field on a rectangular mesh in X-Y, and inverts it
using LaplaceXY, writing the result to file. Only really useful for
checking if an installation runs, and for trying out different solvers
and preconditioners. See the "ksptype" and "pctype" settings in BOUT.inp

Run with

    $ ./test-laplacexy -ksp_monitor

which should print the KSP norms from PETSc:

     0 KSP Residual norm 5.656854249492e+00
     1 KSP Residual norm 4.732163974221e+00
     2 KSP Residual norm 4.084280618934e+00
     3 KSP Residual norm 3.390335900434e+00
     4 KSP Residual norm 2.980304269384e+00
     5 KSP Residual norm 2.583427730146e+00
     6 KSP Residual norm 2.320399960793e+00
     7 KSP Residual norm 2.059145598820e+00
     8 KSP Residual norm 1.832451815744e+00
     9 KSP Residual norm 1.674179696341e+00
    10 KSP Residual norm 1.589376411329e+00
    11 KSP Residual norm 1.549055878503e+00
    12 KSP Residual norm 1.517041587794e+00
    13 KSP Residual norm 1.473466938498e+00
    14 KSP Residual norm 1.382770759212e+00
    15 KSP Residual norm 1.080408049371e+00
    16 KSP Residual norm 4.309526296050e-01
    17 KSP Residual norm 1.115269396077e-01
    18 KSP Residual norm 4.334487475743e-13
