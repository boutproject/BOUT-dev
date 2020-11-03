Alfven wave test case
=====================

Solves equations for the vorticity and the parallel electric field:

    d/dt (Vort) =  Div_par(jpar)
    d/dt (A||) = Grad_par(phi)

where `jpar` = `Laplace_perp(A||)` and `Vort = Laplace_perp(phi)`

Switches in the input file change between forms of these operators:

* `split_n0` If true, the axisymmetric component (n=0) is solved
             separately, using `LaplaceXY`. This includes the poloidal
             (Y) derivatives in the vorticity inversion.
* `newXZsolver`: If true, the LaplaceXZ solver is used for the `n != 0`
               components rather than the older Laplacian solver
* `laplace_perp`: If true, `Laplace_perp` is used to calculate `jpar` from
                  `A||`, otherwise `Delp2` is used.

Test cases
----------

Circular cross-section, large aspect ratio plasma:

    $ ./alfven -d cbm18

Tokamak X-point geometry

    $ mpirun -np 8 ./alfven -d data
