#BOUT++ examples

## Getting started

Some simple examples to get started are

* **conduction** This is discussed in the user manual, and solves a 1D diffusion equation
* **blob2d** A 2D simulation of a plasma blob. This is a more realistic simulation, simple enough
  to start with, but can be used as a starting point for novel studies
* **hasegawa-wakatani** Solves the H-W equations (optionally with modified zonal response) in 2D
* **test-staggered** shows how to use staggered grids in BOUT++

## Test suite

To run the test suite, run the Python code "test_suite" :

```bash
$ ./test_suite
```

This should run through the tests, provided that python is set up correctly.
If you encounter problems, first check that PYTHONPATH includes the BOUT++
pylib directory:

```bash
$ echo $PYTHONPATH
```

The test suite currently includes:

* **test-io**  This just reads data from input, and writes it out to a file. If this test
  fails, then all following tests will probably also fail. Common causes of failure are
  problems with the NetCDF library, or with Python routines for reading data (see user manual).
* **test-fieldfactory** This checks that the FieldFactory class can generate values from analytic
  expressions correctly. Since this functionality is used in nearly all other tests, this is crucial.
* **test-laplace** Tests the Laplacian inversion code
* **test-cyclic** Tests the tridiagonal solver
* **test-invpar** Tests the parallel parabolic solver, which depends on test-cyclic.
* **test-smooth** Tests smoothing operators
* **test-gyro** Tests gyro-averaging operators, mainly used in gyrofluid models
* **test-delp2** Tests the second derivative Laplacian operator
* **MMS/diffusion** is a Method of Manufactured Solutions check of convergence for a diffusion equation in 1D
* **MMS/wave-1d** is an MMS test of a wave equation in X, including Neumann and Dirichlet boundaries
* **MMS/wave-1d-y** is an MMS test of a wave equation in Y, including Neumann and Dirichlet boundaries
* **drift-instability** calculates the frequency and growth-rate of a resistive drift wave in a slab, comparing the result against analytic and reference runs
* **interchange-instability** calculates the growth-rate for an interchange mode in a curved slab, for two different curvature radii. Also compares the result against analytic and reference values.



