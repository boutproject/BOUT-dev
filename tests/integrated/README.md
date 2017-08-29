# BOUT++ Integrated tests

This set of tests are designed to test that different components of
the BOUT++ library work together. These tests are more expensive than
the unit tests, but are expected to be run on at least every pull
request, and the majority on every commit.

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

* [**test-io**][test-io] This just reads data from input, and writes
  it out to a file. If this test fails, then all following tests will
  probably also fail. Common causes of failure are problems with the
  NetCDF library, or with Python routines for reading data (see user
  manual).
* [**test-restarting**][test-restarting] Tests that a simulation can
  be restarted correctly
* [**test-fieldfactory**][test-fieldfactory] This checks that the
  FieldFactory class can generate values from analytic expressions
  correctly. Since this functionality is used in nearly all other
  tests, this is crucial.
* [**test-laplace**][test-laplace] Tests the Laplacian inversion code
* [**test-cyclic**][test-cyclic] Tests the tridiagonal solver
* [**test-invpar**][test-invpar] Tests the parallel parabolic solver,
  which depends on test-cyclic.
* [**test-smooth**][test-smooth] Tests smoothing operators
* [**test-gyro**][test-gyro] Tests gyro-averaging operators, mainly
  used in gyrofluid models
* [**test-delp2**][test-delp2] Tests the second derivative Laplacian
  operator
* [**test-vec**][test-vec] Tests Vector3D communication
* [**test-griddata**][test-griddata] Test if variable references in
  the input file work correctly
* [**test-initial**][test-initial] Test initial conditions are created
  correctly
* [**test-drift-instability**][test-drift-instability] calculates the
  frequency and growth-rate of a resistive drift wave in a slab,
  comparing the result against analytic and reference runs
* [**test-interchange-instability**][test-interchange-instability]
  calculates the growth-rate for an interchange mode in a curved slab,
  for two different curvature radii. Also compares the result against
  analytic and reference values.

[test-io]: test-io/README.md
[test-fieldfactory]: test-fieldfactory/README.md
[test-laplace]: test-laplace/README.md
[test-cyclic]: test-cyclic/README.md
[test-invpar]: test-invpar/README.md
[test-smooth]: test-smooth/README.md
[test-gyro]: test-gyro/README.md
[test-delp2]: test-delp2/README.md
[test-drift-instability]: test-drift-instability/README.md
[test-interchange-instability]: test-interchange-instability/README.md
[test-restarting]: test-restarting/README.md
[test-vec]: test-vec/README.md
[test-griddata]: test-griddata/README.md
[test-initial]: test-initial/README.md
