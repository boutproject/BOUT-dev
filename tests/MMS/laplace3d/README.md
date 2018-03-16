Test of Laplacian solvers in 3d geometry
----------------------------------------
----------------------------------------


Running
-------
To run the test: ``./runtest``

To generate the input for the test copy the output of the mms.py script to
data/BOUT.inp, i.e. ``python3 mms.py > data/BOUT.inp``. This will need doing if
you change something in mms.py, like the expression for the solution or the
metric coefficients.
Some default options are output in the [laplace] section by the script, but
those can be changed in BOUT.inp without re-running the script.


About
-----
This test is designed to check the convergence of BOUT++'s perpendicular
Laplacian solvers in a 3 dimensional geometry.
The geometry is a simple tokamak with circular cross section, including only the
closed field-line region.
While the solutions on each perpendicular plane are independent, the parallel
variation of the metric (and geometrical coefficients which depend on
derivatives of the metric) affects the solutions and so is useful to test.
