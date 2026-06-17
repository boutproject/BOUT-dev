#!/usr/bin/env python3

#
# Run the test, check the error
#

from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect

tol = 1e-9  # Absolute tolerance
numTests = 4  # We test 4 different boundary conditions (with slightly different inputs for each)


def test_multigrid_laplace_unsheared():

    print("Running multigrid Laplacian inversion test (unsheared)")
    success = True

    for nproc in [1, 3]:
        # Make sure we don't use too many cores:
        # Reduce number of OpenMP threads when using multiple MPI processes
        mthread = 2
        if nproc > 1:
            mthread = 1

        # set nxpe on the command line as we only use solution from one point in y,
        # so splitting in y-direction is redundant (and also doesn't help test the multigrid solver)
        inputfile = "BOUT_unsheared.inp"
        cmd = f"./test_multigrid_laplace -f {inputfile} NXPE={nproc}"

        shell(["rm data/BOUT.dmp.*.nc"])

        print("   %d processors..." % nproc)
        s, out = launch_safe(cmd, nproc=nproc, mthread=mthread, pipe=True)
        with open("run.log." + str(nproc), "w") as f:
            f.write(out)

        # Collect errors
        errors = [
            collect("max_error" + str(i), path="data") for i in range(1, numTests + 1)
        ]

        for i, e in enumerate(errors):
            print("Checking test " + str(i))
            if e < 0.0:
                print("Fail, solver did not converge")
                success = False
            if e > tol:
                print("Fail, maximum absolute error = " + str(e))
                success = False
            else:
                print("Pass")

    assert success, " => Some failed tests"
