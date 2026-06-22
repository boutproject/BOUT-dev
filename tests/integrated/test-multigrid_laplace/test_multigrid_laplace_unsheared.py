#!/usr/bin/env python3

#
# Run the test, check the error
#

import pytest
from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect

tol = 1e-9  # Absolute tolerance
numTests = 4  # We test 4 different boundary conditions (with slightly different inputs for each)


@pytest.mark.parametrize("nproc", [1, 3])
def test_multigrid_laplace_unsheared(nproc):
    # Make sure we don't use too many cores:
    # Reduce number of OpenMP threads when using multiple MPI processes
    mthread = 1 if nproc > 1 else 2

    # set nxpe on the command line as we only use solution from one point in y,
    # so splitting in y-direction is redundant (and also doesn't help test the multigrid solver)
    inputfile = "BOUT_unsheared.inp"
    cmd = f"./test_multigrid_laplace -f {inputfile} NXPE={nproc}"

    shell(["rm -f data/BOUT.dmp.*.nc"])

    print(
        f"Running unsheared multigrid Laplacian test: {nproc} procs, {mthread} thread(s)"
    )

    s, out = launch_safe(cmd, nproc=nproc, mthread=mthread, pipe=True)

    with open(f"run.log.{nproc}", "w") as f:
        f.write(out)

    # Collect errors
    failures = []
    for i in range(1, numTests + 1):
        var_name = f"max_error{i}"
        e = collect(var_name, path="data", info=False)

        if e < 0.0:
            failures.append(
                f"Sub-test {i} ({var_name}): Solver did not converge (error = {e})"
            )
        elif e > tol:
            failures.append(
                f"Sub-test {i} ({var_name}): Absolute error {e:.4e} exceeds tolerance {tol:.4e}"
            )

    assert not failures, "Some failed tests:\n" + "\n".join(failures)
