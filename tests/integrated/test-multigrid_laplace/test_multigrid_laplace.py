#!/usr/bin/env python3

#
# Run the test, check the error
#

import pytest
from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect

tol = 2e-6  # Absolute tolerance
numTests = 4  # We test 4 different boundary conditions (with slightly different inputs for each)


@pytest.mark.parametrize("nproc", [1, 2, 4])
@pytest.mark.parametrize(
    "inputfile", ["BOUT_jy4.inp", "BOUT_jy63.inp", "BOUT_jy127.inp"]
)
def test_multigrid_laplace(nproc, inputfile):
    # set nxpe on the command line as we only use solution from one point in y,
    # so splitting in y-direction is redundant (and also doesn't help test the multigrid solver)
    cmd = f"./test_multigrid_laplace -f {inputfile} NXPE={nproc} input:error_on_unused_options=false"

    shell(["rm -f data/BOUT.dmp.*.nc"])

    print("Running multigrid Laplacian inversion test")
    print(f"{nproc} processors, input file is: {inputfile}")

    s, out = launch_safe(cmd, nproc=nproc, pipe=True)

    # Save a log file uniquely named for this specific combination
    input_file_name = inputfile.replace(".inp", "")
    with open(f"run.log.{nproc}.{input_file_name}", "w") as f:
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
                f"Sub-test {i} ({var_name}): Absolute error {e} exceeds tolerance {tol}"
            )

    assert not failures, "Some failed tests:\n" + "\n".join(failures)
