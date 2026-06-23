#!/usr/bin/env python3

#
# Run the test, compare results against the benchmark
#

# requires: all_tests
# requires: petsc
# cores: 4

import pytest
from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect

# Variables to compare
vars = [
    ("max_error1", 2.0e-4),
    ("max_error2", 1.0e-4),
    ("max_error3", 1.0e-4),
    ("max_error4", 1.0e-4),
    ("max_error5", 2.0e-3),
    ("max_error6", 3.0e-4),
    ("max_error7", 2.0e-4),
    ("max_error8", 1.0e-4),
]
# tol = 1e-4                  # Absolute (?) tolerance


@pytest.mark.parametrize("nproc", [1, 2, 4])
@pytest.mark.parametrize("jy", [2, 34, 65, 81, 113])
def test_petsc_laplace_MAST_grid(nproc, jy):
    print(
        "Running PETSc Laplacian inversion test with non-identity metric (taken from grid for MAST SOL)"
    )

    cmd = f"./test_petsc_laplace_MAST_grid mesh:file=grids/grid_MAST_SOL_jyis{jy}.nc"

    shell(["rm -f data/BOUT.dmp.*.nc"])

    s, out = launch_safe(cmd, nproc=nproc, pipe=True)

    with open(f"run.log.{nproc}.jy_{jy}", "w") as f:
        f.write(out)

    # Collect output data
    failures = []
    for var_name, tolerance in vars:
        error = collect(var_name, path="data", info=False)

        if error <= 0:
            failures.append(f"{var_name}: Solver did not converge (error = {error})")
        elif error > tolerance:
            failures.append(
                f"{var_name}: Absolute error {error:.4e} exceeds tolerance {tolerance:.4e}"
            )

    assert not failures, "Some failed tests:\n" + "\n".join(failures)
