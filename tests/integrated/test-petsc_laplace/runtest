#!/usr/bin/env python3

#
# Run the test, compare results against the benchmark
#

# requires: petsc
# requires: all_tests
# cores: 4

from boututils.run_wrapper import build_and_log, shell, launch_safe
from boutdata.collect import collect, create_cache

import pathlib
from sys import exit

errors = [
    "max_error1",
    "max_error2",
    "max_error3",
    "max_error5",
    "max_error6",
    "max_error7",
]
tol = 2e-4  # Absolute (?) tolerance

build_and_log("PETSc Laplacian inversion test")

print("Running PETSc Laplacian inversion test")
success = True

for nproc in [1, 2, 4]:
    cmd = "./test_petsc_laplace"

    shell("rm data/BOUT.dmp.*.nc")

    print(f"   {nproc} processors....")
    s, out = launch_safe(cmd, nproc=nproc, pipe=True, verbose=True)

    pathlib.Path(f"run.log.{nproc}").write_text(out)
    cache = create_cache(path="data", prefix="BOUT.dmp")

    # Collect output data
    for varname in errors:
        print(f"      Checking {varname} ... ", end="")
        error = collect(varname, path="data", info=False, datafile_cache=cache)
        if error <= 0:
            print("Convergence error")
            success = False
        elif error > tol:
            print(f"Fail, maximum error is = {error:e}")
            success = False
        else:
            print("Pass")

if success:
    print(" => All PETSc Laplacian inversion tests passed")
    exit(0)
else:
    print(" => Some failed tests")
    exit(1)
