#!/usr/bin/env python3

#
# Run the test, compare results against the benchmark
#

# requires: not metric_3d
# Requires: netcdf
# Cores: 4

from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect, create_cache
import numpy as np
import numpy.testing as npt
from sys import stdout


# Variables to compare
vars = [
    "flag0",
    "flag3",
    "flagis",
    "flagos",
    "flag0a",
    "flag3a",
    "flagisa",
    "flagosa",
    "flag0ac",
    "flag3ac",
    "flagisac",
    "flagosac",
    "flag0ad",
    "flag3ad",
    "flagisad",
    "flagosad",
]
tol = 1e-6  # Absolute tolerance


# Read benchmark values
print("Reading benchmark data")
bmk = {v: collect(v, path="data", prefix="benchmark", info=False) for v in vars}


def test_laplace():

    # Read benchmark values
    print("Reading benchmark data")
    bmk = {}
    for v in vars:
        bmk[v] = collect(v, path="data", prefix="benchmark", info=False)

    print("Running Laplacian inversion test")
    success = True

    cmd = "./test_laplace NXPE=" + str(nxpe) + " laplace:type=" + solver
       
    print(f"   {solver} solver with {nproc} processors ({nxpe=})....")
    s, out = launch_safe(cmd, nproc=nproc, mthread=1, pipe=True)
    with open(f"run.log.{nproc}", "w") as f:
        f.write(out)

    cache = create_cache(path="data", prefix="BOUT.dmp")

    # Collect output data
    for v in vars:
        print(f"      Checking variable {v} ...", end="")
        result = collect(v, path="data", info=False, datafile_cache=cache)
        # Compare benchmark and output
        try:
        npt.assert_allclose(result, bmk[v], atol=tol, rtol=tol)
        print("Pass")
        except AssertionError as e:
        print(f"Fail: {e}")
        success = False

        # Only check FieldPerps on one processor because reading them in is
        # quite annoying on mutliple cores due to mismatched global y indices
        if nproc == 1:
            for v in ["flag0_perp", "flag3_perp"]:
                print(f"      Checking variable {v} ...", end="")
                result = collect(v, path="data", info=False, datafile_cache=cache)
                # Compare benchmark and output
                try:
                    npt.assert_allclose(
                        result, bmk[v.replace("_perp", "")][:, 0, :], atol=tol, rtol=tol
                    )
                    print("Pass")
                except AssertionError as e:
                    print(f"Fail: {e}")
                    success = False

            # Collect output data
            for v in vars:
                stdout.write("      Checking variable " + v + " ... ")
                result = collect(v, path="data", info=False)
                # Compare benchmark and output
                if np.shape(bmk[v]) != np.shape(result):
                    print("Fail, wrong shape")
                    success = False
                diff = np.max(np.abs(bmk[v] - result))
                if diff > tol:
                    print("Fail, maximum difference = " + str(diff))
                    success = False
                else:
                    print("Pass")

    assert success, " => Some failed tests"
