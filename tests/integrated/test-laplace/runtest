#!/usr/bin/env python3

#
# Run the test, compare results against the benchmark
#

# requires: not metric_3d
# Requires: netcdf
# Cores: 4

from boututils.run_wrapper import build_and_log, shell, launch_safe
from boutdata.collect import collect, create_cache
import numpy.testing as npt
from sys import exit


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

build_and_log("Laplacian inversion test")

# Read benchmark values
print("Reading benchmark data")
bmk = {v: collect(v, path="data", prefix="benchmark", info=False) for v in vars}

print("Running Laplacian inversion test")
success = True

for solver in ["cyclic", "pcr", "pcr_thomas"]:
    for nproc in [1, 2, 4]:
        nxpe = 1
        if nproc > 2:
            nxpe = 2

        cmd = f"./test_laplace NXPE={nxpe} laplace:type={solver}"

        shell("rm data/BOUT.dmp.*.nc")

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

if success:
    print(" => All Laplacian inversion tests passed")
    exit(0)
else:
    print(" => Some failed tests")
    exit(1)
