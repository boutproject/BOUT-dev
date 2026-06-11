#!/usr/bin/env python3
from pathlib import Path

#
# Run the test, compare results against the benchmark
#

# requires: not metric_3d
# Requires: netcdf
# Cores: 4

import pytest
import numpy.testing as npt
from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect, create_cache


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


@pytest.fixture(scope="module")
def benchmark_data():
    """Fixture to load benchmark data once for the entire test session."""

    source_data_dir = str(Path(__file__).parent / "data")
    return {
        v: collect(v, path=source_data_dir, prefix="benchmark", info=False)
        for v in vars
    }


@pytest.mark.parametrize("solver", ["cyclic", "pcr", "pcr_thomas"])
@pytest.mark.parametrize("nproc", [1, 2, 4])
def test_laplace(solver, nproc, benchmark_data):

    nxpe = 2 if nproc > 2 else 1
    cmd = f"./test_laplace NXPE={nxpe} laplace:type={solver}"

    shell("rm data/BOUT.dmp.*.nc")

    s, out = launch_safe(cmd, nproc=nproc, mthread=1, pipe=True)

    with open(f"run.log.{nproc}", "w") as f:
        f.write(out)

    cache = create_cache(path="data", prefix="BOUT.dmp")

    # Collect output data
    for v in vars:
        result = collect(v, path="data", info=False, datafile_cache=cache)
        npt.assert_allclose(
            result,
            benchmark_data[v],
            atol=tol,
            rtol=tol,
            err_msg=f"Failed checking variable: {v}",
        )

    # Only check FieldPerps on one processor because reading them in is
    # quite annoying on multiple cores due to mismatched global y indices
    if nproc == 1:
        for v in ["flag0_perp", "flag3_perp"]:
            result = collect(v, path="data", info=False, datafile_cache=cache)

            # Compare benchmark and output
            npt.assert_allclose(
                result,
                benchmark_data[v.replace("_perp", "")][:, 0, :],
                atol=tol,
                rtol=tol,
                err_msg=f"Failed checking perpendicular variable: {v}",
            )
