#!/usr/bin/env python3

#
# Run the test, compare results against the benchmark
#

# Requires: netcdf
# Cores: 4

import pytest
import numpy as np
import itertools
from pathlib import Path
from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect

# Variables to compare
vars = ["yavg2d", "yavg3d", "sm3d"]
tol = 1e-7  # Absolute tolerance, benchmark values are floats


@pytest.fixture(scope="module")
def benchmark_data():
    """Reads the benchmark data once for all test iterations in this module."""
    bmk = {}
    source_data_dir = str(Path(__file__).parent / "data")

    for v in vars:
        bmk[v] = collect(v, path=source_data_dir, prefix="benchmark", info=False)
    return bmk


# Generate test configurations: [(1, 1), (1, 2), (2, 1), (2, 2)]
PROCESSOR_TOPOLOGIES = list(itertools.product([1, 2], [1, 2]))


@pytest.mark.parametrize("nxpe, nype", PROCESSOR_TOPOLOGIES)
def test_smooth(benchmark_data, nxpe, nype):
    nproc = nxpe * nype
    cmd = "./test_smooth"

    # Clean up old data
    shell(["rm -f data/BOUT.dmp.*.nc"])

    # Run the executable
    s, out = launch_safe(f"{cmd} NXPE={nxpe}", nproc=nproc, pipe=True)

    # Save log
    with open(f"run.log.{nproc}", "w") as f:
        f.write(out)

    # Collect output data
    for v in vars:
        result = collect(v, path="data", info=False)
        bmk = benchmark_data[v]

        # Compare benchmark and output

        assert bmk.shape == result.shape, (
            f"Shape mismatch for variable '{v}' on {nxpe}x{nype} grid"
        )

        np.testing.assert_allclose(
            result,
            bmk,
            atol=tol,
            rtol=0,
            err_msg=f"Data mismatch for variable '{v}' on {nxpe}x{nype} grid",
        )
