#!/usr/bin/env python3
from pathlib import Path
from sys import stdout

import numpy as np

#
# Run the test, compare results against the benchmark
#
import pytest
from boutdata.collect import collect
from boututils.run_wrapper import launch_safe, shell

# Variables to compare
variables = ["pade1", "pade2"]

tol = 1e-7  # Absolute tolerance, benchmark values are floats


@pytest.fixture(scope="module")
def benchmark_data():
    """Reads the benchmark data once for all test iterations in this module."""
    bmk = {}
    source_data_dir = str(Path(__file__).parent / "data")

    for v in variables:
        bmk[v] = collect(
            v, path=source_data_dir, prefix="benchmark", info=False, xguards=False
        )
    return bmk


def test_gyro(benchmark_data):

    print("Running Gyro-average inversion test")

    for nproc in [1, 2, 4]:
        nxpe = 1
        if nproc > 2:
            nxpe = 2

        cmd = f"./test_gyro NXPE={nxpe}"

        shell("rm -rf data/BOUT.dmp.*")

        print(f"   {nproc} processors (nxpe = {nxpe})....")
        _s, out = launch_safe(cmd, nproc=nproc, pipe=True)
        with open("run.log." + str(nproc), "w") as f:
            f.write(out)

        # Collect output data
        for v in variables:
            stdout.write("      Checking variable " + v + " ... ")
            result = collect(v, path="data", info=False, xguards=False)
            # Compare benchmark and output
            np.testing.assert_allclose(
                benchmark_data[v],
                result,
                atol=tol,
                rtol=0,
                err_msg="Gyro-average output mismatch",
            )
