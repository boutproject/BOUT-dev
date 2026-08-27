#!/usr/bin/env python3

from boutdata.collect import collect
from boututils.run_wrapper import launch_safe
import numpy.testing as npt


def test_slepc_solver():

    print("Running SLEPc eigen solver test")
    status, out = launch_safe("./test-slepc-solver", nproc=1, pipe=True, verbose=True)

    with open("run.log", "w") as f:
        f.write(out)

    eigenvalues = collect("t_array", path="data", info=False)

    expected_eigenvalues = [0.0, 1.0]

    npt.assert_allclose(
        expected_eigenvalues,
        eigenvalues,
        err_msg=" => SLEPc test failed\nEigenvalues: {eigenvalues}",
    )
