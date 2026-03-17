#!/usr/bin/env python3

# requires: slepc

from boutdata.collect import collect
from boututils.run_wrapper import launch_safe
from numpy import isclose


def test_slepc_solver():

    print("Running SLEPc eigen solver test")
    status, out = launch_safe("./test-slepc-solver", nproc=1, pipe=True, verbose=True)

    with open("run.log", "w") as f:
        f.write(out)

    eigenvalues = collect("t_array", path="data", info=False)

    expected_eigenvalues = [0.0, 1.0]

    assert isclose(expected_eigenvalues, eigenvalues).all(), " => SLEPc test failed\nEigenvalues: {eigenvalues}"
