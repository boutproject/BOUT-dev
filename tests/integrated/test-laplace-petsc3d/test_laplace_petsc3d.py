#!/usr/bin/env python3

import pytest
from boutdata import collect
from boututils.run_wrapper import launch_safe

test_directories = [
    ("data_slab_core", 1),
    ("data_slab_sol", 1),
    ("data_circular_core", 1),
    ("data_circular_core-sol", 1),
]

tolerance = 1.0e-6


@pytest.mark.parametrize("directory, nproc", test_directories)
def test_laplace_petsc3d(directory, nproc):
    command = f"./test-laplace3d -d {directory}"
    print("running on", nproc, "processors:", command)
    launch_safe(command, nproc=nproc, pipe=True)

    error_max = collect("error_max", path=directory, info=False)
    assert error_max < tolerance
