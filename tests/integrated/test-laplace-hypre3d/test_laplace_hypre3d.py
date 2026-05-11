#!/usr/bin/env python3

# requires: hypre

from boutdata import collect
from boututils.run_wrapper import launch_safe

test_directories = [
    ("data_slab_core", 1),
    ("data_slab_sol", 1),
    ("data_circular_core", 1),
    ("data_circular_core-sol", 1),
]

tolerance = 1.0e-6


def test_laplace_hypre3d():

    success = True
    for directory, nproc in test_directories:
        command = "test-laplace3d -d " + directory
        print("running on", nproc, "processors:", command)
        launch_safe(command, nproc=nproc)

        error_max = collect("error_max", path=directory, info=False)

        if error_max > tolerance:
            print(directory + " failed with maximum error {}".format(error_max))
            success = False
        else:
            print(directory + " passed with maximum error {}".format(error_max))

    assert success
