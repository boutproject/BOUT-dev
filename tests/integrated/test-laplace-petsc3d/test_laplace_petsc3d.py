#!/usr/bin/env python3

# requires: petsc

from boutdata import collect
from boututils.run_wrapper import launch


test_directories = [
    ("data_slab_core", 1),
    ("data_slab_sol", 1),
    ("data_circular_core", 1),
    ("data_circular_core-sol", 1),
]

tolerance = 1.0e-6


def test_laplace_petsc3d():

    errors = {}

    for directory, nproc in test_directories:
        command = f"./test-laplace3d -d {directory}"
        print("running on", nproc, "processors:", command)
        status, output = launch(command, nproc=nproc, pipe=True)

        if status:
            print("FAILED")
            print(output)
            errors[directory] = "<bad exit>"
            continue

        error_max = collect("error_max", path=directory, info=False)
        if error_max > tolerance:
            errors[directory] = error_max

    print("**********")

    if errors:
        print("Some failures:")
        for name, error in errors.items():
            print(f"{name}, max error: {error}")

    assert not errors
