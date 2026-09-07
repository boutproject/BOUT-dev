#!/usr/bin/env python3

# requires: petsc

from boututils.run_wrapper import shell_safe


def test_beuler():
    print("Running solver test")
    status, out = shell_safe("./test_beuler", pipe=True)
    with open("run.log", "w") as f:
        f.write(out)

    assert status == 0, out
