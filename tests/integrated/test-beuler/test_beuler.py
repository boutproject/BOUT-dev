#!/usr/bin/env python3

# requires: petsc

from boututils.run_wrapper import shell_safe, launch_safe
from contextlib import chdir


nthreads = 1
nproc = 1


def test_beuler(test_dir):
    print("Making solver test")
    shell_safe("make > make.log")

    print("Running solver test")
    with chdir(test_dir):
        status, out = launch_safe("./test_beuler", nproc=nproc, mthread=nthreads, pipe=True)
    with open("run.log", "w") as f:
        f.write(out)

    if status:
        print(out)

    assert status, out
