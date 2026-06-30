#!/usr/bin/env python3

from boututils.run_wrapper import launch_safe


def test_beuler():
    print("Running solver test")
    status, out = launch_safe("./test_beuler", nproc=1, mthread=1, pipe=True)
    with open("run.log", "w") as f:
        f.write(out)

    assert status == 0, out
