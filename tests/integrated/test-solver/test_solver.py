#!/usr/bin/env python3

from boututils.run_wrapper import shell_safe


def test_solver():
    print("Running solver test")
    status, out = shell_safe("./test_solver", pipe=True)
    with open("run.log", "w") as f:
        f.write(out)

    assert status == 0, out
