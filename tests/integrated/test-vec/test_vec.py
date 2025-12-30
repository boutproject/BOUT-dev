#!/usr/bin/env python3

# Cores: 4

from boututils.run_wrapper import build_and_log, launch_safe

def test_vec():

    build_and_log("vector communication test")

    status, out = launch_safe("./testVec", nproc=4, pipe=True)

    with open("run.log.4", "w") as f:
        f.write(out)

    assert not status, " => Test failed"
