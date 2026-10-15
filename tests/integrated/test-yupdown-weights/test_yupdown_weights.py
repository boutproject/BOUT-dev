#!/usr/bin/env python3

from sys import stdout

from boutdata.collect import collect
from boututils.run_wrapper import launch_safe, shell
from numpy import abs, max


def test_yupdown_weights():

    failed = False
    for shifttype in ["shiftedinterp"]:
        shell(["rm -rf data/BOUT.dmp.*"])
        _s, out = launch_safe(
            "./test_yupdown_weights mesh:paralleltransform:type=" + shifttype,
            nproc=1,
            pipe=True,
            verbose=True,
        )

        with open("run.log", "w") as f:
            f.write(out)

        vars = [("ddy", "ddy2")]
        for v1, v2 in vars:
            stdout.write(f"Testing {v1} and {v2} ... ")
            ddy = collect(v1, path="data", xguards=False, yguards=False, info=False)
            ddy2 = collect(v2, path="data", xguards=False, yguards=False, info=False)

            diff = max(abs(ddy - ddy2))

            if diff < 1e-8:
                print(shifttype + f" passed (Max difference {diff:e})")
            else:
                print(shifttype + f" failed (Max difference {diff:e})")
                failed = True

    assert not failed
