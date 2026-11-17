#!/usr/bin/env python3

from boututils.run_wrapper import launch_safe, shell


def test_unused_options():

    for nproc in [1, 2]:
        print(f"Running with nproc={nproc}...", end="")
        shell(["rm -rf data/BOUT.dmp.*"])
        _s, out = launch_safe(
            "./test_unused_options",
            nproc=nproc,
            mthread=1,
            pipe=True,
        )

        with open(f"run.log.{nproc}", "w") as f:
            f.write(out)

        # Check for printed message
        assert "SUCCESS" in out, "failed"
