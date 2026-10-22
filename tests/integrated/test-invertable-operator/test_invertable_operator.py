#!/usr/bin/env python3
#
# Run the test, compare results against expected value
#

import pytest
from boutdata.collect import collect
from boututils.run_wrapper import launch_safe, shell

nprocs = [1, 2]  # Number of processors to run on
reltol = 1.0e-3  # Allowed relative tolerance
nthreads = 1


# Delete old output files
shell(["rm -rf data/BOUT.dmp.*"])


def run(path, nproc, log=False):
    pipe = bool(log)
    _s, out = launch_safe(
        "./invertable_operator -d " + path, nproc=nproc, mthread=nthreads, pipe=pipe
    )
    if pipe:
        with open(log, "w") as f:
            f.write(out)

    # Get result of the test
    passVerification = collect("passVerification", path=path)[-1]
    maxRelErrLaplacians = collect("maxRelErrLaplacians", path=path)[-1]

    if passVerification == 0:
        print(f"  => Failed (verification step - value is {passVerification})")
        pytest.fail()

    if maxRelErrLaplacians > reltol:
        print(
            f"  => Failed (relative tolerance step -- difference of {maxRelErrLaplacians})"
        )
        pytest.fail()


def test_invertable_operator():

    for np in nprocs:
        run("data", np)
