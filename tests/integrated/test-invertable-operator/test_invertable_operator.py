#!/usr/bin/env python3

# requires: petsc
# requires: not metric_3d
# cores: 2

#
# Run the test, compare results against expected value
#

import pytest
from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect

nprocs = [1, 2]  # Number of processors to run on
reltol = 1.0e-3  # Allowed relative tolerance
nthreads = 1


# Delete old output files
shell(["rm data/BOUT.dmp.*"])


def run(path, nproc, log=False):
    pipe = False
    if log != False:
        pipe = True
    s, out = launch_safe(
        "./invertable_operator -d " + path, nproc=nproc, mthread=nthreads, pipe=pipe
    )
    if pipe:
        f = open(log, "w")
        f.write(out)
        f.close()

    # Get result of the test
    passVerification = collect("passVerification", path=path)[-1]
    maxRelErrLaplacians = collect("maxRelErrLaplacians", path=path)[-1]

    if passVerification == 0:
        print(
            "  => Failed (verification step - value is {s})".format(s=passVerification)
        )
        pytest.fail()

    if maxRelErrLaplacians > reltol:
        print(
            "  => Failed (relative tolerance step -- difference of {s})".format(
                s=maxRelErrLaplacians
            )
        )
        pytest.fail()


def test_invertable_operator():

    for np in nprocs:
        run("data", np)
