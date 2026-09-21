#!/usr/bin/env python3

import pytest
import numpy as np
from sys import stdout
from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect

tol = 1e-10  # Absolute tolerance


# The command to run
exefile = "./test_delp2"

# List of settings to apply
settings = [
    "MXG=2 mesh:nx=36 diffusion:useFFT=true",
    "MXG=1 mesh:nx=34 diffusion:useFFT=true",
    "MXG=2 mesh:nx=36 diffusion:useFFT=false",
    "MXG=1 mesh:nx=34 diffusion:useFFT=false",
]


@pytest.mark.parametrize("setting", settings)
def test_delp2(setting):

    # Read benchmark values
    print("Args: " + setting)
    cmd = exefile + " " + setting

    s, out = launch_safe(cmd, nproc=1, pipe=True)
    file_suffix = ".".join([x.split("=")[-1] for x in setting.split()])
    with open("run.log." + str(file_suffix) + ".1", "w") as f:
        f.write(out)

    n0 = collect("n", path="data", info=False)

    for nproc in [2, 4]:
        shell(["rm data/BOUT.dmp.*.nc"])

        stdout.write("   %d processor...." % (nproc))
        s, out = launch_safe(cmd, nproc=nproc, mthread=1, pipe=True)
        with open("run.log." + str(file_suffix) + "." + str(nproc), "w") as f:
            f.write(out)

        # Collect output data
        n = collect("n", path="data", info=False)

        assert np.shape(n) == np.shape(n0), "Fail, wrong shape"

        diff = np.max(np.abs(n - n0))
        assert diff <= tol, "Fail, maximum difference = " + str(diff)
