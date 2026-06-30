#!/usr/bin/env python3

#
# Run the test, compare results against the benchmark
#

import numpy as np
import pytest

from boutdata.collect import collect
from boututils.run_wrapper import shell, launch_safe, getmpirun


@pytest.mark.parametrize("nproc", [1, 2, 4])
def test_laplacexz(nproc):

    tol = 1e-10  # Absolute tolerance

    MPIRUN = getmpirun()

    print(f"Running LaplaceXZ test on {nproc} processors...")

    # Unique data directory per processor configuration
    cmd = f"./test-laplacexz nxpe={nproc}"

    shell(["rm -f data/BOUT.dmp.*.nc"])

    print(f"   {nproc} processors (nxpe = {nproc})....")
    s, out = launch_safe(cmd, runcmd=MPIRUN, nproc=nproc, mthread=1, pipe=True)
    with open(f"run.log.{nproc}", "w") as f:
        f.write(out)

    # Collect output data
    f = collect("f", path="data", info=False)
    f2 = collect("f2", path="data", info=False)
    print("      Checking tolerance... ")
    # Compare benchmark and output
    np.testing.assert_allclose(
        f2, f, atol=tol, rtol=0, err_msg="LaplaceXZ output mismatch"
    )
