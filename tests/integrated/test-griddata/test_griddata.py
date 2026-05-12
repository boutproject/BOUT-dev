#!/usr/bin/env python3

from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect
import numpy.testing as npt
from sys import stdout


def test_griddata():
    for nproc in [1]:
        stdout.write("Checking %d processors ... " % (nproc))

        shell("rm ./data*nc")
        s, out = launch_safe("./test_griddata -d screw", nproc=nproc, pipe=True)

        with open("run.log." + str(nproc), "w") as f:
            f.write(out)

        prefix = "data"
        Rxy = collect("Rxy", prefix=prefix, info=False)
        Bpxy = collect("Bpxy", prefix=prefix, info=False)
        dx = collect("dx", prefix=prefix, info=False)

        nx, ny = Rxy.shape

        # Handle 3D metric case
        if len(dx.shape) == 3:
            dx = dx[:, :, 0]

        rwidth = 0.4
        dr = float(rwidth) / nx

        # Test value of dx
        npt.assert_allclose(
            dx, dr * Bpxy * Rxy, atol=1e-7, err_msg="Failed: dx does not match"
        )

        print("Passed")
