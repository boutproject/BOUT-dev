#!/usr/bin/env python3
#
# Python script to run and analyse MMS test

from boututils.run_wrapper import launch_safe
from boututils.datafile import DataFile
from boutdata.collect import collect

import numpy as np


nprocs = [1]
mthread = 2

directory = "data"

with DataFile("grid.fci.nc") as grid:
    MXG = grid.get("MXG", default=1)
    xfwd = grid.read("forward_xt_prime")[MXG:-MXG]
    xbwd = grid.read("backward_xt_prime")[MXG:-MXG]

nx = xfwd.shape[0]

regions = {
    "xin_fwd": xfwd < MXG,
    "xout_fwd": xfwd > nx + MXG - 1,
    "xin_bwd": xbwd < MXG,
    "xout_bwd": xbwd > nx + MXG - 1,
}
regions = {k: v.astype(int) for k, v in regions.items()}

for x in "xout", "xin":
    regions[x] = regions[f"{x}_fwd"] + regions[f"{x}_bwd"]
for x in "fwd", "bwd":
    regions[x] = regions[f"xin_{x}"] + regions[f"xout_{x}"]
regions["all"] = regions["xin"] + regions["xout"]

for nproc in nprocs:
    cmd = "./get_par_bndry"
    _, out = launch_safe(cmd, nproc=nproc, mthread=mthread, pipe=True)

    for k, v in regions.items():
        data = collect(
            f"field_{k}",
            info=False,
            path=directory,
            xguards=False,
            yguards=False,
        )
        assert np.allclose(data, v), (
            f"{k} does not match",
            np.sum(data),
            np.sum(v),
            np.max(data),
        )
