#!/usr/bin/env python3
#
# Python script to run and analyse MMS test
#

# Cores: 2
# only working with cmake
# requires: False
from boututils.run_wrapper import launch_safe
from boututils.datafile import DataFile
from boutdata.collect import collect as _collect

import numpy as np


def collect(var):
    return _collect(
        var,
        info=False,
        path=directory,
        xguards=False,
        yguards=False,
    )


nprocs = [1]  # , 2, 4]
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

# for x in "xout", "xin":
#     regions[x] = np.logical_or(regions[f"{x}_fwd"], regions[f"{x}_bwd"])
# for x in "fwd", "bwd":
#     regions[x] = np.logical_or(regions[f"xin_{x}"], regions[f"xout_{x}"])
# regions["all"] = np.logical_or(regions["xin"], regions["xout"])
for x in "xout", "xin":
    regions[x] = regions[f"{x}_fwd"] + regions[f"{x}_bwd"]
for x in "fwd", "bwd":
    regions[x] = regions[f"xin_{x}"] + regions[f"xout_{x}"]
regions["all"] = regions["xin"] + regions["xout"]

for nproc in nprocs:
    cmd = "./get_par_bndry"

    # Launch using MPI
    _, out = launch_safe(cmd, nproc=nproc, mthread=mthread, pipe=True)

    for k, v in regions.items():
        # Collect data
        data = collect(f"field_{k}")
        assert np.allclose(data, v), (
            k + " does not match",
            np.sum(data),
            np.sum(v),
            np.max(data),
        )
