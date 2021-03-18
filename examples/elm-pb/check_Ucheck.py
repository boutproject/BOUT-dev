#!/usr/bin/env python3

from xbout import open_boutdataset

from matplotlib import pyplot as plt
from sys import argv

if len(argv) > 1:
    directory = argv[1]
else:
    directory = "data"

ds = open_boutdataset(f"{directory}/BOUT.dmp.*.nc", keep_xboundaries=False)

ds["relerr2d"] = (ds["Ucheck_2d"] - ds["U"] ) / ds["U"].max(dim=["x", "y", "z"])
ds["relerr3d"] = (ds["Ucheck_3d"] - ds["U"] ) / ds["U"].max(dim=["x", "y", "z"])

plt.rcParams["figure.figsize"] = [8,8]
ds.isel(t=-1).bout.animate_list(["phi_check", "U", "Ucheck_2d", "Ucheck_3d", "relerr2d", "relerr3d"], animate_over="y", show=True, save_as="Ucheck")
