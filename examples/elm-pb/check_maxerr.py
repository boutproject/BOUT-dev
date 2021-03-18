#!/usr/bin/env python3

from xbout import open_boutdataset

from sys import argv

if len(argv) > 1:
    directory = argv[1]
else:
    directory = "data"

ds = open_boutdataset(f"{directory}/BOUT.dmp.*.nc", keep_xboundaries=False)

ds["relerr2d"] = (ds["Ucheck_2d"] - ds["U"] ) / ds["U"].max(dim=["x", "y", "z"])
ds["relerr3d"] = (ds["Ucheck_3d"] - ds["U"] ) / ds["U"].max(dim=["x", "y", "z"])

print("max relerr2d", abs(ds["relerr2d"]).max(dim=["x", "y", "z"]).values)
print("max relerr3d", abs(ds["relerr3d"]).max(dim=["x", "y", "z"]).values)
