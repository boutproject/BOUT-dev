#!/usr/bin/env python3

from xbout import open_boutdataset

from matplotlib import pyplot as plt
from sys import argv

if len(argv) > 1:
    directory = argv[1]
else:
    directory = "data"

ds = open_boutdataset(f"{directory}/BOUT.dmp.*.nc", keep_xboundaries=False)

f, axes = plt.subplots(3,4)
axes = axes.flatten()

for x, ax in zip(
    ["g11", "g22", "g33", "g12", "g13", "g23", "g_11", "g_22", "g_33", "g_12", "g_13", "g_23"],
    axes,
):
    ds[x].plot(ax=ax)

plt.savefig("metrics.pdf")
plt.show()
