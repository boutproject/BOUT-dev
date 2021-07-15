#!/usr/bin/env python3

from boutdata import collect
from boututils.datafile import DataFile
from matplotlib import pyplot
import numpy
from sys import exit

# Set default colormap
pyplot.rcParams["image.cmap"] = "viridis"
# pyplot.rcParams['image.cmap'] = 'plasma'

g11 = collect("g11", path="data", xguards=False, yguards=False)
g12 = collect("g12", path="data", xguards=False, yguards=False)
g13 = collect("g13", path="data", xguards=False, yguards=False)
g22 = collect("g22", path="data", xguards=False, yguards=False)
g23 = collect("g23", path="data", xguards=False, yguards=False)
g33 = collect("g33", path="data", xguards=False, yguards=False)
g_11 = collect("g_11", path="data", xguards=False, yguards=False)
g_12 = collect("g_12", path="data", xguards=False, yguards=False)
g_13 = collect("g_13", path="data", xguards=False, yguards=False)
g_22 = collect("g_22", path="data", xguards=False, yguards=False)
g_23 = collect("g_23", path="data", xguards=False, yguards=False)
g_33 = collect("g_33", path="data", xguards=False, yguards=False)

grid = DataFile("data/d3d_119919.nc")
R = grid["Rxy"][2:-2, :]
Z = grid["Zxy"][2:-2, :]

pyplot.subplot(231)
pyplot.pcolor(R, Z, g11)
pyplot.colorbar()
pyplot.title("g11")

pyplot.subplot(232)
pyplot.pcolor(R, Z, g12)
pyplot.colorbar()
pyplot.title("g12")

pyplot.subplot(233)
pyplot.pcolor(R, Z, g13)
pyplot.colorbar()
pyplot.title("g13")

pyplot.subplot(234)
pyplot.pcolor(R, Z, g22)
pyplot.colorbar()
pyplot.title("g22")

pyplot.subplot(235)
pyplot.pcolor(R, Z, g23)
pyplot.colorbar()
pyplot.title("g23")

pyplot.subplot(236)
pyplot.pcolor(R, Z, g33)
pyplot.colorbar()
pyplot.title("g33")


pyplot.figure()

pyplot.subplot(231)
pyplot.pcolor(R, Z, g_11)
pyplot.colorbar()
pyplot.title("g_11")

pyplot.subplot(232)
pyplot.pcolor(R, Z, g_12)
pyplot.colorbar()
pyplot.title("g_12")

pyplot.subplot(233)
pyplot.pcolor(R, Z, g_13)
pyplot.colorbar()
pyplot.title("g_13")

pyplot.subplot(234)
pyplot.pcolor(R, Z, g_22)
pyplot.colorbar()
pyplot.title("g_22")

pyplot.subplot(235)
pyplot.pcolor(R, Z, g_23)
pyplot.colorbar()
pyplot.title("g_23")

pyplot.subplot(236)
pyplot.pcolor(R, Z, g_33)
pyplot.colorbar()
pyplot.title("g_33")

pyplot.show()

exit(0)
