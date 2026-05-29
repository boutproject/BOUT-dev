#!/usr/bin/env python3

from boutdata import collect
from matplotlib import pyplot
from sys import argv, exit

datadir = argv[1]

yind = int(argv[2])

xg = 2

f = collect("f", path=datadir)[xg:-xg, :, :]
rhs = collect("rhs", path=datadir)[xg:-xg, :, :]
rhs_check = collect("rhs_check", path=datadir)[xg:-xg, :, :]
error = collect("error", path=datadir)[xg:-xg, :, :]

pyplot.subplot(221)
pyplot.pcolormesh(f[:, yind, :])
pyplot.colorbar()
pyplot.title("f")
pyplot.subplot(222)
pyplot.pcolormesh(rhs[:, yind, :])
pyplot.colorbar()
pyplot.title("rhs")
pyplot.subplot(223)
pyplot.pcolormesh(rhs_check[:, yind, :])
pyplot.colorbar()
pyplot.title("rhs_check")
pyplot.subplot(224)
pyplot.pcolormesh(error[:, yind, :])
pyplot.colorbar()
pyplot.title("error")

pyplot.show()

exit(0)
