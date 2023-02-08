#!/usr/bin/env python3

from boutdata.data import BoutOutputs
from matplotlib import pyplot
from sys import exit

Ntests = 4

d = BoutOutputs("data")
# d = BoutOutputs("data", caching=True, info=False)

for i in range(1, Ntests + 1):
    pyplot.figure(i)
    pyplot.subplot(231)
    pyplot.pcolor(d["f" + str(i)][1:-1, 0, :].T)
    pyplot.colorbar()
    pyplot.title("f" + str(i))
    pyplot.xlabel("X")
    pyplot.ylabel("Z")
    pyplot.subplot(232)
    pyplot.pcolor(d["sol" + str(i)][1:-1, 0, :].T)
    pyplot.colorbar()
    pyplot.title("sol" + str(i))
    pyplot.xlabel("X")
    pyplot.ylabel("Z")
    pyplot.subplot(233)
    pyplot.pcolor(d["absolute_error" + str(i)][1:-1, 0, :].T)
    pyplot.colorbar()
    pyplot.title("absolute_error" + str(i))
    pyplot.xlabel("X")
    pyplot.ylabel("Z")
    pyplot.subplot(234)
    pyplot.pcolor(d["b" + str(i)][2:-2, 0, :].T)
    pyplot.colorbar()
    pyplot.title("b" + str(i))
    pyplot.xlabel("X")
    pyplot.ylabel("Z")
    pyplot.subplot(235)
    pyplot.pcolor(d["bcheck" + str(i)][2:-2, 0, :].T)
    pyplot.colorbar()
    pyplot.title("bcheck" + str(i))
    pyplot.xlabel("X")
    pyplot.ylabel("Z")
    pyplot.subplot(236)
    pyplot.pcolor(d["b" + str(i)][2:-2, 0, :].T - d["bcheck" + str(i)][2:-2, 0, :].T)
    pyplot.colorbar()
    pyplot.title("b" + str(i) + "-bcheck" + str(i))
    pyplot.xlabel("X")
    pyplot.ylabel("Z")

pyplot.show()

exit(0)
