#!/usr/bin/env python3

import boutpp

# requires boutpp
# requires not make

import numpy as np
import itertools
from sys import exit

errorlist = []

boutpp.init("-d data -q -q -q".split(" "))


def runtests(functions, derivatives, directions, stag, msg):
    global errorlist
    for direction in directions:
        direction, fac, guards, diff_func = direction
        locations = ["centre"]
        if stag:
            locations.append(direction.lower() + "low")
        for funcs, derivative, inloc, outloc in itertools.product(
            functions, derivatives, locations, locations
        ):
            infunc, outfunc = funcs
            order, diff = derivative

            errors = []
            for nz in nzs:
                boutpp.setOption(
                    "meshD:nD".replace("D", direction),
                    "%d" % (nz + (2 * guards if direction == "x" else 0)),
                    force=True,
                )
                boutpp.setOption(
                    "meshD:dD".replace(
                        "D",
                        direction,
                    ),
                    "2*pi/(%d)" % (nz),
                    force=True,
                )
                dirnfac = direction + "*" + fac
                mesh = boutpp.Mesh(section="mesh" + direction)
                f = boutpp.create3D(infunc.replace("%s", dirnfac), mesh, outloc=inloc)
                sim = diff_func(f, method=diff, outloc=outloc)
                if sim.getLocation() != outloc:
                    cent = ["CENTRE", "CENTER"]
                    if outloc in cent and sim.getLocation() in cent:
                        pass
                    else:
                        errorlist.append(
                            "Location does not match - expected %s but got %s"
                            % (outloc, sim.getLocation())
                        )
                ana = boutpp.create3D(
                    outfunc.replace("%s", dirnfac), mesh, outloc=outloc
                )
                err = sim - ana
                err = err.getAll().flatten()
                if guards:
                    err = err[guards:-guards]
                err = np.max(np.abs(err))
                errors.append(err)
            errc = np.log(errors[-2] / errors[-1])
            difc = np.log(nzs[-1] / nzs[-2])
            conv = errc / difc
            if order - 0.1 < conv < order + 0.1:
                pass
            else:
                info = "%s - %s - %s - %s -> %s " % (
                    infunc,
                    diff,
                    direction,
                    inloc,
                    outloc,
                )
                error = "%s: %s is not working. Expected %f got %f" % (
                    msg,
                    info,
                    order,
                    conv,
                )
                errorlist.append(error)
                if doPlot:
                    from matplotlib import pyplot as plt

                    plt.plot((ana).getAll().flatten())
                    plt.plot((sim).getAll().flatten())
                    plt.show()


mmax = 7
start = 6
doPlot = False
nzs = np.logspace(start, mmax, num=mmax - start + 1, base=2)

functions = [["sin(%s)", "cos(%s)"], ["cos(%s)", "-sin(%s)"]]

derivatives = [
    [2, "C2"],
    [4, "C4"],
    # [2,"W2"] ,
    # [3,"W3"] ,
    # [2,"S2"] ,
]

directions = [
    ["x", "2*pi", 2, boutpp.DDX],
    ["y", "1", 2, boutpp.DDY],
    #    ["z","1"   ,0 ,boutpp.DDZ]
]

runtests(functions, derivatives, directions, stag=False, msg="DD")

derivatives = [
    [2, "C2"],
    [4, "C4"],
]

runtests(functions, derivatives, directions, stag=True, msg="DD")


functions = [["sin(%s)", "-sin(%s)"], ["cos(%s)", "-cos(%s)"]]

derivatives = [[2, "C2"], [4, "C4"]]

directions = [
    ["x", "2*pi", 2, boutpp.D2DX2],
    ["y", "1", 2, boutpp.D2DY2],
    #    ["z","1"   ,0 ,boutpp.D2DZ2]
]

runtests(functions, derivatives, directions, False, "D2D2")

derivatives = [
    [2, "C2"],
]

runtests(functions, derivatives, directions, True, "D2D2")

boutpp.finalise()

if errorlist:
    for error in errorlist:
        print(error)
    exit(1)
else:
    print("Pass")
    exit(0)
