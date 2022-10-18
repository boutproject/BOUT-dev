#!/usr/bin/env python3
import boutpp
import numpy as np

# requires boutpp
# requires not make

errorlist = ""
boutpp.init("-d data -q -q -q")  # +" -f BOUT.settings")

shapes = []
errors = []
mmax = 7
start = 6
doPlot = False
nzs = np.logspace(start, mmax, num=mmax - start + 1, base=2)
for nz in nzs:
    boutpp.setOption("meshz:nz", "%d" % nz, force=True)
    boutpp.setOption("meshz:dz", "2*pi/%d" % nz, force=True)
    mesh = boutpp.Mesh(section="meshz")
    f = boutpp.create3D("sin(z)", mesh)
    sim = boutpp.DDZ(f)
    ana = boutpp.create3D("cos(z)", mesh)
    err = sim - ana
    err = np.max(np.abs(err.getAll()))
    errors.append(err)

if doPlot:
    from matplotlib import pyplot as plt

    plt.plot(1 / np.array(nzs), errors, "-o")
    plt.show()

errc = np.log(errors[-2] / errors[-1])
difc = np.log(nzs[-1] / nzs[-2])
conv = errc / difc

if 1.9 < conv < 2.1:
    print("The convergence is: %f" % conv)
else:
    errorlist += "DDZ is not working\n"

if errorlist != "":
    print("Found errors:\n%s" % errorlist)
    exit(1)
