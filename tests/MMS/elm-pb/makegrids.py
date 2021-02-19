#!/usr/bin/env python3
#
# Creates the grid files and input options for
# a set of simulations. These can then be submitted
# to a batch queue and then processed with runtest
# (setting "running = False")
#

from __future__ import print_function
from __future__ import division

nxlist = [4, 8, 16, 32, 64, 128]

from boututils.datafile import DataFile
from boututils.run_wrapper import shell

from os.path import isfile

from boutdata.mms import SimpleTokamak, x, y
from sympy import sin, cos

from math import pi

shape = SimpleTokamak()

## Add equilibrium profiles
MU0 = 4.0e-7 * pi

J0 = 1 - x - sin(x * pi) ** 2 * cos(y)  # Parallel current
P0 = 2 + cos(x * pi)  # Pressure pedestal
bxcvz = -((1.0 / shape.Rxy) ** 2) * cos(y)  # Curvature

# Normalisation
Lbar = 1.0
Bbar = 1.0
J0 = -J0 * shape.Bxy / (MU0 * Lbar)  # Turn into A/m^2
P0 = P0 * Bbar ** 2 / (2.0 * MU0)  # Pascals

shape.add(P0, "pressure")
shape.add(J0, "Jpar0")
shape.add(bxcvz, "bxcvz")

for nx in nxlist:
    # Generate a new mesh file

    filename = "grid%d.nc" % nx

    if isfile(filename):
        print("Grid file '%s' already exists" % filename)
    else:
        print("Creating grid file '%s'" % filename)
        f = DataFile(filename, create=True)
        shape.write(nx, nx, f)
        f.close()

    # Generate BOUT.inp file

    directory = "grid%d" % nx
    shell("mkdir " + directory)
    shell("cp data/BOUT.inp " + directory)
    shell("sed -i 's/MZ = 17/MZ = %d/g' %s/BOUT.inp" % (nx, directory))
    shell(
        'sed -i \'s/grid = "grid16.nc"/grid = "%s"/g\' %s/BOUT.inp'
        % (filename, directory)
    )
