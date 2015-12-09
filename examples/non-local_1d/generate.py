#!/usr/bin/env python

#
# Generate an input mesh
#

from __future__ import division
from past.utils import old_div
from boututils.datafile import DataFile # Wrapper around NetCDF4 libraries
from math import pow
from sys import argv

length = 80. # Length of the domain in m

nx = 5   # Minimum is 5: 2 boundary, one evolved
if len(argv)>1:
  ny = int(argv[1])  # Minimum 5. Should be divisible by number of processors (so powers of 2 nice)
else:
  ny = 256  # Minimum 5. Should be divisible by number of processors (so powers of 2 nice)
#dy = [[1.]*ny]*nx # distance between points in y, in m/g22/lengthunit
g22 = [[pow(old_div(float(ny-1),length),2)]*ny]*nx
g_22 = [[pow(old_div(length,float(ny-1)),2)]*ny]*nx
ixseps1 = -1
ixseps2 = 0

f = DataFile()
f.open("conduct_grid.nc", create=True)

f.write("nx", nx)
f.write("ny", ny)
#f.write("dy", dy)
f.write("g22",g22)
f.write("g_22", g_22)
f.write("ixseps1", ixseps1)
f.write("ixseps2", ixseps2)

f.close()
