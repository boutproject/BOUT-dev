#!/usr/bin/env python

#
# Generate an input mesh
#

from boututils.datafile import DataFile # Wrapper around NetCDF4 libraries

nx = 5   # Minimum is 5: 2 boundary, one evolved
ny = 32  # Minimum 5. Should be divisible by number of processors (so powers of 2 nice)
dy = 1. # distance between points in y, in m/g22/lengthunit
ixseps1 = -1
ixseps2 = -1

f = DataFile()
f.open("test-staggered.nc", create=True)

f.write("nx", nx)
f.write("ny", ny)
f.write("dy", dy)
f.write("ixseps1", ixseps1)
f.write("ixseps2", ixseps2)

f.close()
