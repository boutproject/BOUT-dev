#!/usr/bin/env python

#
# Generate an input mesh
#

from boututils.datafile import DataFile # Wrapper around NetCDF4 libraries

nx = 5   # Minimum is 5: 2 boundary, one evolved
ny = 64  # Minimum 5. Should be divisible by number of processors (so powers of 2 nice)

f = DataFile()
f.open("conduct_grid.nc", create=True)

f.write("nx", nx)
f.write("ny", ny)

f.close()
