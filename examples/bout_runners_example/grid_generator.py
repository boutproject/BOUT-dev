#!/usr/bin/env python

"""Generate an input mesh"""

# Wrapper around NetCDF4 libraries
from boututils import DataFile

def grid_generator(nx = 5,\
                   ny = 64,\
                   nz = 1,\
                   file_name = "conduct_grid"):

    # Create the grid file
    f = DataFile()
    f.open(file_name + ".nc", create=True)

    # Write the dimensions to the grid file
    f.write("nx", nx)
    f.write("ny", ny)

    # Close the grid file
    f.close()
