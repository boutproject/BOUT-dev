#!/usr/bin/env python

"""Generate an input mesh"""

# Wrapper around NetCDF4 libraries
from boututils import DataFile

def generate_grid(nx = 20,\
                  ny = 10,\
                  nz = 8,\
                  file_name = "3D_diffusion_grid"):

    # Create the grid file
    f = DataFile()
    f.open(file_name + ".nc", create=True)

    # Write the dimensions to the grid file
    f.write("nx", nx)
    f.write("ny", ny)

    # Close the grid file
    f.close()
