#!/usr/bin/env python

"""Generate an input mesh"""

# Wrapper around NetCDF4 libraries
from boututils import DataFile

def generate_grid(nx = 20,\
                  ny = 10,\
                  nz = 8,\
                  inp_path  = "data"
                  file_name = "3D_diffusion_grid"):

    # Read geom:Lx and geom:Ly (used to calculate dx and dy)


    # Create the grid file
    f = DataFile()
    f.open(file_name + ".nc", create=True)

    # Write the dimensions to the grid file
    f.write("nx", nx)
    f.write("ny", ny)
    f.write("nz", nz)

    # Write the grid sizes to the grid file
    f.write("dx", dx)
    f.write("dy", dy)

    # Write the lengths to the grid file
    f.write("Lx", Lx)
    f.write("Ly", Ly)

    # Close the grid file
    f.close()
