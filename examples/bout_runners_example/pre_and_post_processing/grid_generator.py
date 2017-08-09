#!/usr/bin/env python

"""Generate an input mesh"""

from boututils.datafile import DataFile
from boututils.options import BOUTOptions
import numpy as np
import os

# Define pi, in case it is found in the BOUT.inp file
pi = np.pi

def generate_grid(nx        = 20  ,\
                  ny        = 10  ,\
                  nz        = 8   ,\
                  Lx        = None,\
                  Ly        = None,\
                  MXG       = None,\
                  inp_path  = None,\
                  file_name = None):
    """Generates a grid file based on the input, and on what is being
    found in the BOUT.inp file"""

    # Calculate dx and dy
    #{{{ Use BOUT.inp if Lx, Ly or MXG is not found
    if (Lx == None) or (Ly == None) or (MXG == None):
        # Make a BOUTOption object (see documentation of BOUTOption.py)
        myOpts = BOUTOptions(inp_path)

        # Appendable list
        warnings = []

        # If Lx is not given
        if Lx == None:
            # Read 'Lx' from the 'geom' section
            Lx = myOpts.geom['Lx']
            # Lx is now a string, we let python evaluate the string
            Lx = eval(Lx)
            # Append a warning
            warnings.append('Lx')
        # If Ly is not given
        if Ly == None:
            Ly = myOpts.geom['Ly']
            Ly = eval(Ly)
            warnings.append('Ly')
        # If MXG is not given
        if MXG == None:
            MXG = myOpts.root['MXG']
            MXG = eval(MXG)
            warnings.append('MXG')

        # Print warnings
        for warning in warnings:
            print("\n"*2 + "!"*137)
            print("WARNING!!! " + warning + " not given in generate_grid")
            print("Will use value from BOUT.inp to calculate gridspacing, but"+\
                  " note that this would be inconsistent if '" + warning +\
                  "' is given in a bout_runner object")
            print("!"*137 + "\n"*2)
    #}}}

    # Calculate dx and dy
    internal_x_points = nx - 2*MXG
    internal_y_points = ny
    # Since the grid points lay half between the grid
    # (There is one less line segment than points, and there is one more
    #  "internal" in the grid due to 2 half grid points out to the
    #  boundary)
    dx = Lx / (internal_x_points)
    dy = Ly / (internal_y_points)

    # Set ixseps
    ixseps1 = -1
    ixseps2 = -2

    # Check that folder exists
    grid_folder = os.path.split(file_name)[0]
    if grid_folder != "":
        if not os.path.exists(grid_folder):
                os.makedirs(grid_folder)

    # Write the grid file
    with DataFile(file_name, write=True, create=True) as grid:
        # Write the dimensions to the grid file
        grid.write("nx", nx)
        grid.write("ny", ny)
        grid.write("nz", nz)

        # Write the grid sizes to the grid file
        grid.write("dx", dx)
        grid.write("dy", dy)

        # Write the lengths to the grid file
        grid.write("Lx", Lx)
        grid.write("Ly", Ly)

        # Write the ixseps
        grid.write("ixseps1", ixseps1)
        grid.write("ixseps2", ixseps2)
