#!/usr/bin/env python

"""Generate an input mesh"""

from boututils import DataFile
from boututils.options import BOUTOptions
import numpy as np

# Define pi, in case it is found in the BOUT.inp file
pi = np.pi

def generate_grid(nx = 20,\
                  ny = 10,\
                  nz = 8,\
                  Lx = None,\
                  Ly = None,\
                  MXG = None,\
                  inp_path  = None,\
                  file_name = None):
    """Generates a grid file based on the input, and on what is being
    found in the BOUT.inp file"""

    # Calculate dx and dy
    # Open BOUT.inp if Lx, Ly or MXG is not found
    if (Lx == None) or (Ly == None) or (MXG == None):
        # Make a BOUTOption object (see documentation of BOUTOption.py)
        myOpts = BOUTOptions(inp_path)

    # If Lx is not given
    if Lx == None:
        print("WARNING!!! Lx not given in generate_grid")
        print("Will use value from BOUT.inp to calculate dx, but"+\
              " note that this will be inconsistent if 'Lx' is given"+\
              " in a bout_runner object")
        # Read 'Lx' from the 'geom' section
        Lx = myOpts.geom['Lx']
        # Lx is now a string, we let python evaluate the string
        Lx = eval(Lx)
    # If Ly is not given
    if Ly == None:
        print("WARNING!!! Ly not given in generate_grid")
        print("Will use value from BOUT.inp to calculate dx, but"+\
              " note that this will be inconsistent if 'Ly' is given"+\
              " in a bout_runner object")
        # Read 'Ly' from the 'geom' section
        Ly = myOpts.geom['Ly']
        # Ly is now a string, we let python evaluate the string
        Ly = eval(Ly)
    if MXG == None:
        print("WARNING!!! MXG not given in generate_grid")
        print("Will use value from BOUT.inp to calculate dx, but"+\
              " note that this will be inconsistent if 'MXG' is given"+\
              " in a bout_runner object")
        # Read 'MXG' from the 'geom' section
        MXG = myOpts.root['MXG']
        # MXG is now a string, we let python evaluate the string
        MXG = eval(MXG)

    # Calculate dx and dy
    internal_x_points = nx - 2*MXG
    internal_y_points = ny
    # Since the grid points lay half between the grid
    # (There is one less line segment than points, and there is one more
    #  "internal" in the grid due to 2 half grid points out to the
    #  boundary)
    dx = Lx / (internal_x_points)
    dy = Ly / (internal_y_points)

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
