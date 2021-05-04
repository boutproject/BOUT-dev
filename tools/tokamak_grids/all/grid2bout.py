"""
@brief Convert an equilibrium set of flux surfaces into a BOUT++ grid


"""
from __future__ import print_function

from numpy import max

from boututils.datafile import DataFile


def grid2bout(input, output="bout.grd.nc"):
    # List variables needed and dimensions
    dimensions = {
        "nx": 0,
        "ny": 0,
        "psi": 1,
        "mu0p": 1,
        "mu0pprime": 1,
        "f": 1,
        "ffprime": 1,
        "qsafe": 1,
        "Rxy": 2,
        "Zxy": 2,
    }

    # Open the input file
    infile = Datafile(input)

    # Create an empty dictionary for the variables
    vars = {}

    # Check the required variables
    for k, v in list(dimensions.items()):
        nd = infile.ndims(k)
        if nd == None:
            print("ERROR: Variable missing: " + k)
            return False
        if v != nd:
            print("ERROR: Variable '" + k + "' has wrong number of dimensions")
            return False

        # Read the variables
        vars[k] = f.read(k)

    # Read the size of the grid
    nx = vars["nx"]
    ny = vars["ny"]

    # Location of grid points
    Rxy = vars["Rxy"]
    Zxy = vars["Zxy"]

    # Check if the last point is the same as the first
    dl = (Rxy[:, ny - 1] - Rxy[:, 0]) ** 2 + (Zxy[:, ny - 1] - Zxy[:, 0]) ** 2
    if dl.max() < 1.0e-4:
        print("**Last poloidal point duplicates the first. Removing...")
