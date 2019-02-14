"""Flux surface generator for tokamak grid files

"""
from __future__ import print_function

import numpy as np


def gen_surface(grid):
    """Generator for iterating over flux surfaces

    Parameters
    ----------
    grid : DataFile
        An input grid file to read to find flux surfaces

    Yields
    ------
    tuple : (int, list of int, bool)
        A tuple containing the x index, list of y indices and whether
        this flux surface is periodic

    """
    # Read the grid data
    nx = grid.read("nx")
    ny = grid.read("ny")

    npol = grid.read("npol")
    if npol is None:
        # Domains not stored in file (BOUT style input)
        ixseps1 = grid.read("ixseps1")
        ixseps2 = grid.read("ixseps2")
        jyseps1_1 = grid.read("jyseps1_1")
        jyseps1_2 = grid.read("jyseps1_2")
        jyseps2_1 = grid.read("jyseps2_1")
        jyseps2_2 = grid.read("jyseps2_2")

        if ixseps1 == ixseps2:
            # Single null
            ndomains = 3
        else:
            # Double null
            ndomains = 6

        yup_xsplit = np.zeros(ndomains)
        ydown_xsplit = np.zeros(ndomains)
        yup_xin = np.zeros(ndomains)
        yup_xout = np.zeros(ndomains)
        ydown_xin = np.zeros(ndomains)
        ydown_xout = np.zeros(ndomains)

        ystart = np.zeros(ndomains+1)
        ystart[ndomains] = ny

        # Inner lower leg
        ydown_xsplit[0] = -1
        ydown_xout[0] = -1
        yup_xsplit[0] = ixseps1
        yup_xin[0] = ndomains-1 # Outer lower leg
        yup_xout[0] = 1

        # Outer lower leg
        ydown_xsplit[ndomains-1] = ixseps1
        ydown_xin[ndomains-1] = 0
        ydown_xout[ndomains-1] = ndomains-2
        yup_xsplit[ndomains-1] = -1
        yup_xout[ndomains-1] = -1
        ystart[ndomains-1] = jyseps2_2+1

        if ixseps1 == ixseps2:
            # Single null

            ydown_xsplit[1] = ixseps1
            ydown_xin[1] = 1
            ydown_xout[1] = 0
            yup_xsplit[1] = ixseps1
            yup_xin[1] = 1
            yup_xout[1] = 2
            ystart[1] = jyseps1_1+1
        else:
            # Double null
            raise RuntimeError("SORRY - NO DOUBLE NULL YET")
    else:
        # Use domains stored in the file
        ndomains = npol.size # Number of domains
        yup_xsplit   = grid.read("yup_xsplit")
        ydown_xsplit = grid.read("ydown_xsplit")
        yup_xin    = grid.read("yup_xin")
        yup_xout   = grid.read("yup_xout")
        ydown_xin  = grid.read("ydown_xin")
        ydown_xout = grid.read("ydown_xout")

        # Calculate starting positions
        ystart = np.zeros(ndomains+1)
        for i in np.arange(1,ndomains):
            ystart[i] = ystart[i-1] + npol[i-1]
        ystart[ndomains] = ny

    # Record whether a domain has been visited
    visited = np.zeros(ndomains)

    x = 0      # X index
    while True:
        yinds = None # Y indices result

        # Find a domain which hasn't been visited
        domain = None
        for i in np.arange(ndomains):
            if visited[i] == 0:
                domain = i
                break

        if domain is None:
            # All domains visited
            x = x + 1 # Go to next x surface
            visited = np.zeros(ndomains) # Clear the visited array
            if x == nx:
                break # Finished
            continue

        # Follow surface back until it hits a boundary
        while True:
            if x < ydown_xsplit[domain]:
                d = ydown_xin[domain]
            else:
                d = ydown_xout[domain]
            if d < 0:
                break # Hit boundary
            domain = d # Keep going

        # Starting from domain, follow surface

        periodic = False
        while domain >= 0:
            if visited[domain] == 1:
                # Already visited domain -> periodic
                periodic = True
                break;
            # Get range of y indices in this domain
            yi = np.arange(ystart[domain], ystart[domain+1])
            if yinds is None:
                yinds = yi
            else:
                yinds = np.concatenate((yinds, yi))
            # mark this domain as visited
            visited[domain] = 1
            # Get next domain
            if x < yup_xsplit[domain]:
                domain = yup_xin[domain]
            else:
                domain = yup_xout[domain]

        # Finished this surface
        yield x, yinds, periodic
