"""
Contour calculation routines

https://web.archive.org/web/20140901225541/https://members.bellatlantic.net/~vze2vrva/thesis.html"""
from __future__ import print_function
from __future__ import division
from past.utils import old_div

import numpy as np


def contour(f, level):
    """Return a list of contours matching the given level"""

    if len(f.shape) != 2:
        print("Contour only works on 2D data")
        return None
    nx,ny = f.shape

    # Go through each cell edge and mark which ones contain
    # a level crossing. Approximating function as
    # f = axy + bx + cy + d
    # Hence linear interpolation along edges.

    edgecross = {} # Dictionary: (cell number, edge number) to crossing location

    for i in np.arange(nx-1):
        for j in np.arange(ny-1):
            # Lower-left corner of cell is (i,j)
            if (np.max(f[i:(i+2),j:(j+2)]) < level) or (np.min(f[i:(i+2),j:(j+2)]) > level):
                # not in the correct range - skip
                continue

            # Check each edge
            ncross = 0
            def location(a, b):
                if (a > level) ^ (a > level):
                    # One of the corners is > level, and the other is <= level
                    ncross += 1
                    # Find location
                    return old_div((level - a), (b - a))
                else:
                    return None

            loc = [
                location(f[i,j], f[i+1,j]),
                location(f[i+1,j], f[i+1,j+1]),
                location(f[i+1,j+1], f[i,j+1]),
                location(f[i,j+1], f[i,j])]

            if ncross != 0: # Only put into dictionary if has a crossing
                cellnr = (ny-1)*i + j # The cell number
                edgecross[cellnr] = [loc,ncross] # Tack ncross onto the end

    # Process crossings into contour lines

    while True:
        # Start from an arbitrary location and follow until
        # it goes out of the domain or closes
        try:
            startcell, cross = edgecross.popitem()
        except KeyError:
            # No keys left so finished
            break

        def follow():
            return

        # Follow

    return

def find_opoints(var2d):
    """Find O-points in psi i.e. local minima/maxima"""
    return

def find_xpoints(var2d):
    """Find X-points in psi i.e. inflection points"""
    return

