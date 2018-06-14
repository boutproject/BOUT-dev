from __future__ import print_function
from __future__ import division

from boututils.datafile import DataFile
import numpy as np
from scipy.ndimage import map_coordinates


def pol_slice(var3d, gridfile, n=1, zangle=0.0, nyInterp=None):
    """Takes a 3D variable, and returns a 2D slice at fixed toroidal angle

    Parameters
    ----------
    var3d : array_like
        The input array. Should be 3D
    gridfile : str
        The gridfile containing the coordinate system to used
    n : int, optional
        The number of times the data must be repeated for a full torus,
        e.g. n=2 is half a torus
    zangle : float, optional
        The (real) toroidal angle of the result
    nyInterp : int, optional
        The number of y (theta) points to use in the final result.

    Returns
    -------
    array
        A 2D-slice of var3d interpolated at a fixed toroidal angle
    """
    n = int(n)
    zangle = float(zangle)

    s = np.shape(var3d)
    if len(s) != 3:
        raise ValueError("pol_slice expects a 3D variable (got {} dimensions)"
                         .format(len(s)))

    nx, ny, nz = s

    # Open the grid file
    with DataFile(gridfile) as gf:
        # Check the grid size is correct
        grid_nx = gf.read("nx")
        if grid_nx != nx:
            raise ValueError("Grid X size ({}) is different to the variable ({})"
                             .format(grid_nx, nx))
        grid_ny = gf.read("ny")
        if grid_ny != ny:
            raise ValueError("Grid Y size ({}) is different to the variable ({})"
                             .format(grid_ny, ny))

        # Get the toroidal shift
        zShift = gf.read("qinty")

        if zShift is not None:
            print("Using qinty as toroidal shift angle")
        else:
            zShift = gf.read("zShift")
            if zShift is not None:
                print("Using zShift as toroidal shift angle")
            else:
                raise ValueError("Neither qinty nor zShift found")

    # Decide if we've asked to do interpolation
    if nyInterp is not None and nyInterp != ny:
        varTmp = var3d

        # Interpolate to output positions and make the correct shape
        # np.mgrid gives us an array of indices
        var3d = map_coordinates(varTmp, coordinates=np.mgrid[0:nx, 0:nyInterp, 0:nz],
                                cval=-999)
        zShift = map_coordinates(zShift, np.mgrid[0:nx, 0:nyInterp], cval=-999)

        # Update shape
        ny = nyInterp

    var2d = np.zeros([nx, ny])

    ######################################
    # Perform 2D slice
    dz = 2.*np.pi / float(n * nz)
    zind = (zangle - zShift) / dz
    z0f = np.floor(zind)
    z0 = z0f.astype(int)
    p = zind - z0f

    # Make z0 between 0 and (nz-2)
    z0 = ((z0 % (nz-1)) + (nz-1)) % (nz-1)

    # Get z+ and z-
    zp = (z0 + 1) % (nz-1)
    zm = (z0 - 1 + (nz-1)) % (nz-1)

    # For some reason numpy imposes a limit of 32 entries to choose
    # so if nz>32 we have to use a different approach. This limit may change with numpy version
    if nz >= 32:
        for x in np.arange(nx):
            for y in np.arange(ny):
                var2d[x, y] = (0.5*p[x, y]*(p[x, y]-1.0) * var3d[x, y, zm[x, y]] +
                               (1.0 - p[x, y]*p[x, y]) * var3d[x, y, z0[x, y]] +
                               0.5*p[x, y]*(p[x, y]+1.0) * var3d[x, y, zp[x, y]])
    else:
        var2d = (0.5*p*(p-1.0) * np.choose(zm.T, var3d.T).T +
                 (1.0 - p*p) * np.choose(z0.T, var3d.T).T +
                 0.5*p*(p+1.0) * np.choose(zp.T, var3d.T).T)

    return var2d
