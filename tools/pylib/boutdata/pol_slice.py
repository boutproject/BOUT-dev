from __future__ import print_function
from __future__ import division

# Takes a 3D variable, and returns a 2D slice at fixed toroidal angle
#
# N sets the number of times the data must be repeated for a full
# torus, e.g. n=2 is half a torus
# zangle gives the (real) toroidal angle of the result

try:
    import numpy as np
except ImportError:
    print("ERROR: NumPy module not available")
    raise

try:
    from boututils.datafile import DataFile
except ImportError:
    print("ERROR: boututils.DataFile not available")
    print("=> Set $PYTHONPATH variable to include BOUT++ pylib")
    raise SystemExit

def pol_slice(var3d, gridfile, n=1, zangle=0.0):
    """ data2d = pol_slice(data3d, 'gridfile', n=1, zangle=0.0) """
    n = int(n)
    zangle = float(zangle)

    s = np.shape(var3d)
    if len(s) != 3:
        print("ERROR: pol_slice expects a 3D variable")
        return None

    nx, ny, nz = s

    dz = 2.*np.pi / float(n * nz)

    try:
        # Open the grid file
        gf = DataFile(gridfile)

        # Check the grid size is correct
        if gf.read("nx") != nx:
            print("ERROR: Grid X size is different to the variable")
            return None
        if gf.read("ny") != ny:
            print("ERROR: Grid Y size is different to the variable")
            return None

        # Get the toroidal shift
        zShift = gf.read("qinty")

        if zShift is not None:
            print("Using qinty as toroidal shift angle")
        else:
            zShift = gf.read("zShift")
            if zShift is not None:
                print("Using zShift as toroidal shift angle")
            else:
                print("ERROR: Neither qinty nor zShift found")
                return None

        gf.close()
    except:
        print("ERROR: pol_slice couldn't read grid file")
        return None

    var2d = np.zeros([nx, ny])

    ######################################
    # Perform 2D slice
    zind = (zangle - zShift) / dz
    z0f = np.floor(zind)
    z0 = z0f.astype(int)
    p = zind - z0f

    # Make z0 between 0 and (nz-2)
    z0 = ((z0 % (nz-1)) + (nz-1)) % (nz-1)

    # Get z+ and z-
    zp = (z0 + 1) % (nz-1)
    zm = (z0 - 1 + (nz-1)) % (nz-1)

    # There may be some more cunning way to do this indexing
    for x in np.arange(nx):
        for y in np.arange(ny):
            var2d[x,y] = 0.5*p[x,y]*(p[x,y]-1.0) * var3d[x,y,zm[x,y]] + \
                         (1.0 - p[x,y]*p[x,y])   * var3d[x,y,z0[x,y]] + \
                         0.5*p[x,y]*(p[x,y]+1.0) * var3d[x,y,zp[x,y]]

    return var2d
