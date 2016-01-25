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

def pol_slice(var3d, gridfile, n=1, zangle=0.0, withRepeat=False, nyInterp=None):
    """ data2d = pol_slice(data3d, 'gridfile', n=1, zangle=0.0) 
    By default data3d will not include the repeated z point. If it does
    contain this point one should set withRepeat=True to get the correct result.
    Set nyInterp to the number of y (theta) points you wish to use in the final result.
    """
    n = int(n)
    zangle = float(zangle)

    s = np.shape(var3d)
    if len(s) != 3:
        print("ERROR: pol_slice expects a 3D variable")
        return None

    nx, ny, nz = s

    #Fix to account for repeat point which is now usually dropped
    if not withRepeat:
        nz = nz+1 #Increase the nz value to include the repeated point, which is not present in var3d

    dz = 2.*np.pi / float(n * (nz-1))

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

        if zShift != None:
            print("Using qinty as toroidal shift angle")
        else:
            zShift = gf.read("zShift")
            if zShift != None:
                print("Using zShift as toroidal shift angle")
            else:
                print("ERROR: Neither qinty nor zShift found")
                return None

        gf.close()
    except:
        print("ERROR: pol_slice couldn't read grid file")
        return None

    #Decide if we've asked to do interpolation
    doInterp = False
    if nyInterp is not None:
        if ny != nyInterp:
            doInterp = True

    #Setup interpolation if requested
    if doInterp:
        try:
            from numpy import linspace, meshgrid, append, newaxis
        except ImportError:
            print("ERROR: NumPy module not available")
            raise
        try:
            from scipy.ndimage import map_coordinates
        except ImportError:
            print("ERROR: SciPy module not available")
            raise

        #print("Interpolating in y")
        #Create a copy of first point at end of array, we'd then drop the last y point
        #in the interpolated data. This is in effort to encourage periodicity at the inboard
        #side, though it doesn't seem to work (possibly because it doesn't do the z-shifting
        #that is really required to achieve this).
#        varTmp = append(var3d,var3d[:,0,newaxis,:],axis=1)

        varTmp = var3d
        #These are the index space co-ordinates of the output arrays
        xOut = linspace(0,nx-1,nx)
        yOut = linspace(0,ny-1,nyInterp)
        zOut = linspace(0,nz-1,nz)

        #Use meshgrid to create 3 length N arrays (N=total number of points)
        xx,yy,zz = meshgrid(xOut,yOut,zOut,indexing='ij')
        #Interpolate to output positions and make the correct shape
        var3d = map_coordinates(input=varTmp,coordinates=[xx,yy,zz],cval=-999)

        #As above
        xx,yy = meshgrid(xOut,yOut,indexing='ij')
        zShift = map_coordinates(zShift,[xx,yy],cval=-999)

        #Update shape
        ny = nyInterp

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
    # for x in np.arange(nx):
    #     for y in np.arange(ny):
    #         var2d[x,y] = 0.5*p[x,y]*(p[x,y]-1.0) * var3d[x,y,zm[x,y]] + \
    #                      (1.0 - p[x,y]*p[x,y])   * var3d[x,y,z0[x,y]] + \
    #                      0.5*p[x,y]*(p[x,y]+1.0) * var3d[x,y,zp[x,y]]

    # var2dList = [x.T for x in var3d.T]
    # var2d = 0.5*p*(p-1.0) * np.choose(zm,var2dList) + \
    #         (1.0 - p*p)   * np.choose(z0,var2dList) + \
    #         0.5*p*(p+1.0) * np.choose(zp,var2dList) 

    #For some reason numpy imposes a limit of 32 entries to choose
    #so if nz>32 we have to use a different approach
    if nz >= 32:
        #Untested
        # v3m = np.array([var3d.T[zm.T[I]][I] for I in ndi.ndindex(zm.T.shape)]).T
        # v30 = np.array([var3d.T[z0.T[I]][I] for I in ndi.ndindex(z0.T.shape)]).T
        # v3p = np.array([var3d.T[zp.T[I]][I] for I in ndi.ndindex(zp.T.shape)]).T
        # var2d = 0.5*p*(p-1.0) * v3m + \
        #         (1.0 - p*p) * v30 + \
        #         0.5*p*(p+1.0) * v3p
        for x in np.arange(nx):
            for y in np.arange(ny):
                var2d[x,y] = 0.5*p[x,y]*(p[x,y]-1.0) * var3d[x,y,zm[x,y]] + \
                             (1.0 - p[x,y]*p[x,y])   * var3d[x,y,z0[x,y]] + \
                             0.5*p[x,y]*(p[x,y]+1.0) * var3d[x,y,zp[x,y]]
    else:
        var2d = 0.5*p*(p-1.0) * np.choose(zm.T,var3d.T).T + \
                (1.0 - p*p) * np.choose(z0.T,var3d.T).T + \
                0.5*p*(p+1.0) * np.choose(zp.T,var3d.T).T


    # var2d = 0.5*p*(p-1.0) * var3d[:,:,zm] + \
    #                      (1.0 - p*p)   * var3d[:,:,z0] + \
    #                      0.5*p*(p+1.0) * var3d[:,:,zp]


    return var2d
