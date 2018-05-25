"""Routines for manipulating grid files

"""
from __future__ import print_function

from boututils.datafile import DataFile

from numpy import ndarray, zeros, concatenate


def slice(infile, outfile, region=None, xind=None, yind=None):
    """Copy an X-Y slice from one DataFile to another

    TODO: rename to not clobber builtin `slice`
    TODO: better regions?

    Parameters
    ----------
    infile : str
        Name of DataFile to read slice from
    outfile : str
        Name of DataFile to write slice to. File will be created, and
        will be overwritten if it already exists
    region : {0, 1, 2, 3, 4, 5, None}, optional
        Copy a whole region. The available regions are:
            - 0: Lower inner leg
            - 1: Inner core
            - 2: Upper inner leg
            - 3: Upper outer leg
            - 4: Outer core
            - 5: Lower outer leg
    xind, yind : (int, int), optional
        Index ranges for x and y. Range includes first point, but not
        last point

    """

    # Open input and output files
    indf = DataFile(infile)
    outdf = DataFile(outfile, create=True)
    
    nx = indf["nx"][0]
    ny = indf["ny"][0]

    if region:
        # Select a region of the mesh
        
        xind = [0, nx]
        if region == 0:
            # Lower inner leg
            yind = [0, indf["jyseps1_1"][0]+1]
        elif region == 1:
            # Inner core
            yind = [indf["jyseps1_1"][0]+1, indf["jyseps2_1"][0]+1]
        elif region == 2:
            # Upper inner leg
            yind = [indf["jyseps2_1"][0]+1, indf["ny_inner"][0]]
        elif region == 3:
            # Upper outer leg
            yind = [indf["ny_inner"][0], indf["jyseps1_2"][0]+1]
        elif region == 4:
            # Outer core
            yind = [indf["jyseps1_2"][0]+1, indf["jyseps2_2"][0]+1]
        else:
            # Lower outer leg
            yind = [indf["jyseps2_2"][0]+1, ny]
    else:
        # Use indices
        if not xind:
            xind = [0, nx]
        if not yind:
            yind = [0, ny]
    
    print("Indices: [%d:%d, %d:%d]" % (xind[0], xind[1], yind[0], yind[1]))
    # List of variables requiring special handling
    special = ["nx", "ny", "ny_inner",
               "ixseps1", "ixseps2", 
               "jyseps1_1", "jyseps1_2", "jyseps2_1", "jyseps2_2",
               "ShiftAngle"]
    
    outdf["nx"] = xind[1] - xind[0]
    outdf["ny"] = yind[1] - yind[0]    
    outdf["ny_inner"] = indf["ny_inner"][0] - yind[0]

    outdf["ixseps1"] = indf["ixseps1"][0]
    outdf["ixseps2"] = indf["ixseps2"][0]
    
    outdf["jyseps1_1"] = indf["jyseps1_1"][0] - yind[0]
    outdf["jyseps2_1"] = indf["jyseps2_1"][0] - yind[0]
    outdf["jyseps1_2"] = indf["jyseps1_2"][0] - yind[0]
    outdf["jyseps2_2"] = indf["jyseps2_2"][0] - yind[0]

    outdf["ShiftAngle"] = indf["ShiftAngle"][xind[0]:xind[1]]
    
    # Loop over all variables
    for v in list(indf.keys()):
        if v in special:
            continue # Skip these variables
        
        ndims = indf.ndims(v)
        if ndims == 0:
            # Copy scalars
            print("Copying variable: " + v)
            outdf[v] = indf[v][0]
        elif ndims == 2:
            # Assume [x,y]
            print("Slicing variable: " + v);
            outdf[v] = indf[v][xind[0]:xind[1], yind[0]:yind[1]]
        else:
            # Skip
            print("Skipping variable: " + v)

    indf.close()
    outdf.close()


def rotate(gridfile, yshift, output=None):
    """Shifts a grid file by the specified number of points in y

    This moves the branch cut around, and can be used to change the
    limiter location

    Parameters
    ----------
    gridfile : str
        Name of DataFile to rotate
    yshift : int
        Number of points in y to shift by
    output : str, optional
        Name of DataFile to write to. If None, will write to a new
        file with the same name as `gridfile` + '_rot'

    """

    if output is None:
        output = gridfile + "_rot"

    print("Rotating grid file '%s' -> '%s'" % (gridfile, output))
        
    # Open input grid file
    with DataFile(gridfile) as d:
        # Open output file
        with DataFile(output, write=True, create=True) as out:
            # Loop over variables
            for varname in d.list():
                # Number of dimensions
                ndims = d.ndims(varname)
                
                if ndims == 2:
                    print("Shifting '%s' (x,y)" % (varname,))
                    # 2D, assume X-Y

                    var = d[varname] # Read
                    ny = var.shape[1]
                    
                    # Make sure yshift is positive and in range
                    yshift = ((yshift % ny) + ny) % ny
                    
                    newvar = ndarray(var.shape)

                    # Rotate
                    newvar[:,0:(ny-yshift)] = var[:,yshift:ny]
                    newvar[:,(ny-yshift):] = var[:,:yshift]

                    # Write to output
                    #out[varname] = newvar # Write
                    out.write(varname, newvar)
                elif ndims == 3:
                    print("Shifting '%s' (x,y,z)" % (varname,))
                    # 3D, assume X-Y-Z
                    
                    var = d[varname] # Read
                    ny = var.shape[1]

                    # Make sure yshift is positive and in range
                    yshift = ((yshift % ny) + ny) % ny
                    
                    newvar = ndarray(var.shape)

                    newvar[:,0:(ny-yshift),:] = var[:,yshift:ny,:]
                    newvar[:,(ny-yshift):,:] = var[:,:yshift,:]

                    # Write to output
                    out.write(varname, newvar)
                else:
                    # Just copy
                    print("Copying '%s' (%d dimensions)" % (varname, ndims))
                    out.write(varname, d[varname])

import matplotlib.pyplot as plt        
from numpy import linspace, amin, amax


def gridcontourf(grid, data2d, nlevel=31, show=True,
                 mind=None, maxd=None, symmetric=False,
                 cmap=None, ax=None,
                 xlabel="Major radius [m]", ylabel="Height [m]",
                 separatrix=False):
    """Plots a 2D contour plot, taking into account branch cuts
    (X-points).

    TODO: move into a plotting module

    Parameters
    ----------
    grid : DataFile
        A DataFile object
    data2d : array_like
        A 2D (x,y) NumPy array of data to plot
    nlevel : int, optional
        Number of levels in the contour plot
    show : bool, optional
        If True, will immediately show the plot
    mind : float, optional
        Minimum data level
    maxd : float, optional
        Maximum data level
    symmetric : bool, optional
        Make mind, maxd symmetric about zero
    cmap : Colormap, optional
        A matplotlib colormap to use. If None, use the current default
    ax : Axes, optional
        A matplotlib axes instance to plot to. If None, create a new
        figure and axes, and plot to that
    xlabel, ylabel : str, optional
        Labels for the x/y axes
    separatrix : bool, optional
        Add separatrix

    Returns
    -------
    con
        The contourf instance

    Examples
    --------

    To put a plot into an axis with a color bar:

    >>> fig, axis = plt.subplots()
    >>> c = gridcontourf(grid, data, show=False, ax=axis)
    >>> fig.colorbar(c, ax=axis)
    >>> plt.show()

    """

    if cmap is None:
        cmap = plt.cm.get_cmap("YlOrRd")

    if len(data2d.shape) != 2:
        raise ValueError("data2d must be 2D (x,y)")
    
    j11 = grid["jyseps1_1"]
    j12 = grid["jyseps1_2"]
    j21 = grid["jyseps2_1"]
    j22 = grid["jyseps2_2"]
    ix1 = grid["ixseps1"]
    ix2 = grid["ixseps2"]
    try:
      nin = grid["ny_inner"]
    except:
      nin = j12
     
    nx  = grid["nx"]
    ny  = grid["ny"]
    
    if (data2d.shape[0] != nx) or (data2d.shape[1] != ny):
        raise ValueError("data2d has wrong size: (%d,%d), expected (%d,%d)" % (data2d.shape[0], data2d.shape[1], nx, ny))

    if hasattr(j11, "__len__"):
        # Arrays rather than scalars
        try:
          j11 = j11[0]
          j12 = j12[0]
          j21 = j21[0]
          j22 = j22[0]
          ix1 = ix1[0]
          ix2 = ix2[0]
          nin = nin[0]
          nx  = nx[0]
          ny  = ny[0]
        except:
          pass

    R = grid["Rxy"]
    Z = grid["Zxy"]

    if data2d.shape != (nx, ny):
        raise ValueError("Dimensions do not match")

    add_colorbar = False
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        add_colorbar = True
    
    if mind is None:
      mind = amin(data2d)
    if maxd is None:
      maxd = amax(data2d)    

    if symmetric:
        # Make mind, maxd symmetric about zero
        maxd = max([maxd, abs(mind)])
        mind = -maxd

    levels = linspace(mind, maxd, nlevel, endpoint=True)

    ystart = 0  # Y index to start the next section
    if j11 >= 0:
        # plot lower inner leg
        ax.contourf(R[:,ystart:(j11+1)], Z[:,ystart:(j11+1)], data2d[:,ystart:(j11+1)], levels,cmap=cmap)
    
        yind = [j11, j22+1]
        ax.contourf(R[:ix1, yind].transpose(), Z[:ix1, yind].transpose(), data2d[:ix1, yind].transpose(), levels,cmap=cmap)
        
        ax.contourf(R[ix1:,j11:(j11+2)], Z[ix1:,j11:(j11+2)], data2d[ix1:,j11:(j11+2)], levels,cmap=cmap)
        ystart = j11+1
    
        yind = [j22, j11+1]
        ax.contourf(R[:ix1, yind].transpose(), Z[:ix1, yind].transpose(), data2d[:ix1, yind].transpose(), levels, cmap=cmap)
        
    # Inner SOL
    con = ax.contourf(R[:,ystart:(j21+1)], Z[:,ystart:(j21+1)], data2d[:,ystart:(j21+1)], levels, cmap=cmap)
    ystart = j21+1
    
    if j12 > j21: 
        # Contains upper PF region
        
        # Inner leg
        ax.contourf(R[ix1:,j21:(j21+2)], Z[ix1:,j21:(j21+2)], data2d[ix1:,j21:(j21+2)], levels, cmap=cmap)
        ax.contourf(R[:,ystart:nin], Z[:,ystart:nin], data2d[:,ystart:nin], levels, cmap=cmap)
        
        # Outer leg
        ax.contourf(R[:,nin:(j12+1)], Z[:,nin:(j12+1)], data2d[:,nin:(j12+1)], levels, cmap=cmap)
        ax.contourf(R[ix1:,j12:(j12+2)], Z[ix1:,j12:(j12+2)], data2d[ix1:,j12:(j12+2)], levels, cmap=cmap)
        ystart = j12+1
        
        yind = [j21, j12+1]
        ax.contourf(R[:ix1, yind].transpose(), Z[:ix1, yind].transpose(), data2d[:ix1, yind].transpose(), levels, cmap=cmap)

        yind = [j21+1, j12]
        ax.contourf(R[:ix1, yind].transpose(), Z[:ix1, yind].transpose(), data2d[:ix1, yind].transpose(), levels, cmap=cmap)
    else:
        ystart -= 1
    # Outer SOL
    ax.contourf(R[:,ystart:(j22+1)], Z[:,ystart:(j22+1)], data2d[:,ystart:(j22+1)], levels, cmap=cmap)
    
    ystart = j22+1

    if j22+1 < ny:
        # Outer leg
        ax.contourf(R[ix1:,j22:(j22+2)], Z[ix1:,j22:(j22+2)], data2d[ix1:,j22:(j22+2)], levels, cmap=cmap)
        ax.contourf(R[:,ystart:ny], Z[:,ystart:ny], data2d[:,ystart:ny], levels, cmap=cmap)

        # X-point
        Rx = [ [R[ix1-1,j11], R[ix1,j11], R[ix1,j11+1], R[ix1-1,j11+1]],
               [R[ix1-1,j22+1], R[ix1,j22+1], R[ix1,j22], R[ix1-1,j22]] ]

        
        Zx = [ [Z[ix1-1,j11], Z[ix1,j11], Z[ix1,j11+1], Z[ix1-1,j11+1]],
               [Z[ix1-1,j22+1], Z[ix1,j22+1], Z[ix1,j22], Z[ix1-1,j22]] ]
        Dx = [ [data2d[ix1-1,j11], data2d[ix1,j11], data2d[ix1,j11+1], data2d[ix1-1,j11+1]],
               [data2d[ix1-1,j22+1], data2d[ix1,j22+1], data2d[ix1,j22], data2d[ix1-1,j22]] ]
        ax.contourf(Rx, Zx, Dx, levels, cmap=cmap)

    if add_colorbar:
        fig.colorbar(con)
        
    ax.set_aspect("equal")
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    if separatrix: 
        # Plot separatrix
        
        # Lower X-point location
        Rx = 0.125*(R[ix1-1,j11] + R[ix1,j11] + R[ix1,j11+1] + R[ix1-1,j11+1]
                    + R[ix1-1,j22+1] + R[ix1,j22+1] + R[ix1,j22] + R[ix1-1,j22])
        Zx = 0.125*(Z[ix1-1,j11] + Z[ix1,j11] + Z[ix1,j11+1] + Z[ix1-1,j11+1]
                    + Z[ix1-1,j22+1] + Z[ix1,j22+1] + Z[ix1,j22] + Z[ix1-1,j22])
        # Lower inner leg
        ax.plot( concatenate( (0.5*(R[ix1-1,0:(j11+1)] + R[ix1,0:(j11+1)]), [Rx]) ), concatenate( (0.5*(Z[ix1-1,0:(j11+1)] + Z[ix1,0:(j11+1)]), [Zx]) ), 'k-')
        # Lower outer leg
        ax.plot( concatenate( ([Rx],0.5*(R[ix1-1,(j22+1):] + R[ix1,(j22+1):])) ), concatenate( ([Zx], 0.5*(Z[ix1-1,(j22+1):] + Z[ix1,(j22+1):])) ), 'k-')
        # Core
        
        ax.plot( concatenate( ([Rx], 0.5*(R[ix1-1,(j11+1):(j21+1)] + R[ix1,(j11+1):(j21+1)]), 0.5*(R[ix1-1,(j12+1):(j22+1)] + R[ix1,(j12+1):(j22+1)]), [Rx]) ),
                 concatenate( ([Zx], 0.5*(Z[ix1-1,(j11+1):(j21+1)] + Z[ix1,(j11+1):(j21+1)]), 0.5*(Z[ix1-1,(j12+1):(j22+1)] + Z[ix1,(j12+1):(j22+1)]), [Zx]) ), 'k-')
    if show:
      plt.show()
      
    return con


def bout2sonnet(grdname, outf):
    """Creates a Sonnet format grid from a BOUT++ grid.

    NOTE: Branch cuts are not yet supported

    Parameters
    ----------
    grdname : str
        Filename of BOUT++ grid file
    outf : File
        The file-like object to write to

    Examples
    --------

    >>> with open("output.sonnet", "w") as f:
    ...     bout2sonnet("BOUT.grd.nc", f)

    """

    with DataFile(grdname) as g:
        Rxy = g["Rxy"]
        Zxy = g["Zxy"]
        Bpxy = g["Bpxy"]
        Btxy = g["Btxy"]
        Bxy = g["Bxy"]

    # Now iterate over cells in the order Eirene expects

    nx, ny = Rxy.shape

    # Extrapolate values in Y
    R = zeros([nx,ny+2])
    Z = zeros([nx,ny+2])

    R[:,1:-1] = Rxy
    Z[:,1:-1] = Zxy

    R[:,0] = 2.*R[:,1] - R[:,2]
    Z[:,0] = 2.*Z[:,1] - Z[:,2]

    R[:,-1] = 2.*R[:,-2] - R[:,-3]
    Z[:,-1] = 2.*Z[:,-2] - Z[:,-3]
    
    element = 1  # Element number

    outf.write("BOUT++: "+grdname+"\n\n")
    
    outf.write("=====================================\n")
    
    for i in range(2, nx-2):
        # Loop in X, excluding guard cells
        for j in range(1,ny+1):
            # Loop in Y. Guard cells not in grid file

            # Lower left (low Y, low X)
            ll = ( 0.25*(R[i-1,j-1] + R[i-1,j] + R[i,j-1] + R[i,j]),
                   0.25*(Z[i-1,j-1] + Z[i-1,j] + Z[i,j-1] + Z[i,j]) )

            # Lower right (low Y, upper X)
            lr = ( 0.25*(R[i+1,j-1] + R[i+1,j] + R[i,j-1] + R[i,j]),
                   0.25*(Z[i+1,j-1] + Z[i+1,j] + Z[i,j-1] + Z[i,j]) )

            # Upper left (upper Y, lower X)
            ul = ( 0.25*(R[i-1,j+1] + R[i-1,j] + R[i,j+1] + R[i,j]),
                   0.25*(Z[i-1,j+1] + Z[i-1,j] + Z[i,j+1] + Z[i,j]) )
        
            # Upper right (upper Y, upper X)
            ur = ( 0.25*(R[i+1,j+1] + R[i+1,j] + R[i,j+1] + R[i,j]),
                   0.25*(Z[i+1,j+1] + Z[i+1,j] + Z[i,j+1] + Z[i,j]) )
            
            # Element number
            outf.write("     ELEMENT   %d = ( %d, %d): (%e, %e) (%e, %e)\n" % (
                element,
                j-1, i-2,
                ll[0], ll[1],
                ul[0], ul[1]))

            # Ratio Bt / Bp at cell centre. Note j-1 because
            # Bpxy and Btxy have not had extra points added
            outf.write("     FIELD RATIO  = %e  (%e, %e)\n" % (Bpxy[i,j-1] / Btxy[i,j-1], R[i,j], Z[i,j]) )
            
            outf.write("                         (%e, %e) (%e, %e)\n" % (
                lr[0], lr[1],
                ur[0], ur[1]))

            if (i == nx-3) and (j == ny+1):
                # Last element
                outf.write("=====================================\n")
            else:
                outf.write("-------------------------------------\n")

            element += 1
