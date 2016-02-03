from __future__ import print_function
# Routines for manipulating grid files

try:
    from boututils.datafile import DataFile
except ImportError:
    print("ERROR: restart module needs DataFile")
    raise

from numpy import ndarray


def slice(infile, outfile, region = None, xind=None, yind=None):
    """
    xind, yind - index ranges. Range includes first point, but not last point
    
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
    """
    Shifts a grid file by the specified number of points in y
    
    This moves the branch cut around, and can be used
    to change the limiter location
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
        
