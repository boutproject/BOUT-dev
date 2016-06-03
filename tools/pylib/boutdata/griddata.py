from __future__ import print_function
# Routines for manipulating grid files

try:
    from boututils import DataFile
except ImportError:
    print("ERROR: restart module needs DataFile")
    raise

from numpy import zeros

def slice(infile, outfile, region = None, xind=None, yind=None):
    """
    infile     - Name of the input file
    outfile    - Name of the output file to be created
    region     - The region of the mesh to slice. 
                 0 = Lower inner leg
                 1 = Inner core
                 2 = Upper inner leg
                 3 = Upper outer leg
                 4 = Outer core
                 5 = Lower outer leg
    xind, yind - index ranges. Range includes first point, 
                 but not last point
    
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


def bout2sonnet(grdname, outf):
    """
    Creates a Sonnet format grid from a BOUT++ grid.
    NOTE: Branch cuts are not yet supported 
    
    Inputs
    ------
    
    grdname - Filename of BOUT++ grid file
    
    outf    - The file-like object to write to
    
    Example
    -------
    
    with open("output.sonnet", "w") as f:
      bout2sonnet("BOUT.grd.nc", f)
    
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

