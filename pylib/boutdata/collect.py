# Requires:
#  - netcdf4-python (HDF5, NetCDF-4)
#  - NumPy

try:
    from netCDF4 import Dataset
except ImportError:
    print "ERROR: netcdf4-python module not found"
    raise

try:
    import os
    import sys
    import glob
except ImportError:
    print "ERROR: os, sys or glob modules not available"
    raise

try:
    import numpy as np
except ImportError:
    print "ERROR: NumPy module not available"
    raise

print "    data = collect('variable', path='.')"

def collect(varname, xind=None, yind=None, zind=None, tind=None, path="."):
    """Collect a variable from a set of BOUT++ outputs."""
    
    # Little helper function to read a variable
    def read_var(file, name):
        var = file.variables[name]
        return var[:]
    
    # Search for BOUT++ dump files in NetCDF format
    file_list = glob.glob(os.path.join(path, "BOUT.dmp.*.nc"))
    if file_list == []:
        print "ERROR: No data files found"
        return None
    nfiles = len(file_list)
    print "Number of files: " + str(nfiles)
    
    # Read data from the first file
    f = Dataset(file_list[0], "r")
    print "File format    : " + f.file_format
    try:
        v = f.variables[varname]
    except KeyError:
        print "ERROR: Variable '"+varname+"' not found"
        return None
    dims = v.dimensions
    ndims = len(dims)

    if ndims == 0:
        # Just read from this file and return
        data = v.getValue()
        f.close()
        return data[0]
    elif ndims == 1:
        data = v[:]
        f.close()
        return data
    elif ndims > 4:
        print "ERROR: Too many dimensions"
        f.close()
        raise CollectError

    mxsub = read_var(f, "MXSUB")[0]
    mysub = read_var(f, "MYSUB")[0]
    mz    = read_var(f, "MZ")[0]
    myg   = read_var(f, "MYG")[0]
    t_array = read_var(f, "t_array")
    nt = len(t_array)
    print "Time-points    : " + str(nt)
    
    # Get the version of BOUT++ (should be > 0.6 for NetCDF anyway)
    try:
        v = f.variables["BOUT_VERSION"]
        print "BOUT++ version : " + str(v.getValue()[0])

        # 2D decomposition
        nxpe = read_var(f, "NXPE")[0]
        mxg  = read_var(f, "MXG")[0]
        nype = read_var(f, "NYPE")[0]
        npe = nxpe * nype

        if npe < nfiles:
            print "WARNING: More files than expected (" + str(npe) + ")"
        elif npe > nfiles:
            print "WARNING: Some files missing. Expected " + str(npe)

        nx = nxpe * mxsub + 2*mxg
    except KeyError:
        print "BOUT++ version : Pre-0.2"
        # Assume number of files is correct
        # No decomposition in X
        nx = mxsub
        mxg = 0
        nxpe = 1
        nype = nfiles
    
    ny = mysub * nype
    
    f.close()
    
    # Check ranges
    
    def check_range(r, low, up, name="range"):
        r2 = r
        if r != None:
            if (len(r) < 1) or (len(r) > 2):
                print "WARNING: "+name+" must be [min, max]"
                r2 = None
            else:
                if len(r) == 1:
                    r = [r,r]
                if r2[0] < low:
                    r2[0] = low
                if r2[0] > up:
                    r2[0] = up
                if r2[1] < 0:
                    r2[1] = 0
                if r2[1] > up:
                    r2[1] = up
                if r2[0] > r2[1]:
                    tmp = r2[0]
                    r2[0] = r2[1]
                    r2[1] = tmp
        else:
            r2 = [low, up]
        return r2
    
    xind = check_range(xind, 0, nx-1, "xind")
    yind = check_range(yind, 0, ny-1, "yind")
    zind = check_range(zind, 0, mz-2, "zind")
    tind = check_range(tind, 0, nt-1, "tind")
    
    
    xsize = xind[1] - xind[0] + 1
    ysize = yind[1] - yind[0] + 1
    zsize = zind[1] - zind[0] + 1
    tsize = tind[1] - tind[0] + 1
    
    # Map between dimension names and output size
    sizes = {'x':xsize, 'y':ysize, 'z':zsize, 't':tsize}

    # Create a list with size of each dimension
    ddims = map(lambda d: sizes[d], dims)
    
    # Create the data array
    data = np.zeros(ddims)
    
    for i in range(nfiles):
        # Get X and Y processor indices
        pe_yind = int(i / nxpe)
        pe_xind = i % nxpe

        # Get local ranges
        ymin = yind[0] - pe_yind*mysub + myg
        ymax = yind[1] - pe_yind*mysub + myg

        xmin = xind[0] - pe_xind*mxsub
        xmax = xind[1] - pe_xind*mxsub
        
        inrange = True

        if (ymin >= (mysub + myg)) or (ymax < myg):
            inrange = False # Y out of range

        if ymin < myg:
            ymin = myg
        if ymax >= mysub+myg:
            ymax = myg + mysub - 1

        # Check lower x boundary
        if pe_xind == 0:
            # Keeping inner boundary
            if xmax < 0: inrange = False
            if xmin < 0: xmin = 0
        else:
            if xmax < mxg: inrange = False
            if xmin < mxg: xmin = mxg
        
        # Upper x boundary
        if pe_xind == (nxpe - 1):
            # Keeping outer boundary
            if xmin >= (mxsub + 2*mxg): inrange = False
            if xmax > (mxsub + 2*mxg - 1): xmax = (mxsub + 2*mxg - 1)
        else:
            if xmin >= (mxsub + mxg): inrange = False
            if xmax >= (mxsub + mxg): xmax = (mxsub+mxg-1)

        # Number of local values
        nx_loc = xmax - xmin + 1
        ny_loc = ymax - ymin + 1

        # Calculate global indices
        xgmin = xmin + pe_xind * mxsub
        xgmax = xmax + pe_xind * mxsub

        ygmin = ymin + pe_yind * mysub - myg
        ygmax = ymax + pe_yind * mysub - myg

        if not inrange:
            continue # Don't need this file
        
        filename = os.path.join(path, "BOUT.dmp." + str(i) + ".nc")
        sys.stdout.write("\rReading from " + filename + ": [" + \
                         str(xmin) + "-" + str(xmax) + "][" + \
                         str(ymin) + "-" + str(ymax) + "] -> [" + \
                         str(xgmin) + "-" + str(xgmax) + "][" + \
                         str(ygmin) + "-" + str(ygmax) + "]")

        f = Dataset(filename, "r")
        var = f.variables[varname]

        if ndims == 4:
            d = var[tind[0]:(tind[1]+1), xmin:(xmax+1), ymin:(ymax+1), zind[0]:(zind[1]+1)]
            data[:, (xgmin-xind[0]):(xgmin-xind[0]+nx_loc), (ygmin-yind[0]):(ygmin-yind[0]+ny_loc), :] = d
        elif ndims == 3:
            # Could be xyz or txy
            
            if dims[3] == 'z': # xyz
                d = var[xmin:(xmax+1), ymin:(ymax+1), zind[0]:(zind[1]+1)]
                data[(xgmin-xind[0]):(xgmin-xind[0]+nx_loc), (ygmin-yind[0]):(ygmin-yind[0]+ny_loc), :] = d
            else: # txy
                d = var[tind[0]:(tind[1]+1), xmin:(xmax+1), ymin:(ymax+1)]
                data[:, (xgmin-xind[0]):(xgmin-xind[0]+nx_loc), (ygmin-yind[0]):(ygmin-yind[0]+ny_loc)] = d
        elif ndims == 2:
            # xy
            d = var[xmin:(xmax+1), ymin:(ymax+1)]
            data[(xgmin-xind[0]):(xgmin-xind[0]+nx_loc), (ygmin-yind[0]):(ygmin-yind[0]+ny_loc)] = d
    
    # Finished looping over all files
    sys.stdout.write("\n")
    return data
