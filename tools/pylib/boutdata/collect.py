# Requires:
#  - boututils
#  - NumPy

try:
    from boututils import DataFile
except ImportError:
    print "ERROR: boututils.DataFile couldn't be loaded"
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

def collect(varname, xind=None, yind=None, zind=None, tind=None, path=".",yguards=False, info=True,prefix="BOUT.dmp"):
    """Collect a variable from a set of BOUT++ outputs.
    
    data = collect(name)
    
    name   Name of the variable (string)
    
    Optional arguments:

    xind = [min,max]   Range of X indices to collect
    yind = [min,max]   Range of Y indices to collect
    zind = [min,max]   Range of Z indices to collect
    tind = [min,max]   Range of T indices to collect
    
    path    = "."          Path to data files
    prefix  = "BOUT.dmp"   File prefix
    yguards = False        Collect Y boundary guard cells?
    info    = True         Print information about collect?
    """
    
    # Search for BOUT++ dump files in NetCDF format
    file_list = glob.glob(os.path.join(path, prefix+".nc"))
    if file_list != []:
        print "Single (parallel) data file"
        f = DataFile(file_list[0]) # Open the file
        
        data = f.read(varname)
        return data
    
    file_list = glob.glob(os.path.join(path, prefix+"*.nc"))
    file_list.sort()
    if file_list == []:
        print "ERROR: No data files found"
        return None
    nfiles = len(file_list)
    #print "Number of files: " + str(nfiles)
    
    # Read data from the first file
    f = DataFile(file_list[0])
    
    #print "File format    : " + f.file_format
    try:
        dimens = f.dimensions(varname)
        ndims = len(dimens)
    except KeyError:
        print "ERROR: Variable '"+varname+"' not found"
        return None

    if ndims < 2:
        # Just read from file
        data = f.read(varname)
        f.close()
        return data

    if ndims > 4:
        print "ERROR: Too many dimensions"
        raise CollectError

    mxsub = f.read("MXSUB")
    mysub = f.read("MYSUB")
    mz    = f.read("MZ")
    myg   = f.read("MYG")
    t_array = f.read("t_array")
    nt = len(t_array)
    
    if info:
        print "mxsub = %d mysub = %d mz = %d\n" % (mxsub, mysub, mz)

    # Get the version of BOUT++ (should be > 0.6 for NetCDF anyway)
    try:
        v = f.read("BOUT_VERSION")

        # 2D decomposition
        nxpe = f.read("NXPE")
        mxg  = f.read("MXG")
        nype = f.read("NYPE")
        npe = nxpe * nype
        
        if info:
            print "nxpe = %d, nype = %d, npe = %d\n" % (nxpe, nype, npe)
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

    if yguards:
        ny = mysub * nype + 2*myg
    else:
        ny = mysub * nype
    
    f.close();

    # Check ranges
    
    def check_range(r, low, up, name="range"):
        r2 = r
        if r != None:
            try:
                n = len(r2)
            except:
                # No len attribute, so probably a single number
                r2 = [r2,r2]
            if (len(r2) < 1) or (len(r2) > 2):
                print "WARNING: "+name+" must be [min, max]"
                r2 = None
            else:
                if len(r2) == 1:
                    r2 = [r2,r2]
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
    ddims = map(lambda d: sizes[d], dimens)
    
    # Create the data array
    data = np.zeros(ddims)
    
    for i in range(npe):
        # Get X and Y processor indices
        pe_yind = int(i / nxpe)
        pe_xind = i % nxpe

        # Get local ranges
        if yguards:
            ymin = yind[0] - pe_yind*mysub
            ymax = yind[1] - pe_yind*mysub
        else:
            ymin = yind[0] - pe_yind*mysub + myg
            ymax = yind[1] - pe_yind*mysub + myg
        
        xmin = xind[0] - pe_xind*mxsub
        xmax = xind[1] - pe_xind*mxsub
        
        inrange = True

        if yguards:
            # Check lower y boundary
            if pe_yind == 0:
                # Keeping inner boundary
                if ymax < 0: inrange = False
                if ymin < 0: ymin = 0
            else:
                if ymax < myg: inrange = False
                if ymin < myg: ymin = myg

            # Upper y boundary
            if pe_yind == (nype - 1):
                # Keeping outer boundary
                if ymin >= (mysub + 2*myg): inrange = False
                if ymax > (mysub + 2*myg - 1): ymax = (mysub + 2*myg - 1)
            else:
                if ymin >= (mysub + myg): inrange = False
                if ymax >= (mysub + myg): ymax = (mysub+myg-1)
            
        else:
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

        if yguards:
            ygmin = ymin + pe_yind * mysub
            ygmax = ymax + pe_yind * mysub

        else:
            ygmin = ymin + pe_yind * mysub - myg
            ygmax = ymax + pe_yind * mysub - myg


        if not inrange:
            continue # Don't need this file
        
        filename = os.path.join(path, prefix+"." + str(i) + ".nc")
        if info:
            sys.stdout.write("\rReading from " + filename + ": [" + \
                                 str(xmin) + "-" + str(xmax) + "][" + \
                                 str(ymin) + "-" + str(ymax) + "] -> [" + \
                                 str(xgmin) + "-" + str(xgmax) + "][" + \
                                 str(ygmin) + "-" + str(ygmax) + "]")

        f = DataFile(filename)

        if ndims == 4:
            d = f.read(varname, ranges=[tind[0],tind[1]+1,
                                        xmin, xmax+1, 
                                        ymin, ymax+1, 
                                        zind[0],zind[1]+1])
            data[:, (xgmin-xind[0]):(xgmin-xind[0]+nx_loc), (ygmin-yind[0]):(ygmin-yind[0]+ny_loc), :] = d
        elif ndims == 3:
            # Could be xyz or txy
            
            if dimens[2] == 'z': # xyz
                d = f.read(varname, ranges=[xmin, xmax+1, 
                                            ymin, ymax+1, 
                                            zind[0],zind[1]+1])
                data[(xgmin-xind[0]):(xgmin-xind[0]+nx_loc), (ygmin-yind[0]):(ygmin-yind[0]+ny_loc), :] = d
            else: # txy
                d = f.read(varname, ranges=[tind[0],tind[1]+1,
                                            xmin, xmax+1, 
                                            ymin, ymax+1])
                data[:, (xgmin-xind[0]):(xgmin-xind[0]+nx_loc), (ygmin-yind[0]):(ygmin-yind[0]+ny_loc)] = d
        elif ndims == 2:
            # xy
            d = f.read(varname, ranges=[xmin, xmax+1, 
                                        ymin, ymax+1])
            data[(xgmin-xind[0]):(xgmin-xind[0]+nx_loc), (ygmin-yind[0]):(ygmin-yind[0]+ny_loc)] = d

        f.close()
    
    # Finished looping over all files
    if info:
        sys.stdout.write("\n")
    return data
