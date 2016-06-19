from __future__ import print_function
from __future__ import division
try:
    from builtins import str
except:
    print("Warning: No str in builtins")

try:
    from builtins import range
except:
    print("Warning: No range in builtins")

# Requires:
#  - boututils
#  - NumPy

try:
    from boututils.datafile import DataFile
except ImportError:
    print("ERROR: boututils.DataFile couldn't be loaded")
    raise

try:
    import os
    import sys
    import glob
except ImportError:
    print("ERROR: os, sys or glob modules not available")
    raise

try:
    import numpy as np
except ImportError:
    print("ERROR: NumPy module not available")
    raise

def findVar(varname, varlist):
    """
    Find variable name in a list

    First does case insensitive comparison, then
    checks for abbreviations.

    Returns the matched string, or raises a ValueError

    """
    # Try a variation on the case
    v = [name for name in varlist if name.lower() == varname.lower()]
    if len(v) == 1:
        # Found case match
        print("Variable '%s' not found. Using '%s' instead" % (varname, v[0]))
        return v[0]
    elif len(v) > 1:
        print("Variable '"+varname+"' not found, and is ambiguous. Could be one of: "+str(v))
        raise ValueError("Variable '"+varname+"' not found")

    # None found. Check if it's an abbreviation
    v = [name for name in varlist if name[:len(varname)].lower() == varname.lower()]
    if len(v) == 1:
        print("Variable '%s' not found. Using '%s' instead" % (varname, v[0]))
        return v[0]

    if len(v) > 1:
        print("Variable '"+varname+"' not found, and is ambiguous. Could be one of: "+str(v))
    raise ValueError("Variable '"+varname+"' not found")

def collect(varname, xind=None, yind=None, zind=None, tind=None, path=".",yguards=False, xguards=True, info=True,prefix="BOUT.dmp",strict=False):
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
    xguards = True         Collect X boundary guard cells?
                           (Set to True to be consistent with the
                           definition of nx)
    info    = True         Print information about collect?
    strict  = False        Fail if the exact variable name is not found?
    """

    # Search for BOUT++ dump files in NetCDF format
    file_list_nc = glob.glob(os.path.join(path, prefix+".nc"))
    file_list_h5 = glob.glob(os.path.join(path, prefix+".hdf5"))
    if file_list_nc != [] and file_list_h5 != []:
        raise IOError("Error: Both NetCDF and HDF5 files are present: do not know which to read.")
    elif file_list_h5 != []:
        suffix = ".hdf5"
        file_list = file_list_h5
    else:
        suffix = ".nc"
        file_list = file_list_nc
    if file_list != []:
        print("Single (parallel) data file")
        f = DataFile(file_list[0]) # Open the file

        data = f.read(varname)
        return data

    file_list_nc = glob.glob(os.path.join(path, prefix+".*nc"))
    file_list_h5 = glob.glob(os.path.join(path, prefix+".*hdf5"))
    if file_list_nc != [] and file_list_h5 != []:
        raise IOError("Error: Both NetCDF and HDF5 files are present: do not know which to read.")
    elif file_list_h5 != []:
        suffix = ".hdf5"
        file_list = file_list_h5
    else:
        suffix = ".nc"
        file_list = file_list_nc

    file_list.sort()
    if file_list == []:
        raise IOError("ERROR: No data files found")

    nfiles = len(file_list)

    # Read data from the first file
    f = DataFile(file_list[0])

    try:
        dimens = f.dimensions(varname)
        #ndims = len(dimens)
        ndims = f.ndims(varname)
    except:
        if strict:
            raise
        else:
            # Find the variable
            varname = findVar(varname, f.list())

            dimens = f.dimensions(varname)
            #ndims = len(dimens)
            ndims = f.ndims(varname)

    # ndims is 0 for reals, and 1 for f.ex. t_array
    if ndims < 2:
        # Just read from file
        if varname != 't_array':
            data = f.read(varname)
        elif (varname == 't_array') and (tind is None):
            data = f.read(varname)
        elif (varname == 't_array') and (tind is not None):
            data = f.read(varname, ranges=[tind[0],tind[1]+1])
        f.close()
        return data

    if ndims > 4:
        raise ValueError("ERROR: Too many dimensions")

    mxsub = f.read("MXSUB")
    if mxsub is None:
        raise ValueError("Missing MXSUB variable")
    mysub = f.read("MYSUB")
    mz    = f.read("MZ")
    myg   = f.read("MYG")
    t_array = f.read("t_array")
    if t_array is None:
        nt = 1
        t_array = np.zeros(1)
    else:
        nt = len(t_array)

    if info:
        print("mxsub = %d mysub = %d mz = %d\n" % (mxsub, mysub, mz))

    # Get the version of BOUT++ (should be > 0.6 for NetCDF anyway)
    try:
        v = f.read("BOUT_VERSION")

        # 2D decomposition
        nxpe = f.read("NXPE")
        mxg  = f.read("MXG")
        nype = f.read("NYPE")
        npe = nxpe * nype

        if info:
            print("nxpe = %d, nype = %d, npe = %d\n" % (nxpe, nype, npe))
            if npe < nfiles:
                print("WARNING: More files than expected (" + str(npe) + ")")
            elif npe > nfiles:
                print("WARNING: Some files missing. Expected " + str(npe))

        if xguards:
            nx = nxpe * mxsub + 2*mxg
        else:
            nx = nxpe * mxsub
    except KeyError:
        print("BOUT++ version : Pre-0.2")
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
                print("WARNING: "+name+" must be [min, max]")
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
    ddims = [sizes[d] for d in dimens]

    # Create the data array
    data = np.zeros(ddims)

    for i in range(npe):
        # Get X and Y processor indices
        pe_yind = int(i/nxpe)
        pe_xind = i % nxpe

        inrange = True

        if yguards:
            # Get local ranges
            ymin = yind[0] - pe_yind*mysub
            ymax = yind[1] - pe_yind*mysub

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

            # Calculate global indices
            ygmin = ymin + pe_yind * mysub
            ygmax = ymax + pe_yind * mysub

        else:
            # Get local ranges
            ymin = yind[0] - pe_yind*mysub + myg
            ymax = yind[1] - pe_yind*mysub + myg

            if (ymin >= (mysub + myg)) or (ymax < myg):
                inrange = False # Y out of range

            if ymin < myg:
                ymin = myg
            if ymax >= mysub+myg:
                ymax = myg + mysub - 1

            # Calculate global indices
            ygmin = ymin + pe_yind * mysub - myg
            ygmax = ymax + pe_yind * mysub - myg

        if xguards:
            # Get local ranges
            xmin = xind[0] - pe_xind*mxsub
            xmax = xind[1] - pe_xind*mxsub

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

            # Calculate global indices
            xgmin = xmin + pe_xind * mxsub
            xgmax = xmax + pe_xind * mxsub

        else:
            # Get local ranges
            xmin = xind[0] - pe_xind*mxsub + mxg
            xmax = xind[1] - pe_xind*mxsub + mxg

            if (xmin >= (mxsub + mxg)) or (xmax < mxg):
                inrange = False # X out of range

            if xmin < mxg:
                xmin = mxg
            if xmax >= mxsub+mxg:
                xmax = mxg + mxsub - 1

            # Calculate global indices
            xgmin = xmin + pe_xind * mxsub - mxg
            xgmax = xmax + pe_xind * mxsub - mxg


        # Number of local values
        nx_loc = xmax - xmin + 1
        ny_loc = ymax - ymin + 1

        if not inrange:
            continue # Don't need this file

        filename = os.path.join(path, prefix+"." + str(i) + suffix)
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

    # Force the precision of arrays of dimension>1
    if ndims>1:
        try:
            data = data.astype(t_array.dtype, copy=False)
        except TypeError:
            data = data.astype(t_array.dtype)

    # Finished looping over all files
    if info:
        sys.stdout.write("\n")
    return data
