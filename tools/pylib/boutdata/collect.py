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

from boututils.datafile import DataFile
from boututils.boutarray import BoutArray

import os
import sys
import glob

import numpy as np


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
        print("Variable '"+varname +
              "' not found, and is ambiguous. Could be one of: "+str(v))
        raise ValueError("Variable '"+varname+"' not found")

    # None found. Check if it's an abbreviation
    v = [name for name in varlist
         if name[:len(varname)].lower() == varname.lower()]
    if len(v) == 1:
        print("Variable '%s' not found. Using '%s' instead" % (varname, v[0]))
        return v[0]

    if len(v) > 1:
        print("Variable '"+varname +
              "' not found, and is ambiguous. Could be one of: "+str(v))
    raise ValueError("Variable '"+varname+"' not found")


def collect(varname, xind=None, yind=None, zind=None, tind=None, path=".", yguards=False, xguards=True, info=True, prefix="BOUT.dmp", strict=False, tind_auto=False, datafile_cache=None):
    """Collect a variable from a set of BOUT++ outputs.

    data = collect(name)

    varname   Name of the variable (string)

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
    tind_auto = False      Read all files, to get the shortest length of time_indices
                           useful if writing got interrupted.
    datafile_cache = None  Optional cache of open DataFile instances:
                           namedtuple as returned by create_cache. Used by
                           BoutOutputs to pass in a cache so that we do not
                           have to re-open the dump files to read another
                           variable.
    """

    if datafile_cache is None:
        # Search for BOUT++ dump files
        file_list, parallel, suffix = findFiles(path, prefix)
    else:
        parallel = datafile_cache.parallel
        suffix = datafile_cache.suffix
        file_list = datafile_cache.file_list

    # Get the DataFile from the cache, if present, otherwise open the DataFile
    def getDataFile(i):
        if datafile_cache is not None:
            return datafile_cache.datafile_list[i]
        else:
            return DataFile(file_list[i])

    if parallel:
        print("Single (parallel) data file")
        f = getDataFile(0)  # Get the file
        dimens = f.dimensions(varname)
        ndims = f.ndims(varname)
        try:
            mxg = f["MXG"]
        except KeyError:
            mxg = 0
            print("MXG not found, setting to {}".format(mxg))
        try:
            myg = f["MYG"]
        except KeyError:
            myg = 0
            print("MYG not found, setting to {}".format(myg))
            
        if tind is not None or xind is not None or yind is not None or zind is not None:
            raise ValueError("tind, xind, yind, zind arguments are not implemented yet for single (parallel) data file")
        
        if xguards:
            xstart = 0
            xlim = None
        else:
            xstart = mxg
            if mxg > 0:
                xlim = -mxg
            else:
                xlim = None
        if yguards:
            ystart = 0
            ylim = None
        else:
            ystart = mxg
            if myg > 0:
                ylim = -myg
            else:
                ylim = None

        data = f.read(varname)
        attributes = f.attributes(varname)
        if ndims == 2:
            # Field2D
            data = data[xstart:xlim, ystart:ylim]
        elif ndims == 3:
            if dimens[2] == 'z':
                # Field3D
                data = data[xstart:xlim, ystart:ylim, :]
            else:
                # evolving Field2D
                data = data[:, xstart:xlim, ystart:ylim]
        elif ndims == 4:
            # evolving Field3D
            data = data[:, xstart:xlim, ystart:ylim, :]
        return BoutArray(data, attributes=attributes)
    nfiles = len(file_list)

    # Read data from the first file
    f = getDataFile(0)
    attributes = f.attributes(varname)

    try:
        dimens = f.dimensions(varname)
        ndims = f.ndims(varname)
    except:
        if strict:
            raise
        else:
            # Find the variable
            varname = findVar(varname, f.list())

            dimens = f.dimensions(varname)
            ndims = f.ndims(varname)

    # ndims is 0 for reals, and 1 for f.ex. t_array
    if ndims == 0:
        # Just read from file
        data = f.read(varname)
        if datafile_cache is None:
            # close the DataFile if we are not keeping it in a cache
            f.close()
        return BoutArray(data, attributes=attributes)

    if ndims > 4:
        raise ValueError("ERROR: Too many dimensions")

    def load_and_check(varname):
        var = f.read(varname)
        if var is None:
            raise ValueError("Missing " + varname + " variable")
        return var

    mxsub = load_and_check("MXSUB")
    mysub = load_and_check("MYSUB")
    mz = load_and_check("MZ")
    mxg = load_and_check("MXG")
    myg = load_and_check("MYG")
    t_array = f.read("t_array")
    if t_array is None:
        nt = 1
        t_array = np.zeros(1)
    else:
        try:
            nt = len(t_array)
        except TypeError:
            # t_array is not an array here, which probably means it was a
            # one-element array and has been read as a scalar.
            nt = 1
        if tind_auto:
            for i in range(nfiles):
                t_array_ = getDataFile(i).read("t_array")
                nt = min(len(t_array_), nt)

    if info:
        print("mxsub = %d mysub = %d mz = %d\n" % (mxsub, mysub, mz))

    # Get the version of BOUT++ (should be > 0.6 for NetCDF anyway)
    try:
        version = f["BOUT_VERSION"]
    except KeyError:
        print("BOUT++ version : Pre-0.2")
        version = 0
    if version < 3.5:
        # Remove extra point
        nz = mz-1
    else:
        nz = mz

    # Fallback to sensible (?) defaults
    try:
        nxpe = f["NXPE"]
    except KeyError:
        nxpe = 1
        print("NXPE not found, setting to {}".format(nxpe))
    try:
        nype = f["NYPE"]
    except KeyError:
        nype = nfiles
        print("NYPE not found, setting to {}".format(nype))

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

    if yguards:
        ny = mysub * nype + 2*myg
    else:
        ny = mysub * nype

    # Check ranges

    def check_range(r, low, up, name="range"):
        r2 = r
        if r is not None:
            try:
                n = len(r2)
            except:
                # No len attribute, so probably a single number
                r2 = [r2, r2]
            if (len(r2) < 1) or (len(r2) > 2):
                print("WARNING: "+name+" must be [min, max]")
                r2 = None
            else:
                if len(r2) == 1:
                    r2 = [r2, r2]
                if r2[0] < 0 and low >= 0:
                    r2[0] += (up-low+1)
                if r2[1] < 0 and low >= 0:
                    r2[1] += (up-low+1)
                if r2[0] < low:
                    r2[0] = low
                if r2[0] > up:
                    r2[0] = up
                if r2[1] < low:
                    r2[1] = low
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
    zind = check_range(zind, 0, nz-1, "zind")
    tind = check_range(tind, 0, nt-1, "tind")

    xsize = xind[1] - xind[0] + 1
    ysize = yind[1] - yind[0] + 1
    zsize = zind[1] - zind[0] + 1
    tsize = tind[1] - tind[0] + 1

    if ndims == 1:
        if tind is None:
            data = f.read(varname)
        else:
            data = f.read(varname, ranges=[tind[0], tind[1]+1])
        if datafile_cache is None:
            # close the DataFile if we are not keeping it in a cache
            f.close()
        return BoutArray(data, attributes=attributes)

    if datafile_cache is None:
        # close the DataFile if we are not keeping it in a cache
        f.close()

    # Map between dimension names and output size
    sizes = {'x': xsize, 'y': ysize, 'z': zsize, 't': tsize}

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
                if ymax < 0:
                    inrange = False
                if ymin < 0:
                    ymin = 0
            else:
                if ymax < myg:
                    inrange = False
                if ymin < myg:
                    ymin = myg

            # Upper y boundary
            if pe_yind == (nype - 1):
                # Keeping outer boundary
                if ymin >= (mysub + 2*myg):
                    inrange = False
                if ymax > (mysub + 2*myg - 1):
                    ymax = (mysub + 2*myg - 1)
            else:
                if ymin >= (mysub + myg):
                    inrange = False
                if ymax >= (mysub + myg):
                    ymax = (mysub+myg-1)

            # Calculate global indices
            ygmin = ymin + pe_yind * mysub
            ygmax = ymax + pe_yind * mysub

        else:
            # Get local ranges
            ymin = yind[0] - pe_yind*mysub + myg
            ymax = yind[1] - pe_yind*mysub + myg

            if (ymin >= (mysub + myg)) or (ymax < myg):
                inrange = False  # Y out of range

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
                if xmax < 0:
                    inrange = False
                if xmin < 0:
                    xmin = 0
            else:
                if xmax < mxg:
                    inrange = False
                if xmin < mxg:
                    xmin = mxg

            # Upper x boundary
            if pe_xind == (nxpe - 1):
                # Keeping outer boundary
                if xmin >= (mxsub + 2*mxg):
                    inrange = False
                if xmax > (mxsub + 2*mxg - 1):
                    xmax = (mxsub + 2*mxg - 1)
            else:
                if xmin >= (mxsub + mxg):
                    inrange = False
                if xmax >= (mxsub + mxg):
                    xmax = (mxsub+mxg-1)

            # Calculate global indices
            xgmin = xmin + pe_xind * mxsub
            xgmax = xmax + pe_xind * mxsub

        else:
            # Get local ranges
            xmin = xind[0] - pe_xind*mxsub + mxg
            xmax = xind[1] - pe_xind*mxsub + mxg

            if (xmin >= (mxsub + mxg)) or (xmax < mxg):
                inrange = False  # X out of range

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
            continue  # Don't need this file

        if info:
            sys.stdout.write("\rReading from " + file_list[i] + ": [" +
                             str(xmin) + "-" + str(xmax) + "][" +
                             str(ymin) + "-" + str(ymax) + "] -> [" +
                             str(xgmin) + "-" + str(xgmax) + "][" +
                             str(ygmin) + "-" + str(ygmax) + "]")

        f = getDataFile(i)

        if ndims == 4:
            d = f.read(varname, ranges=[tind[0], tind[1]+1,
                                        xmin, xmax+1,
                                        ymin, ymax+1,
                                        zind[0], zind[1]+1])
            data[:, (xgmin-xind[0]):(xgmin-xind[0]+nx_loc),
                 (ygmin-yind[0]):(ygmin-yind[0]+ny_loc), :] = d
        elif ndims == 3:
            # Could be xyz or txy

            if dimens[2] == 'z':  # xyz
                d = f.read(varname, ranges=[xmin, xmax+1,
                                            ymin, ymax+1,
                                            zind[0], zind[1]+1])
                data[(xgmin-xind[0]):(xgmin-xind[0]+nx_loc),
                     (ygmin-yind[0]):(ygmin-yind[0]+ny_loc), :] = d
            else:  # txy
                d = f.read(varname, ranges=[tind[0], tind[1]+1,
                                            xmin, xmax+1,
                                            ymin, ymax+1])
                data[:, (xgmin-xind[0]):(xgmin-xind[0]+nx_loc),
                     (ygmin-yind[0]):(ygmin-yind[0]+ny_loc)] = d
        elif ndims == 2:
            # xy
            d = f.read(varname, ranges=[xmin, xmax+1,
                                        ymin, ymax+1])
            data[(xgmin-xind[0]):(xgmin-xind[0]+nx_loc),
                 (ygmin-yind[0]):(ygmin-yind[0]+ny_loc)] = d

        if datafile_cache is None:
            # close the DataFile if we are not keeping it in a cache
            f.close()

    # Force the precision of arrays of dimension>1
    if ndims > 1:
        try:
            data = data.astype(t_array.dtype, copy=False)
        except TypeError:
            data = data.astype(t_array.dtype)

    # Finished looping over all files
    if info:
        sys.stdout.write("\n")
    return BoutArray(data, attributes=attributes)


def attributes(varname, path=".", prefix="BOUT.dmp"):
    """
    Returns a dictionary of variable attributes

    varname   Name of the variable (string)

    Optional arguments:

    path    = "."          Path to data files
    prefix  = "BOUT.dmp"   File prefix

    """
    # Search for BOUT++ dump files in NetCDF format
    file_list, _, _ = findFiles(path, prefix)

    # Read data from the first file
    f = DataFile(file_list[0])

    return f.attributes(varname)


def dimensions(varname, path=".", prefix="BOUT.dmp"):
    """
    Returns a list of dimensions

    varname   Name of the variable (string)

    Optional arguments:

    path    = "."          Path to data files
    prefix  = "BOUT.dmp"   File prefix
    """
    file_list, _, _ = findFiles(path, prefix)
    return DataFile(file_list[0]).dimensions(varname)


def findFiles(path, prefix):
    """
    Find files matching prefix in path.

    Netcdf (nc) and HDF5 (hdf5) files are searched.

    Returns the list of files, whether the files are a parallel dump file and
    the file suffix.

    """

    # Make sure prefix does not have a trailing .
    if prefix[-1] == ".":
        prefix = prefix[:-1]

    # Look for parallel dump files
    suffixes = [".nc", ".ncdf", ".cdl", ".h5", ".hdf5", ".hdf"]
    file_list_parallel = None
    suffix_parallel = ""
    for test_suffix in suffixes:
        files = glob.glob(os.path.join(path, prefix+test_suffix))
        if files:
            if file_list_parallel:  # Already had a list of files
                raise IOError("Parallel dump files with both {0} and {1} extensions are present. Do not know which to read.".format(
                    suffix, test_suffix))
            suffix_parallel = test_suffix
            file_list_parallel = files

    file_list = None
    suffix = ""
    for test_suffix in suffixes:
        files = glob.glob(os.path.join(path, prefix+".*"+test_suffix))
        if files:
            if file_list:  # Already had a list of files
                raise IOError("Dump files with both {0} and {1} extensions are present. Do not know which to read.".format(
                    suffix, test_suffix))
            suffix = test_suffix
            file_list = files

    if file_list_parallel and file_list:
        raise IOError("Both regular (with suffix {0}) and parallel (with suffix {1}) dump files are present. Do not know which to read.".format(
            suffix_parallel, suffix))
    elif file_list_parallel:
        return file_list_parallel, True, suffix_parallel
    elif file_list:
        # make sure files are in the right order
        nfiles = len(file_list)
        file_list = [os.path.join(path, prefix+"."+str(i)+suffix)
                     for i in range(nfiles)]
        return file_list, False, suffix
    else:
        raise IOError("ERROR: No data files found in path {0}".format(path))


def create_cache(path, prefix):
    """
    Create a list of DataFile objects to be passed repeatedly to collect.

    Stores the cache in a namedtuple along with the file_list, and parallel and suffix attributes
    """

    # define namedtuple to return as the result
    from collections import namedtuple
    datafile_cache_tuple = namedtuple(
        "datafile_cache", ["file_list", "parallel", "suffix", "datafile_list"])

    file_list, parallel, suffix = findFiles(path, prefix)

    cache = []
    for f in file_list:
        cache.append(DataFile(f))

    return datafile_cache_tuple(file_list=file_list, parallel=parallel, suffix=suffix, datafile_list=cache)
