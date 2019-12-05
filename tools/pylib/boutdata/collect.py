from __future__ import print_function
from __future__ import division

from builtins import str, range

import os
import sys
import glob

import numpy as np

from boututils.datafile import DataFile
from boututils.boutarray import BoutArray


def findVar(varname, varlist):
    """Find variable name in a list

    First does case insensitive comparison, then
    checks for abbreviations.

    Returns the matched string, or raises a ValueError

    Parameters
    ----------
    varname : str
        Variable name to look for
    varlist : list of str
        List of possible variable names

    Returns
    -------
    str
        The closest match to varname in varlist

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


def _convert_to_nice_slice(r, N, name="range"):
    """Convert r to a "sensible" slice in range [0, N]

    If r is None, the slice corresponds to the full range.

    Lists or tuples of one or two ints are converted to slices.

    Slices with None for one or more arguments have them replaced with
    sensible values.

    Private helper function for collect

    Parameters
    ----------
    r : None, int, slice or list of int
        Range-like to check/convert to slice
    N : int
        Size of range
    name : str, optional
        Name of range for error message

    Returns
    -------
    slice
        "Sensible" slice with no Nones for start, stop or step
    """

    if N == 0:
        raise ValueError("No data available in %s"%name)
    if r is None:
        temp_slice = slice(N)
    elif isinstance(r, slice):
        temp_slice = r
    elif isinstance(r, (int, np.integer)):
        if r >= N or r <-N:
            # raise out of bounds error as if we'd tried to index the array with r
            # without this, would return an empty array instead
            raise IndexError(name+" index out of range, value was "+str(r))
        elif r == -1:
            temp_slice = slice(r, None)
        else:
            temp_slice = slice(r, r + 1)
    elif len(r) == 0:
        return _convert_to_nice_slice(None, N, name)
    elif len(r) == 1:
        return _convert_to_nice_slice(r[0], N, name)
    elif len(r) == 2:
        r2 = list(r)
        if r2[0] < 0:
            r2[0] += N
        if r2[1] < 0:
            r2[1] += N
        if r2[0] > r2[1]:
            raise ValueError("{} start ({}) is larger than end ({})"
                             .format(name, *r2))
        # Lists uses inclusive end, we need exclusive end
        temp_slice = slice(r2[0], r2[1] + 1)
    else:
        raise ValueError("Couldn't convert {} ('{}') to slice. Please pass a "
                         "slice(start, stop, step) if you need to set a step."
                         .format(name, r))

    # slice.indices converts None to actual values
    return slice(*temp_slice.indices(N))


def collect(varname, xind=None, yind=None, zind=None, tind=None, path=".",
            yguards=False, xguards=True, info=True, prefix="BOUT.dmp",
            strict=False, tind_auto=False, datafile_cache=None):
    """Collect a variable from a set of BOUT++ outputs.

    Parameters
    ----------
    varname : str
        Name of the variable
    xind, yind, zind, tind : int, slice or list of int, optional
        Range of X, Y, Z or time indices to collect. Either a single
        index to collect, a list containing [start, end] (inclusive
        end), or a slice object (usual python indexing). Default is to
        fetch all indices
    path : str, optional
        Path to data files (default: ".")
    prefix : str, optional
        File prefix (default: "BOUT.dmp")
    yguards : bool, optional
        Collect Y boundary guard cells? (default: False)
    xguards : bool, optional
        Collect X boundary guard cells? (default: True)
        (Set to True to be consistent with the definition of nx)
    info : bool, optional
        Print information about collect? (default: True)
    strict : bool, optional
        Fail if the exact variable name is not found? (default: False)
    tind_auto : bool, optional
        Read all files, to get the shortest length of time_indices.
        Useful if writing got interrupted (default: False)
    datafile_cache : datafile_cache_tuple, optional
        Optional cache of open DataFile instances: namedtuple as returned
        by create_cache. Used by BoutOutputs to pass in a cache so that we
        do not have to re-open the dump files to read another variable
        (default: None)

    Examples
    --------

    >>> collect(name)
    BoutArray([[[[...]]]])

    """

    if datafile_cache is None:
        # Search for BOUT++ dump files
        file_list, parallel, _ = findFiles(path, prefix)
    else:
        parallel = datafile_cache.parallel
        file_list = datafile_cache.file_list

    def getDataFile(i):
        """Get the DataFile from the cache, if present, otherwise open the
        DataFile

        """
        if datafile_cache is not None:
            return datafile_cache.datafile_list[i]
        else:
            return DataFile(file_list[i])

    if parallel:
        if info:
            print("Single (parallel) data file")

        f = getDataFile(0)

        dimensions = f.dimensions(varname)

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

        if xguards:
            nx = f["nx"]
        else:
            nx = f["nx"] - 2*mxg
        if yguards:
            ny = f["ny"] + 2*myg
        else:
            ny = f["ny"]
        nz = f["MZ"]
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

        xind = _convert_to_nice_slice(xind, nx, "xind")
        yind = _convert_to_nice_slice(yind, ny, "yind")
        zind = _convert_to_nice_slice(zind, nz, "zind")
        tind = _convert_to_nice_slice(tind, nt, "tind")

        if not xguards:
            xind = slice(xind.start+mxg, xind.stop+mxg, xind.step)
        if not yguards:
            yind = slice(yind.start+myg, yind.stop+myg, yind.step)

        if len(dimensions) == ():
            ranges = []
        elif dimensions == ('t'):
            ranges = [tind]
        elif dimensions == ('x', 'y'):
            # Field2D
            ranges = [xind, yind]
        elif dimensions == ('x', 'z'):
            # FieldPerp
            ranges = [xind, zind]
        elif dimensions == ('t', 'x', 'y'):
            # evolving Field2D
            ranges = [tind, xind, yind]
        elif dimensions == ('t', 'x', 'z'):
            # evolving FieldPerp
            ranges = [tind, xind, zind]
        elif dimensions == ('x', 'y', 'z'):
            # Field3D
            ranges = [xind, yind, zind]
        elif dimensions == ('t', 'x', 'y', 'z'):
            # evolving Field3D
            ranges = [tind, xind, yind, zind]
        else:
            raise ValueError("Variable has incorrect dimensions ({})"
                             .format(dimensions))

        data = f.read(varname, ranges)
        var_attributes = f.attributes(varname)
        return BoutArray(data, attributes=var_attributes)

    nfiles = len(file_list)

    # Read data from the first file
    f = getDataFile(0)

    dimensions = f.dimensions(varname)

    if varname not in f.keys():
        if strict:
            raise ValueError("Variable '{}' not found".format(varname))
        else:
            varname = findVar(varname, f.list())

    var_attributes = f.attributes(varname)
    ndims = len(dimensions)

    # ndims is 0 for reals, and 1 for f.ex. t_array
    if ndims == 0:
        # Just read from file
        data = f.read(varname)
        if datafile_cache is None:
            # close the DataFile if we are not keeping it in a cache
            f.close()
        return BoutArray(data, attributes=var_attributes)

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

    xind = _convert_to_nice_slice(xind, nx, "xind")
    yind = _convert_to_nice_slice(yind, ny, "yind")
    zind = _convert_to_nice_slice(zind, nz, "zind")
    tind = _convert_to_nice_slice(tind, nt, "tind")

    xsize = xind.stop - xind.start
    ysize = yind.stop - yind.start
    zsize = int(np.ceil(float(zind.stop - zind.start)/zind.step))
    tsize = int(np.ceil(float(tind.stop - tind.start)/tind.step))

    if ndims == 1:
        if tind is None:
            data = f.read(varname)
        else:
            data = f.read(varname, ranges=[tind])
        if datafile_cache is None:
            # close the DataFile if we are not keeping it in a cache
            f.close()
        return BoutArray(data, attributes=var_attributes)

    if datafile_cache is None:
        # close the DataFile if we are not keeping it in a cache
        f.close()

    # Map between dimension names and output size
    sizes = {'x': xsize, 'y': ysize, 'z': zsize, 't': tsize}

    # Create a list with size of each dimension
    ddims = [sizes[d] for d in dimensions]

    # Create the data array
    data = np.zeros(ddims)

    if dimensions == ('t', 'x', 'z') or dimensions == ('x', 'z'):
        yindex_global = None
        # The pe_yind that this FieldPerp is going to be read from
        fieldperp_yproc = None

    for i in range(npe):
        # Get X and Y processor indices
        pe_yind = int(i/nxpe)
        pe_xind = i % nxpe

        inrange = True

        if yguards:
            # Get local ranges
            ystart = yind.start - pe_yind*mysub
            ystop = yind.stop - pe_yind*mysub

            # Check lower y boundary
            if pe_yind == 0:
                # Keeping inner boundary
                if ystop <= 0:
                    inrange = False
                if ystart < 0:
                    ystart = 0
            else:
                if ystop < myg-1:
                    inrange = False
                if ystart < myg:
                    ystart = myg

            # Upper y boundary
            if pe_yind == (nype - 1):
                # Keeping outer boundary
                if ystart >= (mysub + 2*myg):
                    inrange = False
                if ystop > (mysub + 2*myg):
                    ystop = (mysub + 2*myg)
            else:
                if ystart >= (mysub + myg):
                    inrange = False
                if ystop > (mysub + myg):
                    ystop = (mysub + myg)

            # Calculate global indices
            ygstart = ystart + pe_yind * mysub
            ygstop = ystop + pe_yind * mysub

        else:
            # Get local ranges
            ystart = yind.start - pe_yind*mysub + myg
            ystop = yind.stop - pe_yind*mysub + myg

            if (ystart >= (mysub + myg)) or (ystop <= myg):
                inrange = False  # Y out of range

            if ystart < myg:
                ystart = myg
            if ystop > mysub + myg:
                ystop = myg + mysub

            # Calculate global indices
            ygstart = ystart + pe_yind * mysub - myg
            ygstop = ystop + pe_yind * mysub - myg

        if xguards:
            # Get local ranges
            xstart = xind.start - pe_xind*mxsub
            xstop = xind.stop - pe_xind*mxsub

            # Check lower x boundary
            if pe_xind == 0:
                # Keeping inner boundary
                if xstop <= 0:
                    inrange = False
                if xstart < 0:
                    xstart = 0
            else:
                if xstop <= mxg:
                    inrange = False
                if xstart < mxg:
                    xstart = mxg

            # Upper x boundary
            if pe_xind == (nxpe - 1):
                # Keeping outer boundary
                if xstart >= (mxsub + 2*mxg):
                    inrange = False
                if xstop > (mxsub + 2*mxg):
                    xstop = (mxsub + 2*mxg)
            else:
                if xstart >= (mxsub + mxg):
                    inrange = False
                if xstop > (mxsub + mxg):
                    xstop = (mxsub+mxg)

            # Calculate global indices
            xgstart = xstart + pe_xind * mxsub
            xgstop = xstop + pe_xind * mxsub

        else:
            # Get local ranges
            xstart = xind.start - pe_xind*mxsub + mxg
            xstop = xind.stop - pe_xind*mxsub + mxg

            if (xstart >= (mxsub + mxg)) or (xstop <= mxg):
                inrange = False  # X out of range

            if xstart < mxg:
                xstart = mxg
            if xstop > mxsub + mxg:
                xstop = mxg + mxsub

            # Calculate global indices
            xgstart = xstart + pe_xind * mxsub - mxg
            xgstop = xstop + pe_xind * mxsub - mxg

        # Number of local values
        nx_loc = xstop - xstart
        ny_loc = ystop - ystart

        if not inrange:
            continue  # Don't need this file

        if info:
            sys.stdout.write("\rReading from " + file_list[i] + ": [" +
                             str(xstart) + "-" + str(xstop-1) + "][" +
                             str(ystart) + "-" + str(ystop-1) + "] -> [" +
                             str(xgstart) + "-" + str(xgstop-1) + "][" +
                             str(ygstart) + "-" + str(ygstop-1) + "]")

        f = getDataFile(i)

        if dimensions == ('t', 'x', 'y', 'z'):
            d = f.read(varname, ranges=[tind,
                                        slice(xstart, xstop),
                                        slice(ystart, ystop),
                                        zind])
            data[:, (xgstart-xind.start):(xgstart-xind.start+nx_loc),
                 (ygstart-yind.start):(ygstart-yind.start+ny_loc), :] = d
        elif dimensions == ('x', 'y', 'z'):
            d = f.read(varname, ranges=[slice(xstart, xstop),
                                        slice(ystart, ystop),
                                        zind])
            data[(xgstart-xind.start):(xgstart-xind.start+nx_loc),
                 (ygstart-yind.start):(ygstart-yind.start+ny_loc), :] = d
        elif dimensions == ('t', 'x', 'y'):
            d = f.read(varname, ranges=[tind,
                                        slice(xstart, xstop),
                                        slice(ystart, ystop)])
            data[:, (xgstart-xind.start):(xgstart-xind.start+nx_loc),
                 (ygstart-yind.start):(ygstart-yind.start+ny_loc)] = d
        elif dimensions == ('t', 'x', 'z'):
            # FieldPerp should only be defined on processors which contain its yindex_global
            f_attributes = f.attributes(varname)
            temp_yindex = f_attributes["yindex_global"]

            if temp_yindex >= 0:
                if yindex_global is None:
                    yindex_global = temp_yindex

                    # we have found a file with containing the FieldPerp, get the attributes from here
                    var_attributes = f_attributes
                assert temp_yindex == yindex_global

            if temp_yindex >= 0:
                # Check we only read from one pe_yind
                assert fieldperp_yproc is None or fieldperp_yproc == pe_yind

                fieldperp_yproc = pe_yind

                d = f.read(varname, ranges=[tind,
                                            slice(xstart, xstop),
                                            zind])
                data[:, (xgstart-xind.start):(xgstart-xind.start+nx_loc), :] = d
        elif dimensions == ('x', 'y'):
            d = f.read(varname, ranges=[slice(xstart, xstop),
                                        slice(ystart, ystop)])
            data[(xgstart-xind.start):(xgstart-xind.start+nx_loc),
                 (ygstart-yind.start):(ygstart-yind.start+ny_loc)] = d
        elif dimensions == ('x', 'z'):
            # FieldPerp should only be defined on processors which contain its yindex_global
            f_attributes = f.attributes(varname)
            temp_yindex = f_attributes["yindex_global"]

            if temp_yindex >= 0:
                if yindex_global is None:
                    yindex_global = temp_yindex

                    # we have found a file with containing the FieldPerp, get the attributes from here
                    var_attributes = f_attributes
                assert temp_yindex == yindex_global

            if temp_yindex >= 0:
                # Check we only read from one pe_yind
                assert fieldperp_yproc is None or fieldperp_yproc == pe_yind

                fieldperp_yproc = pe_yind

                d = f.read(varname, ranges=[slice(xstart, xstop), zind])
                data[(xgstart-xind.start):(xgstart-xind.start+nx_loc), :] = d
        else:
            raise ValueError('Incorrect dimensions '+str(dimensions)+' in collect')

        if datafile_cache is None:
            # close the DataFile if we are not keeping it in a cache
            f.close()

    # if a step was requested in x or y, need to apply it here
    if xind.step is not None or yind.step is not None:
        if dimensions == ('t', 'x', 'y', 'z'):
            data = data[:, ::xind.step, ::yind.step]
        elif dimensions == ('x', 'y', 'z'):
            data = data[::xind.step, ::yind.step, :]
        elif dimensions == ('t', 'x', 'y'):
            data = data[:, ::xind.step, ::yind.step]
        elif dimensions == ('t', 'x', 'z'):
            data = data[:, ::xind.step, :]
        elif dimensions == ('x', 'y'):
            data = data[::xind.step, ::yind.step]
        elif dimensions == ('x', 'z'):
            data = data[::xind.step, :]
        else:
            raise ValueError('Incorrect dimensions '+str(dimensions)+' applying steps in collect')

    # Force the precision of arrays of dimension>1
    if ndims > 1:
        try:
            data = data.astype(t_array.dtype, copy=False)
        except TypeError:
            data = data.astype(t_array.dtype)

    # Finished looping over all files
    if info:
        sys.stdout.write("\n")
    return BoutArray(data, attributes=var_attributes)


def attributes(varname, path=".", prefix="BOUT.dmp"):
    """Return a dictionary of variable attributes in an output file

    Parameters
    ----------
    varname : str
        Name of the variable
    path : str, optional
        Path to data files (default: ".")
    prefix : str, optional
        File prefix (default: "BOUT.dmp")

    Returns
    -------
    dict
        A dictionary of attributes of varname
    """
    # Search for BOUT++ dump files in NetCDF format
    file_list, _, _ = findFiles(path, prefix)

    # Read data from the first file
    f = DataFile(file_list[0])

    return f.attributes(varname)


def dimensions(varname, path=".", prefix="BOUT.dmp"):
    """Return the names of dimensions of a variable in an output file

    Parameters
    ----------
    varname : str
        Name of the variable
    path : str, optional
        Path to data files (default: ".")
    prefix : str, optional
        File prefix (default: "BOUT.dmp")

    Returns
    -------
    tuple of strs
        The elements of the tuple give the names of corresponding variable
        dimensions

    """
    file_list, _, _ = findFiles(path, prefix)
    return DataFile(file_list[0]).dimensions(varname)


def findFiles(path, prefix):
    """Find files matching prefix in path.

    Netcdf (".nc", ".ncdf", ".cdl") and HDF5 (".h5", ".hdf5", ".hdf")
    files are searched.

    Parameters
    ----------
    path : str
        Path to data files
    prefix : str
        File prefix

    Returns
    -------
    tuple : (list of str, bool, str)
        The first element of the tuple is the list of files, the second is
        whether the files are a parallel dump file and the last element is
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
    """Create a list of DataFile objects to be passed repeatedly to
    collect.

    Parameters
    ----------
    path : str
        Path to data files
    prefix : str
        File prefix

    Returns
    -------
    namedtuple : (list of str, bool, str, list of :py:obj:`~boututils.datafile.DataFile`)
        The cache of DataFiles in a namedtuple along with the file_list,
        and parallel and suffix attributes

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
