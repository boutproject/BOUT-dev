"""
Collect all data from BOUT.dmp.* files and create a single output file.

Output file named BOUT.dmp.nc by default

Useful because this discards ghost cell data (that is only useful for debugging)
and because single files are quicker to download.

"""

from boutdata.data import BoutOutputs
from boututils.datafile import DataFile
from boututils.boutarray import BoutArray
import numpy
import os
import gc
import tempfile
import shutil
import glob


def squashoutput(datadir=".", outputname="BOUT.dmp.nc", format="NETCDF4", tind=None,
                 xind=None, yind=None, zind=None, xguards=True, yguards="include_upper",
                 singleprecision=False, compress=False, least_significant_digit=None,
                 quiet=False, complevel=None, append=False, delete=False):
    """
    Collect all data from BOUT.dmp.* files and create a single output file.

    Parameters
    ----------
    datadir : str
        Directory where dump files are and where output file will be created.
        default "."
    outputname : str
        Name of the output file. File suffix specifies whether to use NetCDF or
        HDF5 (see boututils.datafile.DataFile for suffixes).
        default "BOUT.dmp.nc"
    format : str
        format argument passed to DataFile
        default "NETCDF4"
    tind : slice, int, or [int, int, int]
        tind argument passed to collect
        default None
    xind : slice, int, or [int, int, int]
        xind argument passed to collect
        default None
    yind : slice, int, or [int, int, int]
        yind argument passed to collect
        default None
    zind : slice, int, or [int, int, int]
        zind argument passed to collect
        default None
    xguards : bool
        xguards argument passed to collect
        default True
    yguards : bool or "include_upper"
        yguards argument passed to collect (note different default to collect's)
        default "include_upper"
    singleprecision : bool
        If true convert data to single-precision floats
        default False
    compress : bool
        If true enable compression in the output file
    least_significant_digit : int or None
        How many digits should be retained? Enables lossy
        compression. Default is lossless compression. Needs
        compression to be enabled.
    complevel : int or None
        Compression level, 1 should be fastest, and 9 should yield
        highest compression.
    quiet : bool
        Be less verbose. default False
    append : bool
        Append to existing squashed file
    delete : bool
        Delete the original files after squashing.
    """

    fullpath = os.path.join(datadir, outputname)

    if append:
        datadirnew = tempfile.mkdtemp(dir=datadir)
        for f in glob.glob(datadir + "/BOUT.dmp.*.??"):
            if not quiet:
                print("moving", f, flush=True)
            shutil.move(f, datadirnew)
        oldfile = datadirnew + "/" + outputname
        datadir = datadirnew

    if os.path.isfile(fullpath) and not append:
        raise ValueError(
            fullpath + " already exists. Collect may try to read from this file, which is presumably not desired behaviour.")

    # useful object from BOUT pylib to access output data
    outputs = BoutOutputs(datadir, info=False, xguards=xguards,
                          yguards=yguards, tind=tind, xind=xind, yind=yind, zind=zind)
    outputvars = outputs.keys()
    # Read a value to cache the files
    outputs[outputvars[0]]

    if append:
        # move only after the file list is cached
        shutil.move(fullpath, oldfile)

    t_array_index = outputvars.index("t_array")
    outputvars.append(outputvars.pop(t_array_index))

    kwargs = {}
    if compress:
        kwargs['zlib'] = True
        if least_significant_digit is not None:
            kwargs['least_significant_digit'] = least_significant_digit
        if complevel is not None:
            kwargs['complevel'] = complevel
    if append:
        old = DataFile(oldfile)
        # Check if dump on restart was enabled
        # If so, we want to drop the duplicated entry
        cropnew = 0
        if old['t_array'][-1] == outputs['t_array'][0]:
            cropnew = 1
        # Make sure we don't end up with duplicated data:
        for ot in old['t_array']:
            if ot in outputs['t_array'][cropnew:]:
                raise RuntimeError(
                    "For some reason t_array has some duplicated entries in the new and old file.")
    # Create single file for output and write data
    with DataFile(fullpath, create=True, write=True, format=format, **kwargs) as f:
        for varname in outputvars:
            if not quiet:
                print(varname, flush=True)

            var = outputs[varname]
            if append:
                dims = outputs.dimensions[varname]
                if 't' in dims:
                    var = var[cropnew:, ...]
                    varold = old[varname]
                    var = BoutArray(numpy.append(
                        varold, var, axis=0), var.attributes)

            if singleprecision:
                if not isinstance(var, int):
                    var = BoutArray(numpy.float32(var), var.attributes)

            f.write(varname, var)
            # Write changes, free memory
            f.sync()
            var = None
            gc.collect()

    if delete:
        if append:
            os.remove(oldfile)
        for f in glob.glob(datadir + "/BOUT.dmp.*.??"):
            if not quiet:
                print("Deleting", f, flush=True)
            os.remove(f)
        if append:
            os.rmdir(datadir)
