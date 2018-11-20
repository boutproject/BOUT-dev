"""
Collect all data from BOUT.dmp.* files and create a single output file.

Output file named BOUT.dmp.nc by default

Useful because this discards ghost cell data (that is only useful for debugging)
and because single files are quicker to download.
"""

# imports in function for fast bash completion


def squashoutput(datadir=".", outputname="BOUT.dmp.nc", format="NETCDF4", tind=None,
                 xind=None, yind=None, zind=None, singleprecision=False, compress=False,
                 least_significant_digit=None, quiet=False, complevel=None, append=False,
                 delete=False, progress=False, docontinue=False, earlyExit=None):
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
    progress : bool
        Print a progress bar
    docontinue : bool
        Try to progress a previously interrupted squash
    earlyExit : string
        Exit after variable earlyExit is processed. Mostly for testing.
    """

    from boutdata.data import BoutOutputs
    from boututils.datafile import DataFile
    from boututils.boutarray import BoutArray
    from zoidberg import progress as bar
    import numpy
    import os
    import gc
    import tempfile
    import shutil
    import glob

    if "/" in outputname:
        raise ValueError("only simple filenames are supported for outputname")

    fullpath = os.path.join(datadir, outputname)

    if os.path.isfile(fullpath) and not (append or docontinue):
        raise ValueError(
            fullpath + " already exists. Collect may try to read from this file, which is presumably not desired behaviour.")

    if docontinue or append:
        # move temporary to subdir, so that when the BoutOutputs
        # caches the file list, this file is not found
        datadirtmp = tempfile.mkdtemp(dir=datadir)
        shutil.move(fullpath, datadirtmp)

    # useful object from BOUT pylib to access output data
    outputs = BoutOutputs(datadir, info=False, xguards=True,
                          yguards=True, tind=tind, xind=xind, yind=yind, zind=zind)
    outputvars = outputs.keys()
    # Read a value to cache the files
    outputs[outputvars[0]]

    #?if append:
    #?    # move only after the file list is cached
    #?    shutil.move(fullpath, oldfile)

    if docontinue or append:
        shutil.move(os.path.join(datadirtmp, outputname.split("/")[-1]), datadir)
        os.rmdir(datadirtmp)

    t_array_index = outputvars.index("t_array")
    outputvars.insert(0,outputvars.pop(t_array_index))

    if progress:
        # outputs.sizes() returns for each variable a list with the
        # length in each dimension. We are interrested in the total
        # length of each variable.
        sizes_ = outputs.sizes()
        sizes = {}
        total = 0
        for var, size in sizes_.items():
            cur = 1
            for s in size:
                cur *= s
            total += cur
            sizes[var] = cur
        sizes_ = None
        done = 0

    kwargs = {}
    if compress:
        kwargs['zlib'] = True
        if least_significant_digit is not None:
            kwargs['least_significant_digit'] = least_significant_digit
        if complevel is not None:
            kwargs['complevel'] = complevel
    # Determine where we want to write the data
    toffset = 0
    if append or docontinue:
        # We might be in continue mode
        told = DataFile(fullpath)['t_array'][:]
        # Check if dump on restart was enabled
        # If so, we want to drop the duplicated entry
        toffset = len(told)
        if outputs['t_array'][0] in told:
            toffset = list(told[:]).index(outputs['t_array'][0])
        # Try to cleanup DataFile - is not threadsafe
        gc.collect()
    tmax = toffset + len(outputs['t_array'])

    # Only if nether continue or append mode
    create = not (docontinue or append)

    with DataFile(fullpath, create=create, write=True, **kwargs) as f:
        for varname in outputvars:
            if not quiet:
                print(varname)

            dims = outputs.dimensions[varname]

            if varname in f.keys():
                # Ether not evolved, or already fully read
                if 't' not in dims or f.getRealTlength(varname) == tmax:# and varname != 'tt':
                    if progress:
                        done += sizes[varname]
                        bar.update_progress(done / total, zoidberg=True)
                    continue

            var = outputs[varname]
            if 't' in dims:
                ranges=[slice(toffset,tmax),...]
            else:
                ranges=[...] if len(dims) else None

            if singleprecision:
                if not isinstance(var, int):
                    var = BoutArray(numpy.float32(var), var.attributes)

            f.write(varname, var, ranges=ranges)

            if progress:
                done += sizes[varname]
                bar.update_progress(done / total, zoidberg=True)

            # Write changes, free memory
            f.sync()
            var = None
            gc.collect()

            if earlyExit is not None:
                if earlyExit == varname:
                    return
                    #sys.exit(1)
                    #raise RuntimeError("Trigger error!")

    if delete:
        for f in glob.glob(datadir + "/BOUT.dmp.*.*"):
            if not quiet:
                print("Deleting", f)
            os.remove(f)
