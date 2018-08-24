#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

"""
Collect all data from BOUT.dmp.* files and create a single output file.

Output file named BOUT.dmp.nc by default

Useful because this discards ghost cell data (that is only useful for debugging)
and because single files are quicker to download.

When run as script:
    - first command line argument specifies data directory (default is the
      current directory where the script is run)
    - optional argument "--outputname <name>" can be given to change the name
      of the output file
"""

from boutdata.data import BoutOutputs
from boututils.datafile import DataFile
from boututils.boutarray import BoutArray
import numpy
import os

def squashoutput(datadir=".", outputname="BOUT.dmp.nc", format="NETCDF4", tind=None,
                 xind=None, yind=None, zind=None, singleprecision=False, compress=False,
                 least_significant_digit=None, quiet=False, complevel=None, append=False,
                 delete=False):
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
    """

    fullpath = os.path.join(datadir,outputname)

    if append:
        import tempfile
        import shutil
        import glob
        datadirnew = tempfile.mkdtemp(dir=datadir)
        for f in glob.glob(datadir+"/BOUT.dmp.*.??"):
            if not quiet:
                print("moving",f)
            shutil.move(f,datadirnew)
        oldfile=datadirnew+"/"+outputname
        datadir=datadirnew

    if os.path.isfile(fullpath) and not append:
        raise ValueError(fullpath+" already exists. Collect may try to read from this file, which is presumably not desired behaviour.")

    # useful object from BOUT pylib to access output data
    outputs = BoutOutputs(datadir, info=False, xguards=True, yguards=True, tind=tind, xind=xind, yind=yind, zind=zind)
    outputvars = outputs.keys()
    # Read a value to cache the files
    outputs[outputvars[0]]

    if append:
        # move only after the file list is cached
        shutil.move(fullpath,oldfile)

    t_array_index = outputvars.index("t_array")
    outputvars.append(outputvars.pop(t_array_index))

    kwargs={}
    if compress:
        kwargs['zlib']=True
        if least_significant_digit is not None:
            kwargs['least_significant_digit']=least_significant_digit
        if complevel is not None:
            kwargs['complevel']=complevel
    if append:
        old=DataFile(oldfile)
    # Create single file for output and write data
    with DataFile(fullpath,create=True,write=True,format=format, **kwargs) as f:
        for varname in outputvars:
            if not quiet:
                print(varname)

            var = outputs[varname]
            if append:
                dims=old.dimensions(varname)
                if 't' in dims:
                    varold=old[varname]
                    var=BoutArray(numpy.append(varold,var,axis=0),var.attributes)

            if singleprecision:
                if not isinstance(var, int):
                    var = BoutArray(numpy.float32(var), var.attributes)

            f.write(varname, var)

    if delete:
        if append:
            os.remove(oldfile)
        for f in glob.glob(datadir+"/BOUT.dmp.*.??"):
            if not quiet:
                print("Deleting",f)
            os.remove(f)
        if append:
            os.rmdir(datadir)
if __name__=="__main__":
    # Call the squashoutput function using arguments from
    # command line when this file is called as an executable

    import argparse
    from sys import exit
    try:
        import argcomplete
    except ImportError:
        argcomplete=None

    # Parse command line arguments
    parser = argparse.ArgumentParser(__doc__+"\n\n"+squashoutput.__doc__)
    def str_to_bool(string):
        return string=="True" or string=="true" or string=="T" or string=="t"
    def int_or_none(string):
        try:
            return int(string)
        except ValueError:
            if string=='None' or string=='none':
                return None
            else:
                raise
    parser.add_argument("datadir", nargs='?', default=".")
    parser.add_argument("--outputname",default="BOUT.dmp.nc")
    parser.add_argument("--tind", type=int_or_none, nargs='*', default=[None])
    parser.add_argument("--xind", type=int_or_none, nargs='*', default=[None])
    parser.add_argument("--yind", type=int_or_none, nargs='*', default=[None])
    parser.add_argument("--zind", type=int_or_none, nargs='*', default=[None])
    parser.add_argument("-s","--singleprecision", action="store_true", default=False)
    parser.add_argument("-c","--compress", action="store_true", default=False)
    parser.add_argument("-l","--complevel", type=int_or_none, default=None)
    parser.add_argument("-i","--least-significant-digit", type=int_or_none, default=None)
    parser.add_argument("-q","--quiet", action="store_true", default=False)
    parser.add_argument("-a","--append", action="store_true", default=False)
    parser.add_argument("-d","--delete", action="store_true", default=False)

    if argcomplete:
        argcomplete.autocomplete(parser)

    args = parser.parse_args()

    for ind in "txyz":
        args.__dict__[ind+"ind"]=slice(*args.__dict__[ind+"ind"])
    # Call the function, using command line arguments
    squashoutput(**args.__dict__)

    exit(0)
