#!/usr/bin/env python3

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

def squashoutput(datadir=".", outputname="BOUT.dmp.nc", format="NETCDF4", tslice=[None,None,None], xslice=[None,None,None], yslice=[None,None,None], zslice=[None,None,None], singleprecision=False):
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
    tslice : [int, int, int]
        lower, upper, stride values to slice the t-dimension
        default [None, None, None]
    xslice : [int, int, int]
        lower, upper, stride values to slice the x-dimension
        default [None, None, None]
    yslice : [int, int, int]
        lower, upper, stride values to slice the y-dimension
        default [None, None, None]
    zslice : [int, int, int]
        lower, upper, stride values to slice the z-dimension
        default [None, None, None]
    singleprecision : bool
        If true convert data to single-precision floats
        default False
    """

    fullpath = os.path.join(datadir,outputname)
    if os.path.isfile(fullpath):
        raise ValueError(fullpath+" already exists. Collect may try to read from this file, which is presumably not desired behaviour.")

    # Check *slice arguments, and make sure length is at least 3
    if len(tslice) < 3:
        tslice += [None]*(3-len(tslice))
    elif len(tslice) > 3:
        raise ValueError("Can provide at most 3 arguments for tslice")
    if len(xslice) < 3:
        xslice += [None]*(3-len(xslice))
    elif len(xslice) > 3:
        raise ValueError("Can provide at most 3 arguments for xslice")
    if len(yslice) < 3:
        yslice += [None]*(3-len(yslice))
    elif len(yslice) > 3:
        raise ValueError("Can provide at most 3 arguments for yslice")
    if len(zslice) < 3:
        zslice += [None]*(3-len(zslice))
    elif len(zslice) > 3:
        raise ValueError("Can provide at most 3 arguments for zslice")

    # useful object from BOUT pylib to access output data
    outputs = BoutOutputs(datadir, info=False, xguards=True, yguards=True)
    outputvars = outputs.keys()
    t_array_index = outputvars.index("t_array")
    outputvars.append(outputvars.pop(t_array_index))

    # Create single file for output and write data
    with DataFile(fullpath,create=True,write=True,format=format) as f:
        for varname in outputvars:
            print(varname)

            var = outputs[varname]
            if singleprecision:
                var = BoutArray(numpy.float32(var), var.attributes)
            if var.attributes["bout_type"] == "Field3D_t":
                var = var[tslice[0]:tslice[1]:tslice[2], xslice[0]:xslice[1]:xslice[2], yslice[0]:yslice[1]:yslice[2], zslice[0]:zslice[1]:zslice[2]]
            elif var.attributes["bout_type"] == "Field2D_t":
                var = var[tslice[0]:tslice[1]:tslice[2], xslice[0]:xslice[1]:xslice[2], yslice[0]:yslice[1]:yslice[2]]
            elif var.attributes["bout_type"] == "scalar_t":
                var = var[tslice[0]:tslice[1]:tslice[2]]
            elif var.attributes["bout_type"] == "Field3D":
                var = var[xslice[0]:xslice[1]:xslice[2], yslice[0]:yslice[1]:yslice[2], zslice[0]:zslice[1]:zslice[2]]
            elif var.attributes["bout_type"] == "Field2D":
                var = var[xslice[0]:xslice[1]:xslice[2], yslice[0]:yslice[1]:yslice[2]]
            else:
                var = outputs[varname]

            if varname == "NXPE" or varname == "NYPE":
                f.write(varname, 1)
            else:
                f.write(varname, var)

if __name__=="__main__":
    # Call the squashoutput function using arguments from
    # command line when this file is called as an executable

    import argparse
    from sys import exit

    # Parse command line arguments
    parser = argparse.ArgumentParser()
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
    parser.add_argument("--trange", type=int, nargs='*', default=None)
    parser.add_argument("--tslice", type=int_or_none, nargs='*', default=[None])
    parser.add_argument("--xslice", type=int_or_none, nargs='*', default=[None])
    parser.add_argument("--yslice", type=int_or_none, nargs='*', default=[None])
    parser.add_argument("--zslice", type=int_or_none, nargs='*', default=[None])
    parser.add_argument("--singleprecision", type=str_to_bool, default=True)
    args = parser.parse_args()

    if args.trange is not None:
        if len(args.trange) == 1:
            # just pass an int if only one element
            args.trange = args.trange[0]
        elif len(args.trange) > 2:
            # error
            raise ValueError("Can give at most two arguments for trange")
    # Call the function, using command line arguments
    squashoutput(datadir=args.datadir, outputname=args.outputname, trange=args.trange, tslice=args.tslice, xslice=args.xslice, yslice=args.yslice, zslice=args.zslice, singleprecision=args.singleprecision)

    exit(0)
