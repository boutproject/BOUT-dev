"""Routines for manipulating restart files

TODO
----

- Don't import ``numpy.random.normal`` directly, just the ``random``
  submodule, or sphinx includes the documentation for ``normal``

"""

from __future__ import print_function
from __future__ import division
from builtins import str, range

import os
import glob

from boutdata.collect import collect, create_cache
from boututils.datafile import DataFile
from boututils.boutarray import BoutArray
from boutdata.processor_rearrange import get_processor_layout, create_processor_layout

import multiprocessing
import numpy as np
from numpy import mean, zeros, arange
from numpy.random import normal

from scipy.interpolate import interp1d
try:
    from scipy.interpolate import RegularGridInterpolator
except ImportError:
    pass

def resize3DField(var, data, coordsAndSizesTuple, method, mute):
    """Resize 3D fields

    To be called by resize.

    Written as a function in order to call it using multiprocess. Must
    be defined as a top level function in order to be pickable by the
    multiprocess.

    See the function resize for details

    """

    # Unpack the tuple for better readability
    xCoordOld, yCoordOld, zCoordOld,\
        xCoordNew, yCoordNew, zCoordNew,\
        newNx, newNy, newNz = coordsAndSizesTuple

    if not(mute):
        print("    Resizing "+var +
              ' to (nx,ny,nz) = ({},{},{})'.format(newNx, newNy, newNz))

    # Make the regular grid function (see examples in
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html
    # for details)
    gridInterpolator = RegularGridInterpolator(
        (xCoordOld, yCoordOld, zCoordOld), data, method)

    # Need to fill with one exrta z plane (will only contain zeros)
    newData = np.zeros((newNx, newNy, newNz))

    # Interpolate to the new values
    for xInd, x in enumerate(xCoordNew):
        for yInd, y in enumerate(yCoordNew):
            for zInd, z in enumerate(zCoordNew):
                newData[xInd, yInd, zInd] = gridInterpolator([x, y, z])

    return var, newData


def resize(newNx, newNy, newNz, mxg=2, myg=2,
           path="data", output="./", informat="nc", outformat=None,
           method='linear', maxProc=None, mute=False):
    """Increase/decrease the number of points in restart files.

    NOTE: Can't overwrite
    WARNING: Currently only implemented with uniform BOUT++ grid

    Parameters
    ----------
    newNx, newNy, newNz : int
        nx, ny, nz for the new file (including ghost points)
    mxg, myg : int, optional
        Number of ghost points in x, y (default: 2)
    path : str, optional
        Input path to data files
    output : str, optional
        Path to write new files
    informat : str, optional
        File extension of input
    outformat : {None, str}, optional
        File extension of output (default: use the same as `informat`)
    method : {'linear', 'nearest'}, optional
        What interpolation method to be used
    maxProc : {None, int}, optional
        Limits maximum processors to use when interpolating if set
    mute : bool, optional
        Whether or not output should be printed from this function

    Returns
    -------
    return : bool
        True on success, else False

    TODO
    ----
    - Add 2D field interpolation
    - Replace printing errors with raising `ValueError`
    - Make informat work like `redistribute`

    """

    if method is None:
        # Make sure the method is set
        method = 'linear'

    if outformat is None:
        outformat = informat

    if path == output:
        print("ERROR: Can't overwrite restart files when expanding")
        return False

    def is_pow2(x):
        """Returns true if x is a power of 2"""
        return (x > 0) and ((x & (x-1)) == 0)

    if not is_pow2(newNz):
        print("ERROR: New Z size {} must be a power of 2".format(newNz))
        return False

    file_list = glob.glob(os.path.join(path, "BOUT.restart.*."+informat))
    file_list.sort()
    nfiles = len(file_list)

    if nfiles == 0:
        print("ERROR: No data found in {}".format(path))
        return False

    if not(mute):
        print("Number of files found: " + str(nfiles))

    for f in file_list:
        new_f = os.path.join(output, f.split('/')[-1])
        if not(mute):
            print("Changing {} => {}".format(f, new_f))

        # Open the restart file in read mode and create the new file
        with DataFile(f) as old, DataFile(new_f, write=True, create=True) as new:

            # Find the dimension
            for var in old.list():
                # Read the data
                data = old.read(var)
                # Find 3D variables
                if old.ndims(var) == 3:
                    break

            nx, ny, nz = data.shape
            # Make coordinates
            # NOTE: The max min of the coordinates are irrelevant when
            #       interpolating (as long as old and new coordinates
            #       are consistent), so we just choose all variable to
            #       be between 0 and 1 Calculate the old coordinates
            xCoordOld = np.linspace(0, 1, nx)
            yCoordOld = np.linspace(0, 1, ny)
            zCoordOld = np.linspace(0, 1, nz)

            # Calculate the new coordinates
            xCoordNew = np.linspace(xCoordOld[0], xCoordOld[-1], newNx)
            yCoordNew = np.linspace(yCoordOld[0], yCoordOld[-1], newNy)
            zCoordNew = np.linspace(zCoordOld[0], zCoordOld[-1], newNz)

            # Make a pool of workers
            pool = multiprocessing.Pool(maxProc)
            # List of jobs and results
            jobs = []
            # Pack input to resize3DField together
            coordsAndSizesTuple = (xCoordOld, yCoordOld, zCoordOld,
                                   xCoordNew, yCoordNew, zCoordNew,
                                   newNx, newNy, newNz)

            # Loop over the variables in the old file
            for var in old.list():
                # Read the data
                data = old.read(var)
                attributes = old.attributes(var)

                # Find 3D variables
                if old.ndims(var) == 3:

                    # Asynchronous call (locks first at .get())
                    jobs.append(pool.apply_async(resize3DField,
                                                 args=(var, data, coordsAndSizesTuple, method, mute, )))

                else:
                    if not(mute):
                        print("    Copying "+var)
                        newData = data.copy()
                    if not(mute):
                        print("Writing "+var)
                    new.write(var, newData)

            for job in jobs:
                var, newData = job.get()
                newData = BoutArray(newData, attributes=attributes)
                if not(mute):
                    print("Writing "+var)
                new.write(var, newData)

            # Close the pool of workers
            pool.close()
            # Wait for all processes to finish
            pool.join()

    return True


def resizeZ(newNz, path="data", output="./", informat="nc", outformat=None):
    """Increase the number of Z points in restart files

    NOTE:
        * Can't overwrite
        * Will not yield a result close to the original if there are
          asymmetries in the z-direction

    Parameters
    ----------
    newNz : int
        nz for the new file
    path : str, optional
        Path to original restart files (default: "data")
    output : str, optional
        Path to write new restart files (default: current directory)
    informat : str, optional
        File extension of original files (default: "nc")
    outformat : str, optional
        File extension of new files (default: use the same as `informat`)

    Returns
    -------
    True on success, else False

    TODO
    ----
    - Replace printing errors with raising `ValueError`
    - Make informat work like `redistribute`

    """

    if outformat is None:
        outformat = informat

    if path == output:
        print("ERROR: Can't overwrite restart files when expanding")
        return False

    def is_pow2(x):
        """Returns true if x is a power of 2"""
        return (x > 0) and ((x & (x-1)) == 0)

    if not is_pow2(newNz):
        print("ERROR: New Z size must be a power of 2")
        return False

    file_list = glob.glob(os.path.join(path, "BOUT.restart.*."+informat))
    file_list.sort()
    nfiles = len(file_list)

    if nfiles == 0:
        print("ERROR: No data found")
        return False

    print("Number of files found: " + str(nfiles))

    for f in file_list:
        new_f = os.path.join(output, f.split('/')[-1])
        print("Changing {} => {}".format(f, new_f))

        # Open the restart file in read mode and create the new file
        with DataFile(f) as old,\
                DataFile(new_f, write=True, create=True) as new:
            # Loop over the variables in the old file
            for var in old.list():
                # Read the data
                data = old.read(var)
                attributes = old.attributes(var)

                # Find 3D variables
                if old.ndims(var) == 3:
                    print("    Resizing "+var)

                    nx, ny, nz = data.shape

                    newdata = np.zeros((nx, ny, newNz))
                    for x in range(nx):
                        for y in range(ny):
                            f_old = np.fft.fft(data[x, y, :])

                            # Number of points in f is power of 2
                            f_new = np.zeros(newNz)

                            # Copy coefficients across (ignoring Nyquist)
                            f_new[0] = f_old[0]  # DC
                            for m in range(1, int(nz/2)):
                                # + ve frequencies
                                f_new[m] = f_old[m]
                                # - ve frequencies
                                f_new[newNz-m] = f_old[nz-m]

                            # Invert fft
                            newdata[x, y, :] = np.fft.ifft(f_new).real
                            newdata[x, y, :] = newdata[x, y, 0]

                    # Multiply with the ratio of newNz/nz
                    # This is not needed in the IDL routine as the
                    # forward transfrom has the scaling factor 1/N in
                    # the forward transform, whereas the scaling factor
                    # 1/N is the inverse transform in np.fft
                    # Note that ifft(fft(a)) = a for the same number of
                    # points in both IDL and np.ftt
                    newdata *= (newNz/nz)
                else:
                    print("    Copying "+var)
                    newdata = data.copy()

                newdata = BoutArray(newdata, attributes=attributes)

                new.write(var, newdata)

    return True


def addnoise(path=".", var=None, scale=1e-5):
    """Add random noise to restart files

    .. warning:: Modifies restart files in place! This is in contrast
                 to most of the functions in this module!

    Parameters
    ----------
    path : str, optional
        Path to restart files (default: current directory)
    var : str, optional
        The variable to modify. By default all 3D variables are modified
    scale : float
        Amplitude of the noise. Gaussian noise is used, with zero mean
        and this parameter as the standard deviation

    """
    file_list = glob.glob(os.path.join(path, "BOUT.restart.*"))
    nfiles = len(file_list)

    print("Number of restart files: %d" % (nfiles,))

    for file in file_list:
        print(file)
        with DataFile(file, write=True) as d:
            if var is None:
                for v in d.list():
                    if d.ndims(v) == 3:
                        print(" -> "+v)
                        data = d.read(v, asBoutArray=True)
                        data += normal(scale=scale, size=data.shape)
                        d.write(v, data)
            else:
                # Modify a single variable
                print(" -> "+var)
                data = d.read(var)
                data += normal(scale=scale, size=data.shape)
                d.write(var, data)


def scalevar(var, factor, path="."):
    """Scales a variable by a given factor, modifying restart files in
    place

    .. warning:: Modifies restart files in place! This is in contrast
                 to most of the functions in this module!

    Parameters
    ----------
    var : str
        Name of the variable
    factor : float
        Factor to multiply
    path : str, optional
        Path to the restart files (default: current directory)

    """

    file_list = glob.glob(os.path.join(path, "BOUT.restart.*"))
    nfiles = len(file_list)

    print("Number of restart files: %d" % (nfiles,))
    for file in file_list:
        print(file)
        with DataFile(file, write=True) as d:
            d[var] = d[var] * factor


def create(averagelast=1, final=-1, path="data", output="./", informat="nc", outformat=None):
    """Create restart files from data (dmp) files.

    Parameters
    ----------
    averagelast : int, optional
        Number of time points (counting from `final`, inclusive) to
        average over (default is 1 i.e. just take last time-point)
    final : int, optional
        The last time point to use (default is last, -1)
    path : str, optional
        Path to original restart files (default: "data")
    output : str, optional
        Path to write new restart files (default: current directory)
    informat : str, optional
        File extension of original files (default: "nc")
    outformat : str, optional
        File extension of new files (default: use the same as `informat`)

    """

    if outformat is None:
        outformat = informat

    file_list = glob.glob(os.path.join(path, "BOUT.dmp.*."+informat))
    nfiles = len(file_list)

    print(("Number of data files: ", nfiles))

    for i in range(nfiles):
        # Open each data file
        infname = os.path.join(path, "BOUT.dmp."+str(i)+"."+informat)
        outfname = os.path.join(output, "BOUT.restart."+str(i)+"."+outformat)

        print((infname, " -> ", outfname))

        infile = DataFile(infname)
        outfile = DataFile(outfname, create=True)

        # Get the data always needed in restart files
        hist_hi = infile.read("iteration")
        print(("hist_hi = ", hist_hi))
        outfile.write("hist_hi", hist_hi)

        t_array = infile.read("t_array")
        tt = t_array[final]
        print(("tt = ", tt))
        outfile.write("tt", tt)

        tind = final
        if tind < 0.0:
            tind = len(t_array) + final

        NXPE = infile.read("NXPE")
        NYPE = infile.read("NYPE")
        print(("NXPE = ", NXPE, " NYPE = ", NYPE))
        outfile.write("NXPE", NXPE)
        outfile.write("NYPE", NYPE)

        # Get a list of variables
        varnames = infile.list()

        for var in varnames:
            if infile.ndims(var) == 4:
                # Could be an evolving variable

                print((" -> ", var))

                data = infile.read(var)

                if averagelast == 1:
                    slice = data[final, :, :, :]
                else:
                    slice = mean(data[(final - averagelast)
                                 :final, :, :, :], axis=0)

                print(slice.shape)

                outfile.write(var, slice)

        infile.close()
        outfile.close()


def redistribute(npes, path="data", nxpe=None, output=".", informat=None, outformat=None, mxg=2, myg=2):
    """Resize restart files across NPES processors.

    Does not check if new processor arrangement is compatible with the
    branch cuts. In this respect :py:func:`restart.split` is
    safer. However, BOUT++ checks the topology during initialisation
    anyway so this is not too serious.

    Parameters
    ----------
    npes : int
        Number of processors for the new restart files
    path : str, optional
        Path to original restart files (default: "data")
    nxpe : int, optional
        Number of processors to use in the x-direction (determines
        split: npes = nxpe * nype). Default is None which uses the
        same algorithm as BoutMesh (but without topology information)
        to determine a suitable value for nxpe.
    output : str, optional
        Location to save new restart files (default: current directory)
    informat : str, optional
        Specify file format of old restart files (must be a suffix
        understood by DataFile, e.g. 'nc'). Default uses the format of
        the first 'BOUT.restart.*' file listed by glob.glob.
    outformat : str, optional
        Specify file format of new restart files (must be a suffix
        understood by DataFile, e.g. 'nc'). Default is to use the same
        as informat.

    Returns
    -------
    True on success

    TODO
    ----
    - Replace printing errors with raising `ValueError`

    """

    if npes <= 0:
        print("ERROR: Negative or zero number of processors")
        return False

    if path == output:
        print("ERROR: Can't overwrite restart files")
        return False

    if informat is None:
        file_list = glob.glob(os.path.join(path, "BOUT.restart.*"))
    else:
        file_list = glob.glob(os.path.join(path, "BOUT.restart.*."+informat))

    nfiles = len(file_list)

    # Read old processor layout
    f = DataFile(file_list[0])

    # Get list of variables
    var_list = f.list()
    if len(var_list) == 0:
        print("ERROR: No data found")
        return False

    old_processor_layout = get_processor_layout(f, has_t_dimension=False)
    print("Grid sizes: ", old_processor_layout.nx,
          old_processor_layout.ny, old_processor_layout.mz)

    if nfiles != old_processor_layout.npes:
        print("WARNING: Number of restart files inconsistent with NPES")
        print("Setting nfiles = " + str(old_processor_layout.npes))
        nfiles = old_processor_layout.npes

    if nfiles == 0:
        print("ERROR: No restart files found")
        return False

    informat = file_list[0].split(".")[-1]
    if outformat is None:
        outformat = informat

    try:
        new_processor_layout = create_processor_layout(
            old_processor_layout, npes, nxpe=nxpe)
    except ValueError as e:
        print("Could not find valid processor split. " + e.what())

    nx = old_processor_layout.nx
    ny = old_processor_layout.ny
    mz = old_processor_layout.mz
    mxg = old_processor_layout.mxg
    myg = old_processor_layout.myg
    old_npes = old_processor_layout.npes
    old_nxpe = old_processor_layout.nxpe
    old_nype = old_processor_layout.nype
    old_mxsub = old_processor_layout.mxsub
    old_mysub = old_processor_layout.mysub

    nxpe = new_processor_layout.nxpe
    nype = new_processor_layout.nype
    mxsub = new_processor_layout.mxsub
    mysub = new_processor_layout.mysub
    mzsub = new_processor_layout.mz

    outfile_list = []
    for i in range(npes):
        outpath = os.path.join(output, "BOUT.restart."+str(i)+"."+outformat)
        outfile_list.append(DataFile(outpath, write=True, create=True))

    DataFileCache = create_cache(path, "BOUT.restart")

    for v in var_list:
        dimensions = f.dimensions(v)
        ndims = len(dimensions)

        # collect data
        data = collect(v, xguards=True, yguards=True, info=False,
                datafile_cache=DataFileCache)

        # write data
        for i in range(npes):
            ix = i % nxpe
            iy = int(i/nxpe)
            outfile = outfile_list[i]
            if v == "NPES":
                outfile.write(v, npes)
            elif v == "NXPE":
                outfile.write(v, nxpe)
            elif v == "NYPE":
                outfile.write(v, nype)
            elif v == "MXSUB":
                outfile.write(v, mxsub)
            elif v == "MYSUB":
                outfile.write(v, mysub)
            elif v == "MZSUB":
                outfile.write(v, mzsub)
            elif dimensions == ():
                # scalar
                outfile.write(v, data)
            elif dimensions == ('x', 'y'):
                # Field2D
                outfile.write(
                    v, data[ix*mxsub:(ix+1)*mxsub+2*mxg, iy*mysub:(iy+1)*mysub+2*myg])
            elif dimensions == ('x', 'z'):
                # FieldPerp
                yindex_global = data.attributes['yindex_global']
                if yindex_global + myg >= iy*mysub and yindex_global + myg < (iy+1)*mysub+2*myg:
                    outfile.write(v, data[ix*mxsub:(ix+1)*mxsub+2*mxg, :])
                else:
                    nullarray = BoutArray(np.zeros([mxsub+2*mxg, mysub+2*myg]), attributes={"bout_type":"FieldPerp", "yindex_global":-myg-1})
                    outfile.write(v, nullarray)
            elif dimensions == ('x', 'y', 'z'):
                # Field3D
                outfile.write(
                    v, data[ix*mxsub:(ix+1)*mxsub+2*mxg, iy*mysub:(iy+1)*mysub+2*myg, :])
            else:
                print(
                    "ERROR: variable found with unexpected dimensions,", dimensions, v)

    f.close()
    for outfile in outfile_list:
        outfile.close()

    return True


def resizeY(newy, path="data", output=".", informat="nc", outformat=None, myg=2):
    """Increase the number of Y points in restart files

    NOTE:
        * Can't overwrite

    Parameters
    ----------
    newy : int
        ny for the new file
    path : str, optional
        Path to original restart files (default: "data")
    output : str, optional
        Path to write new restart files (default: current directory)
    informat : str, optional
        File extension of original files (default: "nc")
    outformat : str, optional
        File extension of new files (default: use the same as `informat`)
    myg : int, optional
        Number of ghost points in y (default: 2)

    Returns
    -------
    True on success, else False

    TODO
    ----
    - Replace printing errors with raising `ValueError`
    - Make informat work like `redistribute`

    """

    if outformat is None:
        outformat = informat

    file_list = glob.glob(os.path.join(path, "BOUT.restart.*."+informat))

    nfiles = len(file_list)

    if nfiles == 0:
        print("ERROR: No restart files found")
        return False

    for i in range(nfiles):
        # Open each data file
        infname = os.path.join(path, "BOUT.restart."+str(i)+"."+informat)
        outfname = os.path.join(output, "BOUT.restart."+str(i)+"."+outformat)

        print("Processing %s -> %s" % (infname, outfname))

        infile = DataFile(infname)
        outfile = DataFile(outfname, create=True)

        # Copy basic information
        for var in ["hist_hi", "NXPE", "NYPE", "tt"]:
            data = infile.read(var)
            try:
                # Convert to scalar if necessary
                data = data[0]
            except:
                pass
            outfile.write(var, data)

        # Get a list of variables
        varnames = infile.list()

        for var in varnames:
            if infile.ndims(var) == 3:
                # Could be an evolving variable [x,y,z]

                print(" -> Resizing " + var)

                # Read variable from input
                indata = infile.read(var)

                nx, ny, nz = indata.shape

                # y coordinate in input and output data
                iny = (arange(ny) - myg + 0.5) / (ny - 2*myg)
                outy = (arange(newy) - myg + 0.5) / (newy - 2*myg)

                outdata = zeros([nx, newy, nz])

                for x in range(nx):
                    for z in range(nz):
                        f = interp1d(
                            iny, indata[x, :, z], bounds_error=False, fill_value=0.0)
                        outdata[x, :, z] = f(outy)

                outfile.write(var, outdata)
            elif infile.ndims(var) == 2:
                # Assume evolving variable [x,y]
                print(" -> Resizing " + var)

                # Read variable from input
                indata = infile.read(var)

                nx, ny = indata.shape

                # y coordinate in input and output data
                iny = (arange(ny) - myg + 0.5) / (ny - 2*myg)
                outy = (arange(newy) - myg + 0.5) / (newy - 2*myg)

                outdata = zeros([nx, newy])

                for x in range(nx):
                    f = interp1d(iny, indata[x, :],
                                 bounds_error=False, fill_value=0.0)
                    outdata[x, :] = f(outy)

                outfile.write(var, outdata)
            else:
                # Copy variable
                print(" -> Copying " + var)

                # Read variable from input
                data = infile.read(var)
                try:
                    # Convert to scalar if necessary
                    data = data[0]
                except:
                    pass
                outfile.write(var, data)

        infile.close()
        outfile.close()


def addvar(var, value, path="."):
    """Adds a variable with constant value to all restart files.

    .. warning:: Modifies restart files in place! This is in contrast
                 to most of the functions in this module!

    This is useful for restarting simulations whilst turning on new
    equations. By default BOUT++ throws an error if an evolving
    variable is not in the restart file. By setting an option the
    variable can be set to zero. This allows it to start with a
    non-zero value.

    Parameters
    ----------
    var : str
        The name of the variable to add
    value : float
        Constant value for the variable
    path : str, optional
        Input path to data files (default: current directory)

    """

    file_list = glob.glob(os.path.join(path, "BOUT.restart.*"))
    nfiles = len(file_list)

    print("Number of restart files: %d" % (nfiles,))
    # Loop through all the restart files
    for filename in file_list:
        print(filename)
        # Open the restart file for writing (modification)
        with DataFile(filename, write=True) as df:
            size = None
            # Find a 3D variable and get its size
            for varname in df.list():
                size = df.size(varname)
                if len(size) == 3:
                    break
            if size is None:
                raise Exception("no 3D variables found")

            # Create a new 3D array with input value
            data = np.zeros(size) + value

            # Set the variable in the NetCDF file
            df.write(var, data)
