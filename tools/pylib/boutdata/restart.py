"""Routines for manipulating restart files"""

from __future__ import print_function
from __future__ import division
try:
    from builtins import str
    from builtins import range
except:
    pass

try:
    from boututils.datafile import DataFile
except ImportError:
    raise ImportError("ERROR: restart module needs DataFile")

import multiprocessing
import numpy as np
from numpy import mean, zeros, arange
from math import sqrt
from numpy.random import normal

from scipy.interpolate import interp1d
try:
    from scipy.interpolate import RegularGridInterpolator
except:
    pass


try:
    import os
    import sys
    import glob
except ImportError:
    raise ImportError("ERROR: os, sys or glob modules not available")

def split(nxpe, nype, path="data", output="./", informat="nc", outformat=None, mxg=2, myg=2):
    """Split restart files across NXPE x NYPE processors.

    Returns True on success
    """

    if outformat is None:
        outformat = informat

    npes = nxpe * nype

    if npes <= 0:
        print("ERROR: Negative or zero number of processors")
        return False

    if path == output:
        print("ERROR: Can't overwrite restart files")
        return False

    file_list = glob.glob(os.path.join(path, "BOUT.restart.*."+informat))
    nfiles = len(file_list)

    if nfiles == 0:
        print("ERROR: No restart files found")
        return False

    # Read old processor layout
    f = DataFile(os.path.join(path, file_list[0]))

    # Get list of variables
    var_list = f.list()
    if len(var_list) == 0:
        print("ERROR: No data found")
        return False

    old_npes = f.read('NPES')
    old_nxpe = f.read('NXPE')

    if nfiles != old_npes:
        print("WARNING: Number of restart files inconsistent with NPES")
        print("Setting nfiles = " + str(old_npes))
        nfiles = old_npes

    if old_npes % old_nxpe != 0:
        print("ERROR: Old NPES is not a multiple of old NXPE")
        return False

    old_nype = int(old_npes/old_nxpe)

    if nype % old_nype != 0:
        print("SORRY: New nype must be a multiple of old nype")
        return False

    if nxpe % old_nxpe != 0:
        print("SORRY: New nxpe must be a multiple of old nxpe")
        return False

    # Get dimension sizes

    old_mxsub = 0
    old_mysub = 0
    mz = 0

    for v in var_list:
        if f.ndims(v) == 3:
            s = f.size(v)
            old_mxsub = s[0] - 2*mxg
            old_mysub = s[1] - 2*myg
            mz = s[2]
            break

    f.close()

    # Calculate total size of the grid
    nx = old_mxsub * old_nxpe
    ny = old_mysub * old_nype
    print(("Grid sizes: ", nx, ny, mz))

    # Create the new restart files
    for mype in range(npes):
        # Calculate X and Y processor numbers
        pex = mype % nxpe
        pey = int(mype / nxpe)

        old_pex = int(pex / xs)
        old_pey = int(pey / ys)

        old_x = pex % xs
        old_y = pey % ys

        # Old restart file number
        old_mype = old_nxpe * old_pey + old_pex

        # Calculate indices in old restart file
        xmin = old_x*mxsub
        xmax = xmin + mxsub - 1 + 2*mxg
        ymin = old_y*mysub
        ymax = ymin + mysub - 1 + 2*myg

        print("New: "+str(mype)+" ("+str(pex)+", "+str(pey)+")")
        print(" =>  "+str(old_mype)+" ("+str(old_pex)+", "+str(old_pey)+") : ("+str(old_x)+", "+str(old_y)+")")

        #

def resize3DField(var, data, coordsAndSizesTuple, method, mute):
    """
    Resizing of the 3D fields.

    To be called by resize.

    Written as a function in order to call it using multiprocesse. Must
    be defined as a top level function in order to be pickable by the
    multiprocess.

    See the function resize for details.
    """

    # Unpack the tuple for better readability
    xCoordOld, yCoordOld, zCoordOld,\
    xCoordNew, yCoordNew, zCoordNew,\
    newNx, newNy, newNz = coordsAndSizesTuple

    if not(mute):
        print("    Resizing "+var + ' to (nx,ny,nz) = ({},{},{})'.format(newNx, newNy, newNz) )


    # Make the regular grid function (see examples in
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html
    # for details)
    gridInterpolator = RegularGridInterpolator((xCoordOld, yCoordOld, zCoordOld), data, method)

    # Need to fill with one exrta z plane (will only contain zeros)
    newData = np.zeros((newNx, newNy, newNz))

    # Interpolate to the new values
    for xInd, x in enumerate(xCoordNew):
        for yInd, y in enumerate(yCoordNew):
            for zInd, z in enumerate(zCoordNew):
                newData[xInd, yInd, zInd] = gridInterpolator([x, y, z])

    return var, newData



def resize(newNx, newNy, newNz, mxg=2, myg=2,\
           path="data", output="./", informat="nc", outformat=None,\
           method='linear', maxProc=None, mute=False):
    """
    Increase/decrease the number of points in restart files.

    NOTE: Can't overwrite
    WARNING: Currently only implemented with uniform BOUT++ grid

    Parameters
    -----
    newNx : int
        nx for the new file (including ghost points)
    newNy : int
        ny for the new file (including ghost points)
    newNz : int
        nz for the new file (including last unused z-plane)
    mxg : int
        Number of ghost points in x. **NOTE:** Default is 2
    myg : int
        Number of ghost points in y. **NOTE:** Default is 2
    path : str
        Input path
    output : str
        Output path
    informat : str
        File extension of input
    outformat : [None|str]
        File extension of output
    method : ['linear'|'nearest']
        What interpolation method to be used
    maxProc: [None|int]
        Limits maximum processors to use when interpolating if set
    mute : [True|False]
        Whether or not output should be printed from this function

    Returns
    -------
    return : [True|False]
        True on success, else False

    Todo
    ----
    Add 2D field interpolation
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
            coordsAndSizesTuple = (xCoordOld, yCoordOld, zCoordOld,\
                                   xCoordNew, yCoordNew, zCoordNew,\
                                   newNx, newNy, newNz)

            # Loop over the variables in the old file
            for var in old.list():
                # Read the data
                data = old.read(var)

                # Find 3D variables
                if old.ndims(var) == 3:

                    # Asynchronous call (locks first at .get())
                    jobs.append(pool.apply_async(resize3DField,\
                               args = (var, data, coordsAndSizesTuple, method, mute, )))

                else:
                    if not(mute):
                        print("    Copying "+var)
                        newData = data.copy()
                    if not(mute):
                        print("Writing "+var)
                    new.write(var, newData)

            for job in jobs:
                var, newData = job.get()
                if not(mute):
                    print("Writing "+var)
                new.write(var, newData)

            # Close the pool of workers
            pool.close()
            # Wait for all processes to finish
            pool.join()

    return True



def resizeZ(newNz, path="data", output="./", informat="nc", outformat=None):
    """
    Increase the number of Z points in restart files.

    The python equivalent of ../../idllib/expand_restarts.pro

    NOTE:
        * Can't overwrite
        * Will not yield a result close to the original if there are
          aymmetires in the z-direction

    Input
    -----
    path       Input path
    output     Output path
    informat   File extension of input
    outformat  File extension of output

    Returns
    -------
    True on success, else False
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

                # Find 3D variables
                if old.ndims(var) == 3:
                    print("    Resizing "+var)

                    nx, ny, nz = data.shape

                    newdata = np.zeros((nx, ny, newNz))
                    for x in range(nx):
                        for y in range(ny):
                            f_old = np.fft.fft(data[x, y,:])

                            # Number of points in f is power of 2
                            f_new = np.zeros(newNz )

                            # Copy coefficients across (ignoring Nyquist)
                            f_new[0] = f_old[0] # DC
                            for m in range(1, int(nz/2)):
                                # + ve frequencies
                                f_new[m] = f_old[m]
                                # - ve frequencies
                                f_new[newNz-m] = f_old[nz-m]

                            # Invert fft
                            newdata[x, y,:] = np.fft.ifft(f_new).real
                            newdata[x, y,:] = newdata[x, y, 0]

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

                new.write(var, newdata)

    return True



def addnoise(path=".", var=None, scale=1e-5):
    """
    Add random noise to restart files

    Inputs
    ------

    path   Path to the restart files
    var    The variable to modify. By default all 3D variables are modified
    scale  Amplitude of the noise. Gaussian noise is used, with zero mean
           and this parameter as the standard deviation

    Returns
    -------
    None
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
                        data = d.read(v)
                        data += normal(scale=scale, size=data.shape)
                        d.write(v, data)
            else:
                # Modify a single variable
                print(" -> "+var)
                data = d.read(var)
                data += normal(scale=scale, size=data.shape)
                d.write(var, data)

def scalevar(var, factor, path="."):
    """
    Scales a variable by a given factor, modifying
    restart files in place

    Inputs
    ------

    var      Name of the variable  (string)
    factor   Factor to multiply    (float)
    path     Path to the restart files

    Returns
    -------
    None
    """

    file_list = glob.glob(os.path.join(path, "BOUT.restart.*"))
    nfiles = len(file_list)

    print("Number of restart files: %d" % (nfiles,))
    for file in file_list:
        print(file)
        with DataFile(file, write=True) as d:
            d[var] = d[var] * factor



def create(averagelast=1, final=-1, path="data", output="./", informat="nc", outformat=None):
    """
    Create restart files from data (dmp) files.

    Inputs
    ======

    averagelast   Number of time points to average over.
                  Default is 1 i.e. just take last time-point

    final         The last time point to use. Default is last (-1)

    path          Path to the input data files

    output        Path where the output restart files should go

    informat      Format of the input data files

    outformat     Format of the output restart files

    """

    if outformat is None:
        outformat = informat

    file_list = glob.glob(os.path.join(path, "BOUT.dmp.*."+informat))
    nfiles = len(file_list)

    print(("Number of data files: ", nfiles))

    for i in range(nfiles):
        # Open each data file
        infname  = os.path.join(path, "BOUT.dmp."+str(i)+"."+informat)
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
        NPES = NXPE * NYPE
        print(("NPES = ", NPES, " NXPE = ", NXPE))
        outfile.write("NPES", NPES)
        outfile.write("NXPE", NXPE)

        # Get a list of variables
        varnames = infile.list()

        for var in varnames:
            if infile.ndims(var) == 4:
                # Could be an evolving variable

                print((" -> ", var))

                data = infile.read(var)

                if averagelast == 1:
                    slice = data[final,:,:,:]
                else:
                    slice = mean(data[(final - averagelast):final,:,:,:], axis=0)

                print(slice.shape)

                outfile.write(var, slice)

        infile.close()
        outfile.close()

def redistribute(npes, path="data", nxpe=None, output=".", informat=None, outformat=None, mxg=2, myg=2):
    """Resize restart files across NPES processors.

    Does not check if new processor arrangement is compatible with the branch cuts. In this respect restart.split is safer. However, BOUT++ checks the topology during initialisation anyway so this is not too serious.

    Parameters
    ----------
    npes : int
        number of processors for the new restart files
    path : string, optional
        location of old restart files
    nxpe : int, optional
        number of processors to use in the x-direction (determines split: npes = nxpe * nype). Default is None which uses the same algorithm as BoutMesh (but without topology information) to determine a suitable value for nxpe.
    output : string, optional
        location to save new restart files
    informat : string, optional
        specify file format of old restart files (must be a suffix understood by DataFile, e.g. 'nc'). Default uses the format of the first 'BOUT.restart.*' file listed by glob.glob.
    outformat : string, optional
        specify file format of new restart files (must be a suffix understood by DataFile, e.g. 'nc'). Default is to use the same as informat.

    Returns
    -------
    True on success
    """

    try:
        from boutdata.processor_rearrange import get_processor_layout, create_processor_layout
    except ImportError:
        raise ImportError("ERROR: restart.redistribute needs boutdata.processor_rearrange")

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

    old_nxpe, old_nype, old_npes, old_mxsub, old_mysub, nx, ny, mz, mxg, myg = get_processor_layout(f, has_t_dimension=False)
    print("Grid sizes: ", nx, ny, mz)

    if nfiles != old_npes:
        print("WARNING: Number of restart files inconsistent with NPES")
        print("Setting nfiles = " + str(old_npes))
        nfiles = old_npes

    if nfiles == 0:
        print("ERROR: No restart files found")
        return False

    informat = file_list[0].split(".")[-1]
    if outformat is None:
        outformat = informat

    try:
        nxpe, nype, mxsub, mysub = create_processor_layout(npes, nx=nx, ny=ny, nxpe=nxpe, mxg=mxg, myg=myg)
    except ValueError as e:
        print("Could not find valid processor split. " + e.what())

    outfile_list = []
    for i in range(npes):
        outpath = os.path.join(output, "BOUT.restart."+str(i)+"."+outformat)
        outfile_list.append(DataFile(outpath, write=True, create=True))
    infile_list = []
    for i in range(old_npes):
        inpath = os.path.join(path, "BOUT.restart."+str(i)+"."+outformat)
        infile_list.append(DataFile(inpath))

    for v in var_list:
        ndims = f.ndims(v)

        # collect data
        if ndims == 0:
            # scalar
            data = f.read(v)
        elif ndims == 2:
            data = np.zeros( (nx+2*mxg, ny+2*myg) )
            for i in range(old_npes):
                ix = i%old_nxpe
                iy = int(i/old_nxpe)
                ixstart = mxg
                if ix == 0:
                    ixstart = 0
                ixend = -mxg
                if ix == old_nxpe-1:
                    ixend = 0
                iystart = myg
                if iy == 0:
                    iystart = 0
                iyend = -myg
                if iy == old_nype-1:
                    iyend = 0
                data[ix*old_mxsub+ixstart:(ix+1)*old_mxsub+2*mxg+ixend, iy*old_mysub+iystart:(iy+1)*old_mysub+2*myg+iyend] = infile_list[i].read(v)[ixstart:old_mxsub+2*mxg+ixend, iystart:old_mysub+2*myg+iyend]
        elif ndims == 3:
            data = np.zeros( (nx+2*mxg, ny+2*myg, mz) )
            for i in range(old_npes):
                ix = i%old_nxpe
                iy = int(i/old_nxpe)
                ixstart = mxg
                if ix == 0:
                    ixstart = 0
                ixend = -mxg
                if ix == old_nxpe-1:
                    ixend = 0
                iystart = myg
                if iy == 0:
                    iystart = 0
                iyend = -myg
                if iy == old_nype-1:
                    iyend = 0
                data[ix*old_mxsub+ixstart:(ix+1)*old_mxsub+2*mxg+ixend, iy*old_mysub+iystart:(iy+1)*old_mysub+2*myg+iyend,:] = infile_list[i].read(v)[ixstart:old_mxsub+2*mxg+ixend, iystart:old_mysub+2*myg+iyend,:]
        else:
            print("ERROR: variable found with unexpected number of dimensions,", ndims, v)
            return False

        # write data
        for i in range(npes):
            ix = i%nxpe
            iy = int(i/nxpe)
            outfile = outfile_list[i]
            if v == "NPES":
                outfile.write(v, npes)
            elif v == "NXPE":
                outfile.write(v, nxpe)
            elif v == "NYPE":
                outfile.write(v, nype)
            elif ndims == 0:
                # scalar
                outfile.write(v, data)
            elif ndims == 2:
                # Field2D
                outfile.write(v, data[ix*mxsub:(ix+1)*mxsub+2*mxg, iy*mysub:(iy+1)*mysub+2*myg])
            elif ndims == 3:
                # Field3D
                outfile.write(v, data[ix*mxsub:(ix+1)*mxsub+2*mxg, iy*mysub:(iy+1)*mysub+2*myg,:])
            else:
                print("ERROR: variable found with unexpected number of dimensions,", f.ndims(v))

    f.close()
    for infile in infile_list:
        infile.close()
    for outfile in outfile_list:
        outfile.close()

    return True

def resizeY(newy, path="data", output=".", informat="nc", outformat=None,myg=2):
    """
    Resize all the restart files in Y
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
        infname  = os.path.join(path, "BOUT.restart."+str(i)+"."+informat)
        outfname = os.path.join(output, "BOUT.restart."+str(i)+"."+outformat)

        print("Processing %s -> %s" % (infname, outfname))

        infile = DataFile(infname)
        outfile = DataFile(outfname, create=True)

        # Copy basic information
        for var in ["hist_hi", "NPES", "NXPE", "tt"]:
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
                        f = interp1d(iny, indata[x,:, z], bounds_error=False, fill_value=0.0)
                        outdata[x,:, z] = f(outy)

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
                    f = interp1d(iny, indata[x,:], bounds_error=False, fill_value=0.0)
                    outdata[x,:] = f(outy)

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

    This is useful for restarting simulations whilst turning on new
    equations. By default BOUT++ throws an error if an evolving
    variable is not in the restart file. By setting an option the
    variable can be set to zero. This allows it to start with a
    non-zero value.

    Input
    -----
    var      The variable to add
    value    Constant value for the variable
    path     Path to directory containing restart files
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
