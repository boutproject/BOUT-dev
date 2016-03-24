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

import numpy
from numpy import mean, zeros, arange
from math import sqrt
from numpy.random import normal

from scipy.interpolate import interp1d

try:
    import os
    import sys
    import glob
except ImportError:
    raise ImportError("ERROR: os, sys or glob modules not available")

def split(nxpe, nype, path="data", output="./", informat="nc", outformat=None):
    """Split restart files across NXPE x NYPE processors.

    Returns True on success
    """

    if outformat == None:
        outformat = informat

    mxg = 2
    myg = 2

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

    old_nype =int(old_npes/old_nxpe)

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
        pey = int(mype/ nxpe)

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

def expand(newz, path="data", output="./", informat="nc", outformat=None):
    """Increase the number of Z points in restart files

    """
    if outformat == None:
        outformat = informat

    if path == output:
        print("ERROR: Can't overwrite restart files when expanding")
        return False

    def is_pow2(x):
        """Returns true if x is a power of 2"""
        return (x > 0) and ((x & (x-1)) == 0)

    if not is_pow2(newz-1):
        print("ERROR: New Z size must be a power of 2 + 1")
        return False

    file_list = glob.glob(os.path.join(path, "BOUT.restart.*."+informat))
    nfiles = len(file_list)

    # Get the file extension
    ind = file_list[0].rfind(".")


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
            if var == None:
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

    if outformat == None:
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

    if npes <= 0:
        print("ERROR: Negative or zero number of processors")
        return False

    if path == output:
        print("ERROR: Can't overwrite restart files")
        return False

    if informat == None:
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

    old_npes = f.read('NPES')
    old_nxpe = f.read('NXPE')
    old_nype = int(old_npes/old_nxpe)

    if nfiles != old_npes:
        print("WARNING: Number of restart files inconsistent with NPES")
        print("Setting nfiles = " + str(old_npes))
        nfiles = old_npes

    if nfiles == 0:
        print("ERROR: No restart files found")
        return False

    informat = file_list[0].split(".")[-1]
    if outformat == None:
        outformat = informat

    old_mxsub = 0
    old_mysub = 0
    mz = 0

    for v in var_list:
        if f.ndims(v) == 3:
            s = f.size(v)
            old_mxsub = s[0] - 2*mxg
            if old_mxsub < 0:
                if s[0] == 1:
                    old_mxsub = 1
                    mxg = 0
                elif s[0] == 3:
                    old_mxsub = 1
                    mxg = 1
                else:
                    print("Number of x points is wrong?")
                    return False

            old_mysub = s[1] - 2*myg
            if old_mysub < 0:
                if s[1] == 1:
                    old_mysub = 1
                    myg = 0
                elif s[1] == 3:
                    old_mysub = 1
                    myg = 1
                else:
                    print("Number of y points is wrong?")
                    return False

            mz = s[2]
            break

    # Calculate total size of the grid
    nx = old_mxsub * old_nxpe
    ny = old_mysub * old_nype
    print("Grid sizes: ", nx, ny, mz)

    if nxpe == None: # Copy algorithm from BoutMesh for selecting nxpe
        ideal = sqrt(float(nx) * float(npes) / float(ny)) # Results in square domain

        for i in range(1,npes+1):
            if npes%i == 0 and nx%i == 0 and int(nx/i) >= mxg and ny%(npes/i) == 0:
                # Found an acceptable value
                # Warning: does not check branch cuts!

                if nxpe==None or abs(ideal - i) < abs(ideal - nxpe):
                    nxpe = i # Keep value nearest to the ideal

        if nxpe == None:
            print("ERROR: could not find a valid value for nxpe")
            return False

    nype = int(npes/nxpe)

    outfile_list = []
    for i in range(npes):
        outpath = os.path.join(output, "BOUT.restart."+str(i)+"."+outformat)
        outfile_list.append(DataFile(outpath, write=True, create=True))
    infile_list = []
    for i in range(old_npes):
        inpath = os.path.join(path, "BOUT.restart."+str(i)+"."+outformat)
        infile_list.append(DataFile(inpath))

    old_mxsub = int(nx/old_nxpe)
    old_mysub = int(ny/old_nype)
    mxsub = int(nx/nxpe)
    mysub = int(ny/nype)
    for v in var_list:
          ndims = f.ndims(v)

          #collect data
          if ndims == 0:
              #scalar
              data = f.read(v)
          elif ndims == 2:
              data = numpy.zeros( (nx+2*mxg,ny+2*nyg) )
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
              data = numpy.zeros( (nx+2*mxg,ny+2*myg,mz) )
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
                  data[ix*old_mxsub+ixstart:(ix+1)*old_mxsub+2*mxg+ixend, iy*old_mysub+iystart:(iy+1)*old_mysub+2*myg+iyend, :] = infile_list[i].read(v)[ixstart:old_mxsub+2*mxg+ixend, iystart:old_mysub+2*myg+iyend, :]
          else:
              print("ERROR: variable found with unexpected number of dimensions,",ndims,v)
              return False

          # write data
          for i in range(npes):
              ix = i%nxpe
              iy = int(i/nxpe)
              outfile = outfile_list[i]
              if v == "NPES":
                  outfile.write(v,npes)
              elif v == "NXPE":
                  outfile.write(v,nxpe)
              elif ndims == 0:
                  # scalar
                  outfile.write(v,data)
              elif ndims == 2:
                  # Field2D
                  outfile.write(v,data[ix*mxsub:(ix+1)*mxsub+2*mxg, iy*mysub:(iy+1)*mysub+2*myg])
              elif ndims == 3:
                  # Field3D
                  outfile.write(v,data[ix*mxsub:(ix+1)*mxsub+2*mxg, iy*mysub:(iy+1)*mysub+2*myg, :])
              else:
                  print("ERROR: variable found with unexpected number of dimensions,",f.ndims(v))

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

    if outformat == None:
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

        print("Processing %s -> %s", infname, outfname)

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

                print(" -> " + var)

                # Read variable from input
                indata = infile.read(var)

                nx,ny,nz = indata.shape

                # y coordinate in input and output data
                iny = (arange(ny) - myg + 0.5) / (ny - 2*myg)
                outy = (arange(newy) - myg + 0.5) / (newy - 2*myg)

                outdata = zeros([nx, newy, nz])

                for x in range(nx):
                    for z in range(nz):
                        f = interp1d(iny, indata[x,:,z], bounds_error=False, fill_value=0.0)
                        outdata[x,:,z] = f(outy)

                outfile.write(var, outdata)
        infile.close()
        outfile.close()
