# Routines for manipulating restart files

try:
    from boututils import DataFile
except ImportError:
    print "ERROR: restart module needs DataFile"
    raise

from numpy import mean

try:
    import os
    import sys
    import glob
except ImportError:
    print "ERROR: os, sys or glob modules not available"
    raise

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
        print "ERROR: Negative or zero number of processors"
        return False
    
    if path == output:
        print "ERROR: Can't overwrite restart files"
        return False

    file_list = glob.glob(os.path.join(path, "BOUT.restart.*."+informat))
    nfiles = len(file_list)

    if nfiles == 0:
        print "ERROR: No restart files found"
        return False

    # Read old processor layout
    f = DataFile(os.path.join(path, file_list[0]))

    # Get list of variables
    var_list = f.list()
    if len(var_list) == 0:
        print "ERROR: No data found"
        return False
    
    old_npes = f.read('NPES')
    old_nxpe = f.read('NXPE')

    if nfiles != old_npes:
        print "WARNING: Number of restart files inconsistent with NPES"
        print "Setting nfiles = " + str(old_npes)
        nfiles = old_npes

    if old_npes % old_nxpe != 0:
        print "ERROR: Old NPES is not a multiple of old NXPE"
        return False

    old_nype = old_npes / old_nxpe

    if nype % old_nype != 0:
        print "SORRY: New nype must be a multiple of old nype"
        return False

    if nxpe % old_nxpe != 0:
        print "SORRY: New nxpe must be a multiple of old nxpe"
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
    print "Grid sizes: ", nx, ny, mz
    
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

        print "New: "+str(mype)+" ("+str(pex)+", "+str(pey)+")"
        print " =>  "+str(old_mype)+" ("+str(old_pex)+", "+str(old_pey)+") : ("+str(old_x)+", "+str(old_y)+")"

        # 

def expand(newz, path="data", output="./", informat="nc", outformat=None):
    """Increase the number of Z points in restart files

    """
    if outformat == None:
        outformat = informat
    
    if path == output:
        print "ERROR: Can't overwrite restart files when expanding"
        return False
    
    def is_pow2(x):
        """Returns true if x is a power of 2"""
        return (x > 0) and ((x & (x-1)) == 0)
    
    if not is_pow2(newz-1):
        print "ERROR: New Z size must be a power of 2 + 1"
        return False
    
    file_list = glob.glob(os.path.join(path, "BOUT.restart.*."+informat))
    nfiles = len(file_list)
    
    # Get the file extension
    ind = file_list[0].rfind(".")
    

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
    
    print "Number of data files: ", nfiles

    for i in range(nfiles):
        # Open each data file
        infname  = os.path.join(path, "BOUT.dmp."+str(i)+"."+informat)
        outfname = os.path.join(output, "BOUT.restart."+str(i)+"."+outformat)

        print infname, " -> ", outfname
        
        infile = DataFile(infname)
        outfile = DataFile(outfname, create=True)

        # Get the data always needed in restart files
        hist_hi = infile.read("iteration")
        print "hist_hi = ", hist_hi
        outfile.write("hist_hi", hist_hi)
        
        t_array = infile.read("t_array")
        tt = t_array[final]
        print "tt = ", tt
        outfile.write("tt", tt)

        NXPE = infile.read("NXPE")
        NYPE = infile.read("NYPE")
        NPES = NXPE * NYPE
        print "NPES = ", NPES, " NXPE = ", NXPE
        outfile.write("NPES", NPES)
        outfile.write("NXPE", NXPE)
        
        # Get a list of variables
        varnames = infile.list()

        for var in varnames:
            if infile.ndims(var) == 4:
                # Could be an evolving variable

                print " -> ", var
                
                data = infile.read(var)

                if averagelast == 1:
                    slice = data[final,:,:,:]
                else:
                    slice = mean(data[(final - averagelast):final,:,:,:], axis=0)

                print slice.shape
                
                outfile.write(var, slice)
        
        infile.close()
        outfile.close()
