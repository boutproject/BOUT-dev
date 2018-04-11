"""Routines for redistributing files over different numbers of processors"""

from math import sqrt

def get_processor_layout(boutfile, has_t_dimension=True, mxg=2, myg=2):
    """
    Given a BOUT.restart.* or BOUT.dmp.* file (as a DataFile object), return the processor layout for its data
    """
    nxpe = boutfile.read('NXPE')
    nype = boutfile.read("NYPE")
    npes = nxpe * nype

    # Get list of variables
    var_list = boutfile.list()
    if len(var_list) == 0:
        raise ValueError("ERROR: No data found")

    mxsub = 0
    mysub = 0
    mz = 0
    if has_t_dimension:
        maxdims = 4
    else:
        maxdims = 3
    for v in var_list:
        if boutfile.ndims(v) == maxdims:
            s = boutfile.size(v)
            mxsub = s[maxdims-3] - 2*mxg
            if mxsub < 0:
                if s[maxdims-3] == 1:
                    mxsub = 1
                    mxg = 0
                elif s[maxdims-3] == 3:
                    mxsub = 1
                    mxg = 1
                else:
                    print("Number of x points is wrong?")
                    return False

            mysub = s[maxdims-2] - 2*myg
            if mysub < 0:
                if s[maxdims-2] == 1:
                    mysub = 1
                    myg = 0
                elif s[maxdims-2] == 3:
                    mysub = 1
                    myg = 1
                else:
                    print("Number of y points is wrong?")
                    return False

            mz = s[maxdims-1]
            break

    # Calculate total size of the grid
    nx = mxsub * nxpe
    ny = mysub * nype

    return nxpe, nype, npes, mxsub, mysub, nx, ny, mz, mxg, myg

def create_processor_layout(npes, nx, ny, nxpe=None, mxg=2, myg=2):
    """
    If nxpe==None, use algorithm from BoutMesh to select optimal nxpe.
    Otherwise, check nxpe is valid (divides npes)
    """

    if nxpe is None: # Copy algorithm from BoutMesh for selecting nxpe
        ideal = sqrt(float(nx) * float(npes) / float(ny)) # Results in square domain

        for i in range(1,npes+1):
            if npes%i == 0 and nx%i == 0 and int(nx/i) >= mxg and ny%(npes/i) == 0:
                # Found an acceptable value
                # Warning: does not check branch cuts!

                if nxpe is None or abs(ideal - i) < abs(ideal - nxpe):
                    nxpe = i # Keep value nearest to the ideal

        if nxpe is None:
            raise ValueError("ERROR: could not find a valid value for nxpe")
    elif npes%nxpe != 0:
        raise ValueError("ERROR: requested nxpe is invalid, it does not divide npes")

    nype = int(npes/nxpe)

    mxsub = int(nx/nxpe)
    mysub = int(ny/nype)

    return nxpe, nype, mxsub, mysub
