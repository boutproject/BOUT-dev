"""Routines for redistributing files over different numbers of
processors

"""

from math import sqrt
from collections import namedtuple

processor_layout = namedtuple("BOUT_processor_layout", ["nxpe", "nype", "npes", "mxsub", "mysub", "nx", "ny", "mz", "mxg", "myg"])


def get_processor_layout(boutfile, has_t_dimension=True, mxg=2, myg=2):
    """Given a BOUT.restart.* or BOUT.dmp.* file (as a DataFile object),
    return the processor layout for its data

    Parameters
    ----------
    boutfile : DataFile
        Restart or dump file to read
    has_t_dimension : bool, optional
        Does this file have a time dimension?
    mxg, myg : int, optional
        Number of x, y guard cells

    Returns
    -------
    processor_layout : (int, int, int, int, int, int, int, int, int, int)
        A namedtuple containing the number of processors in x and y;
        the total number of procesors; the size of the grid in x and y
        on a single processor; the total size of the grid in x, y and
        z; and the number of guard cells in x and y

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
            mxsub = s[maxdims - 3] - 2 * mxg
            if mxsub < 0:
                if s[maxdims - 3] == 1:
                    mxsub = 1
                    mxg = 0
                elif s[maxdims - 3] == 3:
                    mxsub = 1
                    mxg = 1
                else:
                    print("Number of x points is wrong?")
                    return False

            mysub = s[maxdims - 2] - 2 * myg
            if mysub < 0:
                if s[maxdims - 2] == 1:
                    mysub = 1
                    myg = 0
                elif s[maxdims - 2] == 3:
                    mysub = 1
                    myg = 1
                else:
                    print("Number of y points is wrong?")
                    return False

            mz = s[maxdims - 1]
            break

    # Calculate total size of the grid
    nx = mxsub * nxpe
    ny = mysub * nype

    result = processor_layout(nxpe=nxpe, nype=nype, npes=npes, mxsub=mxsub, mysub=mysub, nx=nx, ny=ny, mz=mz, mxg=mxg, myg=myg)

    return result


def create_processor_layout(old_processor_layout, npes, nxpe=None):
    """Convert one processor layout into another one with a different
    total number of processors

    If nxpe is None, use algorithm from BoutMesh to select optimal nxpe.
    Otherwise, check nxpe is valid (divides npes)

    Parameters
    ----------
    old_processor_layout : processor_layout
        The processor layout to convert
    npes : int
        The new total number of procesors
    nxpe : int, optional
        The number of procesors in x to use

    Returns
    -------
    processor_layout : (int, int, int, int, int, int, int, int, int, int)
        A namedtuple containing the number of processors in x and y;
        the total number of procesors; the size of the grid in x and y
        on a single processor; the total size of the grid in x, y and
        z; and the number of guard cells in x and y

    """

    if nxpe is None:  # Copy algorithm from BoutMesh for selecting nxpe
        ideal = sqrt(float(old_processor_layout.nx) * float(npes) / float(old_processor_layout.ny))
                     # Results in square domain

        for i in range(1, npes + 1):
            if npes % i == 0 and old_processor_layout.nx % i == 0 and int(old_processor_layout.nx / i) >= old_processor_layout.mxg and old_processor_layout.ny % (npes / i) == 0:
                # Found an acceptable value
                # Warning: does not check branch cuts!

                if nxpe is None or abs(ideal - i) < abs(ideal - nxpe):
                    nxpe = i  # Keep value nearest to the ideal

        if nxpe is None:
            raise ValueError("ERROR: could not find a valid value for nxpe")
    elif npes % nxpe != 0:
        raise ValueError(
            "ERROR: requested nxpe is invalid, it does not divide npes")

    nype = int(npes / nxpe)

    mxsub = int(old_processor_layout.nx / nxpe)
    mysub = int(old_processor_layout.ny / nype)

    result = processor_layout(nxpe=nxpe, nype=nype, npes=npes, mxsub=mxsub, mysub=mysub, nx=old_processor_layout.nx, ny=old_processor_layout.ny, mz=old_processor_layout.mz, mxg=old_processor_layout.mxg, myg=old_processor_layout.myg)

    return result
