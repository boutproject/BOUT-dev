# Generates an input mesh for circular, large aspect-ratio
# simulations:
#
# o Constant magnetic field
# o Curvature output as a 3D logB variable
# o Z is poloidal direction
# o Y is parallel (toroidal)
#
# NOTE: This reverses the standard BOUT/BOUT++ convention
#       so here Bt and Bp are reversed
#

from __future__ import division
from __future__ import print_function
from builtins import range

from numpy import zeros, ndarray, pi, cos, sin, outer, linspace, sqrt

from boututils.datafile import DataFile  # Wrapper around NetCDF4 libraries


def generate(
    nx,
    ny,
    R=2.0,
    r=0.2,  # Major & minor radius
    dr=0.05,  # Radial width of domain
    Bt=1.0,  # Toroidal magnetic field
    q=5.0,  # Safety factor
    mxg=2,
    file="circle.nc",
):

    # q = rBt / RBp
    Bp = r * Bt / (R * q)

    # Minor radius as function of x. Choose so boundary
    # is half-way between grid points

    h = dr / (nx - 2.0 * mxg)  # Grid spacing in r
    rminor = linspace(
        r - 0.5 * dr - (mxg - 0.5) * h, r + 0.5 * dr + (mxg - 0.5) * h, nx
    )

    # mesh spacing in x and y
    dx = ndarray([nx, ny])
    dx[:, :] = r * Bt * h  # NOTE: dx is toroidal flux

    dy = ndarray([nx, ny])
    dy[:, :] = 2.0 * pi / ny

    # LogB = log(1/(1+r/R cos(theta))) =(approx) -(r/R)*cos(theta)
    logB = zeros([nx, ny, 3])  # (constant, n=1 real, n=1 imag)

    # At y = 0, Rmaj = R + r*cos(theta)
    logB[:, 0, 1] = -(rminor / R)

    # Moving in y, phase shift by (toroidal angle) / q
    for y in range(1, ny):
        dtheta = y * 2.0 * pi / ny / q  # Change in poloidal angle

        logB[:, y, 1] = -(rminor / R) * cos(dtheta)
        logB[:, y, 2] = -(rminor / R) * sin(dtheta)

    # Shift angle from one end of y to the other
    ShiftAngle = ndarray([nx])
    ShiftAngle[:] = 2.0 * pi / q

    Rxy = ndarray([nx, ny])
    Rxy[:, :] = r  # NOTE  : opposite to standard BOUT convention

    Btxy = ndarray([nx, ny])
    Btxy[:, :] = Bp

    Bpxy = ndarray([nx, ny])
    Bpxy[:, :] = Bt

    Bxy = ndarray([nx, ny])
    Bxy[:, :] = sqrt(Bt ** 2 + Bp ** 2)

    hthe = ndarray([nx, ny])
    hthe[:, :] = R

    print("Writing to file '" + file + "'")

    f = DataFile()
    f.open(file, create=True)

    # Mesh size
    f.write("nx", nx)
    f.write("ny", ny)

    # Mesh spacing
    f.write("dx", dx)
    f.write("dy", dy)

    # Metric components
    f.write("Rxy", Rxy)
    f.write("Btxy", Btxy)
    f.write("Bpxy", Bpxy)
    f.write("Bxy", Bxy)
    f.write("hthe", hthe)

    # Shift
    f.write("ShiftAngle", ShiftAngle)

    # Curvature
    f.write("logB", logB)

    # Input parameters
    f.write("R", R)
    f.write("r", r)
    f.write("dr", dr)
    f.write("Bt", Bt)
    f.write("q", q)
    f.write("mxg", mxg)

    f.close()


def coordinates(
    nx,
    ny,
    nz,
    R=2.0,
    r=0.2,  # Major & minor radius
    dr=0.05,  # Radial width of domain
    Bt=1.0,  # Toroidal magnetic field
    q=5.0,  # Safety factor
    mxg=2,
):
    """
    Returns coordinates (R,Z) as a pair of arrays

    """

    h = dr / (nx - 2.0 * mxg)  # Grid spacing in r
    rminor = linspace(
        r - 0.5 * dr - (mxg - 0.5) * h, r + 0.5 * dr + (mxg - 0.5) * h, nx
    )

    print("Grid spacing: Lx = %e, Lz = %e" % (h, 2.0 * pi * r / nz))

    Rxyz = ndarray([nx, ny, nz])
    Zxyz = ndarray([nx, ny, nz])

    for y in range(0, ny):
        dtheta = y * 2.0 * pi / ny / q  # Change in poloidal angle
        theta = linspace(0, 2.0 * pi, nz, endpoint=False) + dtheta

        Rxyz[:, y, :] = R + outer(rminor, cos(theta))
        Zxyz[:, y, :] = outer(rminor, sin(theta))

    return Rxyz, Zxyz
