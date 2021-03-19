# Create an input file for testing

from netCDF4 import Dataset
import numpy
from sys import exit

# number of components
nx = 12
ny = 12
nz = 5

ivar = 1
rvar = numpy.pi
f2d = numpy.random.rand(nx, ny)
fperp = numpy.random.rand(nx, nz)
fperp2 = numpy.random.rand(nx, nz)
f3d = numpy.random.rand(nx, ny, nz)

with Dataset("test_io.grd.nc", "w") as f:
    f.createVariable("nx", numpy.int64)
    f["nx"][...] = nx
    f.createVariable("ny", numpy.int64)
    f["ny"][...] = nx
    f.createVariable("ivar", numpy.int64)
    f["ivar"][...] = ivar
    f.createVariable("rvar", numpy.float64)
    f["rvar"][...] = rvar

    f.createDimension("x", nx)
    f.createDimension("y", ny)
    f.createDimension("z", nz)
    f.createVariable("f2d", numpy.float64, ("x", "y"))
    f["f2d"][...] = f2d
    f.createVariable("fperp", numpy.float64, ("x", "z"))
    f["fperp"][...] = fperp
    f.createVariable("fperp2", numpy.float64, ("x", "z"))
    f["fperp2"][...] = fperp2
    f["fperp2"].yindex_global = 11
    f.createVariable("f3d", numpy.float64, ("x", "y", "z"))
    f["f3d"][...] = f3d

exit(0)
