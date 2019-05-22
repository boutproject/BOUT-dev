# Create an input file for testing

from h5py import File
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

with File('test_io.grd.hdf5', 'w') as f:
    f['nx'] = nx
    f['ny'] = nx
    f['ivar'] = ivar
    f['rvar'] = rvar

    f['f2d'] = f2d
    f['fperp'] = fperp
    f['fperp2'] = fperp2
    f['fperp2'].attrs['yindex_global'] = 11
    f['f3d'] = f3d

exit(0)
