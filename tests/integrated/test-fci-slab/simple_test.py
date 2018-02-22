from builtins import range
from numpy import zeros, linspace, concatenate
import boututils.datafile as bdata
from boutdata.input import transform3D

# Parameters
nx = 10
ny = 20
nz = 8

shape = [nx, ny, nz]

xt_prime = zeros(shape)
zt_prime = zeros(shape)

for x in range(nx):
    # No interpolation in x
    xt_prime[x,:,:] = x

    # Each y slice scans between neighbouring z points
    for z in range(nz):
        zt_prime[x,:,z] = z + concatenate([linspace(-1, 1, ny-1), [0]])


with bdata.DataFile('simple_test.nc', write=True, create=True) as f:
    f.write('nx',nx)
    f.write('ny',ny)
    
    for direction_name in ['forward', 'backward']:
        f.write(direction_name + '_xt_prime', transform3D(xt_prime))
        f.write(direction_name + '_zt_prime', transform3D(zt_prime))
    
