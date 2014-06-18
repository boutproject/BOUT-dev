import numpy as np
from scipy.integrate import odeint
import boututils.datafile as bdata
from boutdata.input import transform3D

# Parameters
nx = 5
ny = 5
nz = 5
Lx = 1.
Ly = 2*np.pi
Lz = 2*np.pi
delta_x = Lx/(nx-1)
delta_y = Ly/(ny-1)
delta_z = Lz/(nz-1)

# Coord arrays
x = np.linspace(0,Lx,nx)
z = np.linspace(0,Lz,nz)

class Mappoint():
    def __init__(self, xt, zt):
        self.xt = xt
        self.zt = zt

        self.xt_prime = xt/delta_x
        self.zt_prime = zt/delta_z

        self.i_corner = int(xt/delta_x)
        self.k_corner = int(zt/delta_z)

        t_x = (xt - x[self.i_corner])/delta_x
        t_z = (zt - z[self.k_corner])/delta_z

        self.a_x = (1.-t_x)*(1.+t_x-2.*(t_x**2))
        self.a_z = (1.-t_z)*(1.+t_z-2.*(t_z**2))

        self.a_1mx = (1.-(1.-t_x))*(1.+(1.-t_x)-2.*(1.-t_x)**2)
        self.a_1mz = (1.-(1.-t_z))*(1.+(1.-t_z)-2.*(1.-t_z)**2)

        self.b_x = t_x*(1.-t_x)**2
        self.b_z = t_z*(1.-t_z)**2

        self.b_1mx = (1.-t_x)*(1.-(1.-t_x))**2
        self.b_1mz = (1.-t_z)*(1.-(1.-t_z))**2

def unroll_map_coeff(map_list, coeff):
    coeff_array = np.transpose(np.resize(np.array([getattr(f, coeff) for f in map_list]).reshape( (nx,nz) ), (ny, nx, nz) ), (1, 0, 2) )
    return coeff_array

def b_field(vector, y):
    x0 = 0.5                    # Centre of box, where bz = 0.
    x, z = vector;
    bx = 0.
    bz = (x-x0)*(-1)

    return [bx, bz]

def field_line_tracer(direction, map_list):

    result = np.zeros( (nx, nz, 2) )

    for i in np.arange(0,nx):
        for k in np.arange(0,nz):
            result[i,k,:] = odeint(b_field, [x[i], z[k]], [0, delta_y*direction])[1,:]
            result[i,k,1] = np.mod(result[i,k,1], Lz)
            map_list.append(Mappoint(result[i,k,0],result[i,k,1]))

    return result

def output_coefs(direction, map_list):

    if direction > 0:
        filename = 'forward_coefs.nc'
    else:
        filename = 'backward_coefs.nc'

    with bdata.DataFile(filename, write=True, create=True) as f:
        f.write('nx',nx)
        f.write('ny',ny)

        xt_prime = unroll_map_coeff(map_list, 'xt_prime')
        f.write('xt_prime', (xt_prime))
        zt_prime = unroll_map_coeff(map_list, 'zt_prime')
        f.write('zt_prime', (zt_prime))

        # i_corner = unroll_map_coeff(map_list, 'i_corner')
        # f.write('i_corner', transform3D(i_corner))
        # k_corner = unroll_map_coeff(map_list, 'k_corner')
        # f.write('k_corner', transform3D(k_corner))

        # a_x = unroll_map_coeff(map_list, 'a_x')
        # f.write('a_x', transform3D(a_x))
        # b_x = unroll_map_coeff(map_list, 'b_x')
        # f.write('b_x', transform3D(b_x))

        # a_z = unroll_map_coeff(map_list, 'a_z')
        # f.write('a_z', transform3D(a_z))
        # b_z = unroll_map_coeff(map_list, 'b_z')
        # f.write('b_z', transform3D(b_z))

        # a_1mx = unroll_map_coeff(map_list, 'a_1mx')
        # f.write('a_1mx', transform3D(a_1mx))
        # b_1mx = unroll_map_coeff(map_list, 'b_1mx')
        # f.write('b_1mx', transform3D(b_1mx))

        # a_1mz = unroll_map_coeff(map_list, 'a_1mz')
        # f.write('a_1mz', transform3D(a_1mz))
        # b_1mz = unroll_map_coeff(map_list, 'b_1mz')
        # f.write('b_1mz', transform3D(b_1mz))

if __name__ == "__main__":

    forward_map = []
    forward_coords = field_line_tracer(+1, forward_map)
    output_coefs(+1, forward_map)

    backward_map = []
    backward_coords = field_line_tracer(-1, backward_map)
    output_coefs(-1, backward_map)
