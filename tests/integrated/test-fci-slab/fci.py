from __future__ import division
from builtins import object
from past.utils import old_div
import numpy as np
from math import pi
from scipy.integrate import odeint
import boututils.datafile as bdata
from boutdata.input import transform3D

# Parameters
nx = 34
########## y is toroidal!
ny = 64
########## z is poloidal!
nz = 64

Lx = 0.1     # Radial domain size [m]
Ltor = 10.   # "Toroidal" length [m]
Lpol = 1.    # "Poloidal" length [m]

delta_x = old_div(Lx,(nx))
delta_pol = old_div(Lpol,(nz))
delta_tor = old_div(Ltor,(ny))

Bt  = 1.0   # Magnetic field [T]
Bp  = 0.1   # Poloidal field at the middle of the domain [T]
Bpprime = 1.0  # Bp gradient [T/m] Bp(x) = Bp + Bpprime * x 

# Coord arrays
x = np.linspace(0,Lx,nx)
y = np.linspace(0,Ltor,ny)
z = np.linspace(0,Lpol,nz,endpoint=False)

############################################################

# Effective major radius
R = old_div(Ltor, (2.*pi))

# Set poloidal magnetic field

Bpx = Bp + (x-old_div(Lx,2)) * Bpprime

Bpxy = np.transpose(np.resize(Bpx, (nz, ny, nx)), (2,1,0))

Bxy = np.sqrt(Bpxy**2 + Bt**2)

############################################################

class Mappoint(object):
    def __init__(self, xt, zt):
        self.xt = xt
        self.zt = zt

        self.xt_prime = old_div(xt,delta_x)
        self.zt_prime = old_div(zt,delta_pol)

def unroll_map_coeff(map_list, coeff):
    coeff_array = np.transpose(np.resize(np.array([getattr(f, coeff) for f in map_list]).reshape( (nx,nz) ), (ny, nx, nz) ), (1, 0, 2) )
    return coeff_array

def b_field(vector, y):
    x0 = 0.05                    # Centre of box, where bz = 0.
    x, z = vector;
    bx = 0.
    bz = Bp + (x-x0) * Bpprime

    return [bx, bz]

def field_line_tracer(direction, map_list):

    result = np.zeros( (nx, nz, 2) )

    for i in np.arange(0,nx):
        for k in np.arange(0,nz):
            result[i,k,:] = odeint(b_field, [x[i], z[k]], [0, delta_tor*direction])[1,:]
            result[i,k,1] = np.mod(result[i,k,1], Lpol)

            map_list.append(Mappoint(result[i,k,0],result[i,k,1]))

    return result

if __name__ == "__main__":

    forward_map = []
    forward_coords = field_line_tracer(+1, forward_map)
    backward_map = []
    backward_coords = field_line_tracer(-1, backward_map)

    X,Y = np.meshgrid(x,y,indexing='ij')
    x0 = 0.5
    g_22 = np.sqrt(((Bp + (X-x0) * Lx * Bpprime)**2 + 1))

    with bdata.DataFile('fci.grid.nc', write=True, create=True) as f:
        f.write('nx', nx)
        f.write('ny', ny)
        f.write('nz', nz)
        f.write("dx", delta_x)
        f.write("dy", delta_tor)
        f.write("g_22", g_22)
        f.write("Bxy", transform3D(Bxy))
    
        xt_prime = unroll_map_coeff(forward_map, 'xt_prime')
        f.write('forward_xt_prime', transform3D(xt_prime))
        zt_prime = unroll_map_coeff(forward_map, 'zt_prime')
        f.write('forward_zt_prime', transform3D(zt_prime))

        xt_prime = unroll_map_coeff(backward_map, 'xt_prime')
        f.write('backward_xt_prime', transform3D(xt_prime))
        zt_prime = unroll_map_coeff(backward_map, 'zt_prime')
        f.write('backward_zt_prime', transform3D(zt_prime))
