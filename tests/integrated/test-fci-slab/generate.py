from __future__ import division
from builtins import object
from past.utils import old_div
#
# Routines to generate slab meshes for FCI
#

import numpy as np
from math import pi
from scipy.integrate import odeint
import boututils.datafile as bdata
from boutdata.input import transform3D

def slab(nx, ny, nz,
         filename="fci.grid.nc",
         Lx=0.1, Ly=10., Lz = 1.,
         Bt=1.0, Bp = 0.1, Bpprime = 1.0):
    """
    nx  - Number of radial points
    ny  - Number of toroidal points (NOTE: Different to BOUT++ standard)
    nz  - Number of poloidal points

    Lx  - Radial domain size  [m]
    Ly  - Toroidal domain size [m]
    Lz  - Poloidal domain size [m]

    Bt  - Toroidal magnetic field [T]
    Bp  - Poloidal magnetic field [T]
    Bpprime - Gradient of Bp [T/m]  Bp(x) = Bp + Bpprime * x
    """

    MXG = 2

    # Make sure input types are sane
    nx = int(nx)
    ny = int(ny)
    nz = int(nz)

    Lx = float(Lx)
    Ly = float(Ly)
    Lz = float(Lz)

    delta_x = old_div(Lx,(nx-2.*MXG))
    delta_pol = old_div(Lz,(nz))
    delta_tor = old_div(Ly,(ny))

    # Coord arrays
    x = Lx * (np.arange(nx) - MXG + 0.5)/(nx - 2.*MXG)  # 0 and 1 half-way between cells
    y = np.linspace(0,Ly,ny)
    z = np.linspace(0,Lz,nz,endpoint=False)

    ############################################################

    # Effective major radius
    R = old_div(Ly, (2.*pi))

    # Set poloidal magnetic field

    Bpx = Bp + (x-old_div(Lx,2)) * Bpprime

    Bpxy = np.transpose(np.resize(Bpx, (nz, ny, nx)), (2,1,0))

    Bxy = np.sqrt(Bpxy**2 + Bt**2)[:,:,0]

    class Mappoint(object):
        def __init__(self, xt, zt):
            self.xt = xt
            self.zt = zt

            self.xt_prime = old_div(xt,delta_x) + MXG - 0.5
            self.zt_prime = old_div(zt,delta_pol)

    def unroll_map_coeff(map_list, coeff):
        coeff_array = np.transpose(np.resize(np.array([getattr(f, coeff) for f in map_list]).reshape( (nx,nz) ), (ny, nx, nz) ), (1, 0, 2) )
        return coeff_array

    def b_field(vector, y):
        x0 = old_div(Lx,2.)                    # Centre of box, where bz = 0.
        x, z = vector;
        bx = 0.
        bz = Bp + (x-x0) * Bpprime

        return [bx, bz]

    def field_line_tracer(direction, map_list):

        result = np.zeros( (nx, nz, 2) )

        for i in np.arange(0,nx):
            for k in np.arange(0,nz):
                result[i,k,:] = odeint(b_field, [x[i], z[k]], [0, delta_tor*direction])[1,:]
                map_list.append(Mappoint(result[i,k,0],result[i,k,1]))

        return result

    forward_map = []
    forward_coords = field_line_tracer(+1, forward_map)
    backward_map = []
    backward_coords = field_line_tracer(-1, backward_map)

    X,Y = np.meshgrid(x,y,indexing='ij')
    x0 = 0.5
    g_22 = old_div(((Bp + (X-x0) * Lx * Bpprime)**2 + Bt**2), Bt**2)

    with bdata.DataFile(filename, write=True, create=True) as f:
        f.write('nx', nx)
        f.write('ny', ny)
        f.write('nz', nz)
        f.write("dx", delta_x)
        f.write("dy", delta_tor)
        f.write("g_22", g_22)
        f.write("Bxy", (Bxy))

        xt_prime = unroll_map_coeff(forward_map, 'xt_prime')
        f.write('forward_xt_prime', (xt_prime))
        zt_prime = unroll_map_coeff(forward_map, 'zt_prime')
        f.write('forward_zt_prime', (zt_prime))

        xt_prime = unroll_map_coeff(backward_map, 'xt_prime')
        f.write('backward_xt_prime', (xt_prime))
        zt_prime = unroll_map_coeff(backward_map, 'zt_prime')
        f.write('backward_zt_prime', (zt_prime))


if __name__ == "__main__":
    slab(34, 64, 64, filename="fci.grid.nc")
