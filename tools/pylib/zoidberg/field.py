from builtins import object
from past.utils import old_div

# from math import pi, atan, cos, sin

import numpy as np
from sympy import Symbol, Derivative, atan, atan2, cos, sin, log, pi, sqrt, lambdify

# from . import grid

class MagneticField(object):
    def __default_Byfunc(self, x,z,phi):
        return 1.

    def __init__(self, grid, Bxfunc, Bzfunc, Byfunc=None):

        if Byfunc is None:
            Byfunc = self.__default_Byfunc

        self.grid = grid
        self.Bxfunc = Bxfunc
        self.Byfunc = Byfunc
        self.Bzfunc = Bzfunc
        # Specifiy the magnetic field somehow
        # Two choices:
        #   1. analytic equations
        #   2. from data


    def field_direction(self, pos, phi, flatten=False):
        """
        Calculate the direction of the magnetic field
        Returns the change in x with phi and change in z with phi

        Inputs
        ------
        pos = [x,z]  with x and z in meters
        phi = toroidal angle in radians

        Returns
        -------

        (dx/dphi, dz/dphi) = ( R*Bx/Bphi, R*Bz/Bphi )
        """

        if flatten:
            position = pos.reshape((-1, 2))
            x = position[:,0]
            z = position[:,1]
        else:
            x,z = pos

        # Rate of change of x location [m] with y angle [radians]
        dxdphi =  self.grid.Rmaj * self.Bxfunc(x,z,phi) / self.Byfunc(x,z,phi)
        # Rate of change of z location [m] with y angle [radians]
        dzdphi =  self.grid.Rmaj * self.Bzfunc(x,z,phi) / self.Byfunc(x,z,phi)

        if flatten:
            result = np.column_stack((dxdphi, dzdphi)).flatten()
        else:
            result = [dxdphi, dzdphi]

        return result

    def smooth_field_line(xa,za):
        """Linearly damp the field to be parallel to the edges of the box

        Should take some parameters to adjust rate of smoothing, etc.
        """
        xr_inner = xarray[bi[-1]]
        xr_outer = xarray[-3]
        xl_inner = xarray[bi[0]]
        xl_outer = xarray[2]

        zt_inner = zarray[bk[-1]]
        zt_outer = zarray[-3]
        zb_inner = zarray[bk[0]]
        zb_outer = zarray[2]

        x_left = (xa - xl_inner) / (xl_inner - xl_outer) + 1
        x_right = (xa - xr_inner) / (xr_inner - xr_outer) + 1
        z_top = (za - zt_inner) / (zt_inner - zt_outer) + 1
        z_bottom = (za - zb_inner) / (zb_inner - zb_outer) + 1

        if (xa < xarray[bi[0]]):
            if (za <= zarray[bk[0]]):
                if (xa > za - zarray[1]):
                    P = z_bottom
                elif (xa <= za - zarray[1]):
                    P = x_left
                else:
                    P = 0
            elif (za > zarray[max(bk)]):
                if (np.abs(xa - xarray[bi[0]]) > za-zarray[max(bk)]):
                    P = x_left
                elif (np.abs(xa - xarray[bi[0]]) <= za-zarray[max(bk)]):
                    P = z_top
                else:
                    P = 0
            else:
                P = x_left
        elif (xa > xarray[max(bi)]):
            if (za <= zarray[bk[0]]):
                if (xa-xarray[max(bi)] <= np.abs(za - zarray[bk[0]])):
                    P = z_bottom
                elif (xa-xarray[max(bi)] >= np.abs(za - zarray[bk[0]])):
                    P = x_right
                else:
                    P = 0
            elif (za >= zarray[max(bk)]):
                if (xa > za +zarray[1]):
                    P = x_right
                elif (xa <= za + zarray[1]):
                    P = z_top
                else:
                    P = 0
            else:
                P = x_right

        elif (za <= zarray[bk[0]]):
            P = z_bottom
        elif (za >= zarray[max(bk)]):
            P = z_top
        else:
            P=1.
        if (P<0.):
            P=0.
        return P


class Slab(object):
    def __init__(self, grid, Bt=1.0, Bp = 0.1, Bpprime = 1.0):

        self.grid = grid
        self.Bt = Bt
        self.Bp = Bp
        self.Bpprime = Bpprime

        # Effective major radius
        R = old_div(self.grid.Ly, (2.*np.pi))

        # Set poloidal magnetic field
        Bpx = self.Bp + (self.grid.xarray-old_div(self.grid.Lx,2)) * self.Bpprime

        self.Bpxy = np.resize(Bpx, (self.grid.nz, self.grid.ny, self.grid.nx))
        self.Bpxy = np.transpose(self.Bpxy, (2,1,0))

        self.Bxy = np.sqrt(self.Bpxy**2 + self.Bt**2)

        def Bxfunc(x, z, phi):
            return np.zeros(x.shape)
        def Bzfunc(x, z, phi):
            return self.Bp + (x - self.grid.xcentre)*self.Bpprime

        self.magnetic_field = MagneticField(grid, Bxfunc, Bzfunc)


class Stellarator(object):
    def coil(self, radius, angle, iota, I):
        """Defines a single coil

        Inputs
        ------
        radius - radius to coil
        angle - initial angle of coil
        iota - rotational transform of coil
        I - current through coil

        Returns
        -------
        (x, z) - x, z coordinates of coils along phi
        """

        return (self.grid.xcentre + radius * cos(angle + iota * self.phi),
                self.grid.zcentre + radius * sin(angle + iota * self.phi), I)

    def __init__(self, grid, radius=0.8, iota=1, I_coil=0.05):

        self.x = Symbol('x')
        self.z = Symbol('z')
        self.y = Symbol('y')
        self.r = Symbol('r')
        self.r = (self.x**2 + self.z**2)**(0.5)
        self.phi = Symbol('phi')

        self.radius = radius * ((grid.Lx + grid.Lz)*0.5)
        self.grid = grid

        # Four coils equally spaced, alternating direction for current
        self.coil_list = [self.coil(self.radius, n*pi, iota, ((-1)**np.mod(i,2))*I_coil)
                          for i, n in enumerate(np.arange(4)/2.)]

        A = 0.0
        Bx = 0.0
        Bz = 0.0

        for c in self.coil_list:
            xc, zc, Ic = c
            rc = (xc**2 + zc**2)**(0.5)
            r2 = (self.x - xc)**2 + (self.z - zc)**2
            theta = atan2(self.z - zc, self.x - xc) # Angle relative to coil

            A -= Ic * 0.1 * log(r2)

            B = Ic * 0.2/sqrt(r2)

            Bx += B * sin(theta)
            Bz -= B * cos(theta)

        self.Afunc  = lambdify((self.x, self.z, self.phi), A, "numpy")
        self.Bxfunc = lambdify((self.x, self.z, self.phi), Bx, "numpy")
        self.Bzfunc = lambdify((self.x, self.z, self.phi), Bz, "numpy")

        self.magnetic_field = MagneticField(self.grid, self.Bxfunc, self.Bzfunc)
