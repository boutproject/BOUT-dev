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


class Slab(MagneticField):
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

        super().__init__(grid, Bxfunc, Bzfunc)


class Stellarator(MagneticField):
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

        super().__init__(self.grid, self.Bxfunc, self.Bzfunc)


class VMEC(object):
    """Read a VMEC equilibrium file
    """

    def __weird_average(self, field):
        result = np.zeros_like(field)
        result[:,0]    = field[:,1] - 0.5*field[:,2]
        result[:,2:-2] = 0.5 * (field[:,2:-2] + field[:,3:-1])
        result[:,-1]   = 2*field[:,-2] - field[:,-3]
        return result

    def cfunct(self, field):
        ns = field.shape[0]
        lt = self.theta.size
        lz = self.zeta.size
        # Create mode x angle arrays
        mt = self.xm[:,np.newaxis]*self.theta
        nz = self.xn[:,np.newaxis]*self.zeta
        # Create Trig Arrays
        cosmt = np.cos(mt)
        sinmt = np.sin(mt)
        cosnz = np.cos(nz)
        sinnz = np.sin(nz)
        # Calculate the transform
        f = np.zeros((ns,lt,lz))
        for k, field_slice in enumerate(field):
            rmn = np.repeat(field_slice[:,np.newaxis], lt, axis=1)
            a = rmn*cosmt
            b = np.dot(a.T, cosnz)
            # print("a: {}, b: {}".format(a.shape, b.shape))
            c = rmn*sinmt
            d = np.dot(c.T, sinnz)
            # print("c: {}, d: {}".format(c.shape, d.shape))
            f[k,:,:] = b - d
        return f

    def sfunct(self, field):
        ns = field.shape[0]
        lt = self.theta.size
        lz = self.zeta.size
        # Create mode x angle arrays
        mt = self.xm[:,np.newaxis]*self.theta
        nz = self.xn[:,np.newaxis]*self.zeta
        # Create Trig Arrays
        cosmt = np.cos(mt)
        sinmt = np.sin(mt)
        cosnz = np.cos(nz)
        sinnz = np.sin(nz)
        # Calculate the transform
        f = np.zeros((ns,lt,lz))
        for k, field_slice in enumerate(field):
            rmn = np.repeat(field_slice[:,np.newaxis], lt, axis=1)
            a = rmn*sinmt
            b = np.dot(a.T, cosnz)
            c = rmn*cosmt
            d = np.dot(c.T, sinnz)
            f[k,:,:] = b + d
        return f

    def __init__(self, grid, vmec_file, ntheta=None, nzeta=None, nr=32, nz=32):
        # Only needed here
        from boututils.datafile import DataFile
        from scipy.interpolate import griddata, RegularGridInterpolator
        
        # Read necessary stuff
        with DataFile(vmec_file, write=False) as f:
            self.xm = f['xm'].T
            self.xn = f['xn'].T
            ns = int(f['ns'])
            xm_big = np.repeat(self.xm[:,np.newaxis], ns, axis=1)
            xn_big = np.repeat(self.xn[:,np.newaxis], ns, axis=1)
            # s and c seem to swap meanings here...
            rumns = -f['rmnc'].T*xm_big
            rvmns = -f['rmnc'].T*xn_big
            zumnc = f['zmns'].T*xm_big
            zvmnc = f['zmns'].T*xn_big

            try:
                iasym = f['iasym']
            except KeyError:
                iasym = 0

            if iasym:
                rumnc = -f['rmns'].T*xm_big
                rvmnc = -f['rmns'].T*xn_big
                zumns = f['zmnc'].T*xm_big 
                zvmns = f['zmnc'].T*xn_big 

            bsupumnc = self.__weird_average(f['bsupumnc'].T).T
            bsupvmnc = self.__weird_average(f['bsupvmnc'].T).T

            if ntheta is None:
                self.ntheta = int(f['mpol'])
            else:
                self.ntheta = ntheta

            if nzeta is None:
                self.nzeta = int(f['ntor']) + 1
            else:
                self.nzeta = nzeta

            self.theta = np.linspace(0, 2*np.pi, self.ntheta)
            self.zeta  = np.linspace(0, 2*np.pi, self.nzeta)
            # R, Z on (s, theta, zeta)
            self.r_stz = self.cfunct(f['rmnc'])
            self.z_stz = self.sfunct(f['zmns'])

            bu = self.cfunct(bsupumnc)
            bv = self.cfunct(bsupvmnc)

            drdu = self.sfunct(rumns.T)
            drdv = self.sfunct(rvmns.T)
            dzdu = self.cfunct(zumnc.T)
            dzdv = self.cfunct(zvmnc.T)
            if iasym:
                self.r_stz = self.r_stz + self.sfunct(f['rmnc'])
                drdu = drdu + self.cfunct(rumnc.T)
                drdv = drdv + self.cfunct(rvmnc.T)
                self.z_stz = self.z_stz + self.cfunct(f['zmnc'])
                dzdu = dzdu + self.sfunct(zumns.T)
                dzdv = dzdv + self.sfunct(zvmns.T)

        # Convert to bR, bZ, bphi
        self.br = bu*drdu + bv*drdv
        self.bphi = self.r_stz*bv
        self.bz = bu*dzdu + bv*dzdv
        # self.br_pol = bu*drdu;
        # self.bz_pol = bu*dzdu;

        # Make a new rectangular grid in (R,Z)
        self.r_1D = np.linspace(self.r_stz.min(), self.r_stz.max(), nr)
        self.z_1D = np.linspace(self.z_stz.min(), self.z_stz.max(), nr)
        self.R_2D, self.Z_2D = np.meshgrid(self.r_1D, self.z_1D, indexing='ij')

        # Interpolate 
        self.br_rz = np.zeros( (nr, nz, self.nzeta) )
        self.bz_rz = np.zeros( (nr, nz, self.nzeta) )
        self.bphi_rz = np.zeros( (nr, nz, self.nzeta) )

        for k, (br, bz, bphi, r, z) in enumerate(zip(self.br.T, self.bz.T, self.bphi.T, self.r_stz.T, self.z_stz.T)):
            points = np.column_stack( (r.flatten(), z.flatten()) )
            self.br_rz[...,k] = griddata(points, br.flatten(), (self.R_2D, self.Z_2D),
                                         method='linear', fill_value=0.0)
            self.bz_rz[...,k] = griddata(points, bz.flatten(), (self.R_2D, self.Z_2D),
                                         method='linear', fill_value=0.0)
            self.bphi_rz[...,k] = griddata(points, bphi.flatten(), (self.R_2D, self.Z_2D),
                                           method='linear', fill_value=1.0)

        # zeta_grid = np.repeat(self.zeta[np.newaxis, :], self.R[1], axis=0)
        # zeta_grid = np.repeat(zeta_grid[np.newaxis, :, :], ns, axis=0)
        points = ( self.r_1D, self.z_1D, self.zeta )
        self.br_interp = RegularGridInterpolator(points, self.br_rz, bounds_error=False, fill_value=0.0)
        self.bz_interp = RegularGridInterpolator(points, self.bz_rz, bounds_error=False, fill_value=0.0)
        self.bphi_interp = RegularGridInterpolator(points, self.bphi_rz, bounds_error=False, fill_value=1.0)

        def Bxfunc(x, z, phi):
            phi = np.mod(phi, 2.*np.pi)
            return self.br_interp((x, z, phi))

        def Bzfunc(x, z, phi):
            phi = np.mod(phi, 2.*np.pi)
            return self.bz_interp((x, z, phi))

        def bphifunc(x, z, phi):
            phi = np.mod(phi, 2.*np.pi)
            return self.bphi_interp((x, z, phi))

        self.magnetic_field = MagneticField(grid, Bxfunc, Bzfunc)
======= end
