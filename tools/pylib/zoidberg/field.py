from builtins import object

# from math import pi, atan, cos, sin

import numpy as np

# from . import grid

class MagneticField(object):
    def __default_Byfunc(self, x,z,phi):
        return np.ones(x.shape)

    def __init__(self, grid, Bxfunc, Bzfunc, Byfunc=None, smooth=False, smooth_args={}):
        """

        Inputs
        ------
        grid = Zoidberg grid
        Bxfunc = Function for magnetic field in x
        Bzfunc = Function for magnetic field in z
        Byfunc = Function for magnetic field in y (default = 1.)
        smooth = Damp Bx, Bz to 0 at edges of grid (default off)
        smooth_args = Dict containing arguments to smooth_field_line
        """

        if Byfunc is None:
            Byfunc = self.__default_Byfunc

        self.grid = grid
        self.Bxfunc = Bxfunc
        self.Byfunc = Byfunc
        self.Bzfunc = Bzfunc

        self.smooth = smooth
        self.smooth_args = smooth_args

        if self.smooth:
            self.smooth_func = np.vectorize(self.smooth_field_line)

            self.x_width = smooth_args["x_width"] if "z_width" in smooth_args else 4
            self.z_width = smooth_args["z_width"] if "z_width" in smooth_args else 4

            self.xr_inner = smooth_args["xr_inner"] if "xr_inner" in smooth_args else self.grid.xarray[-self.x_width-1]
            self.xr_outer = smooth_args["xr_outer"] if "xr_outer" in smooth_args else self.grid.xarray[-1]
            self.xl_inner = smooth_args["xl_inner"] if "xl_inner" in smooth_args else self.grid.xarray[self.x_width]
            self.xl_outer = smooth_args["xl_outer"] if "xl_outer" in smooth_args else self.grid.xarray[0]
            self.zt_inner = smooth_args["zt_inner"] if "zt_inner" in smooth_args else self.grid.zarray[-self.z_width-1]
            self.zt_outer = smooth_args["zt_outer"] if "zt_outer" in smooth_args else self.grid.zarray[-1]
            self.zb_inner = smooth_args["zb_inner"] if "zb_inner" in smooth_args else self.grid.zarray[self.z_width]
            self.zb_outer = smooth_args["zb_outer"] if "zb_outer" in smooth_args else self.grid.zarray[0]

        self.bx = self.Bxfunc(self.grid.x_3d, self.grid.z_3d, self.grid.y_3d)
        self.by = self.Byfunc(self.grid.x_3d, self.grid.z_3d, self.grid.y_3d)
        self.bz = self.Bzfunc(self.grid.x_3d, self.grid.z_3d, self.grid.y_3d)

        if self.smooth:
            P = self.smooth_func(self.grid.x_3d, self.grid.z_3d)
            self.bx *= P
            self.bz *= P

        self.b_mag = np.sqrt(self.bx**2 + self.by**2 + self.bz**2)

    def field_direction(self, pos, phi, flatten=False):
        """Calculate the direction of the magnetic field
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
        dxdphi =  self.grid.Rmaj(x,z,phi) * self.Bxfunc(x,z,phi) / self.Byfunc(x,z,phi)
        # Rate of change of z location [m] with y angle [radians]
        dzdphi =  self.grid.Rmaj(x,z,phi) * self.Bzfunc(x,z,phi) / self.Byfunc(x,z,phi)

        if self.smooth:
            P = self.smooth_func(x, z, **self.smooth_args)
            dxdphi *= P
            dzdphi *= P

        if flatten:
            result = np.column_stack((dxdphi, dzdphi)).flatten()
        else:
            result = [dxdphi, dzdphi]

        return result

    def smooth_field_line(self, xa, za):
        """Linearly damp the field to be parallel to the edges of the box

        Should take some parameters to adjust rate of smoothing, etc.
        """

        x_left = (xa - self.xl_inner) / (self.xl_inner - self.xl_outer) + 1
        x_right = (xa - self.xr_inner) / (self.xr_inner - self.xr_outer) + 1
        z_top = (za - self.zt_inner) / (self.zt_inner - self.zt_outer) + 1
        z_bottom = (za - self.zb_inner) / (self.zb_inner - self.zb_outer) + 1

        if (xa < self.xl_inner):
            if (za < self.zb_inner):
                P = np.min([x_left, z_bottom])
            elif (za >= self.zt_inner):
                P = np.min([x_left, z_top])
            else:
                P = x_left
        elif (xa >= self.xr_inner):
            if (za < self.zb_inner):
                P = np.min([x_right, z_bottom])
            elif (za >= self.zt_inner):
                P = np.min([x_right, z_top])
            else:
                P = x_right

        elif (za < self.zb_inner):
            P = z_bottom
        elif (za > self.zt_inner):
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
        R = self.grid.Ly/float(2.*np.pi)

        # Set poloidal magnetic field
        Bpx = self.Bp + (self.grid.xarray-self.grid.Lx/2.) * self.Bpprime

        self.Bpxy = np.resize(Bpx, (self.grid.nz, self.grid.ny, self.grid.nx))
        self.Bpxy = np.transpose(self.Bpxy, (2,1,0))

        self.Bxy = np.sqrt(self.Bpxy**2 + self.Bt**2)

        def Bxfunc(x, z, phi):
            return np.zeros(x.shape)
        def Bzfunc(x, z, phi):
            return self.Bp + (x - self.grid.xcentre)*self.Bpprime

        super().__init__(grid, Bxfunc, Bzfunc)

try:
    from sympy import Symbol, Derivative, atan, atan2, cos, sin, log, pi, sqrt, lambdify
    
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
            
        def __init__(self, grid, radius=0.8, iota=1, I_coil=0.05, smooth=False, smooth_args={}):

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
            
            super().__init__(self.grid, self.Bxfunc, self.Bzfunc, smooth=smooth, smooth_args=smooth_args)

except ImportError:
    print("No Sympy module: Can't generate Stellarator fields")


class VMEC(object):
    """Read a VMEC equilibrium file
    """

    def __rolling_average(self, field):
        result = np.zeros_like(field)
        result[:,0]    = field[:,1] - 0.5*field[:,2]
        result[:,2:-2] = 0.5 * (field[:,2:-2] + field[:,3:-1])
        result[:,-1]   = 2*field[:,-2] - field[:,-3]
        return result

    def cfunct(self, field):
        """VMEC DCT
        """
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
        """VMEC DST
        """
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

    def read_vmec_file(self, vmec_file, ntheta=None, nzeta=None):
        """Read a VMEC equilibrium file
        """
        from boututils.datafile import DataFile
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

            bsupumnc = self.__rolling_average(f['bsupumnc'].T).T
            bsupvmnc = self.__rolling_average(f['bsupvmnc'].T).T

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

    def adjust_grid(self, grid):
        """Adjust grid to be consistent with the VMEC grid
        """
        grid.nx = self.nr
        grid.ny = self.nzeta
        grid.nz = self.nz

        grid.xarray = self.r_1D
        grid.yarray = self.zeta
        grid.zarray = self.z_1D

        grid.delta_x = grid.xarray[1] - grid.xarray[0]
        grid.delta_y = grid.yarray[1] - grid.yarray[0]
        grid.delta_z = grid.zarray[1] - grid.zarray[0]

        grid.Lx = grid.xarray[-1] - grid.xarray[0]
        grid.Ly = 2.*np.pi
        grid.Lz = grid.zarray[-1] - grid.zarray[0]

        grid.xcentre = grid.xarray[0] + 0.5*grid.Lx
        grid.zcentre = grid.zarray[0] + 0.5*grid.Lz

        grid.x_3d, grid.y_3d, grid.z_3d = np.meshgrid(grid.xarray, grid.yarray, grid.zarray,
                                                      indexing='ij')

        def Rmaj(x,z,phi):
            return x

        grid.Rmaj = Rmaj


    def __init__(self, grid, vmec_file, ntheta=None, nzeta=None, nr=32, nz=32):
        # Only needed here
        from scipy.interpolate import griddata, RegularGridInterpolator
        from scipy import ndimage

        self.read_vmec_file(vmec_file, ntheta, nzeta)

        self.nr = nr
        self.nz = nz

        # Make a new rectangular grid in (R,Z)
        self.r_1D = np.linspace(self.r_stz.min(), self.r_stz.max(), nr)
        self.z_1D = np.linspace(self.z_stz.min(), self.z_stz.max(), nr)
        self.R_2D, self.Z_2D = np.meshgrid(self.r_1D, self.z_1D, indexing='ij')

        # First, interpolate the magnetic field components onto (R,Z)
        self.br_rz = np.zeros( (nr, nz, self.nzeta) )
        self.bz_rz = np.zeros( (nr, nz, self.nzeta) )
        self.bphi_rz = np.zeros( (nr, nz, self.nzeta) )
        # No need to interpolate in zeta, so do this one slice at a time
        for k, (br, bz, bphi, r, z) in enumerate(zip(self.br.T, self.bz.T, self.bphi.T, self.r_stz.T, self.z_stz.T)):
            points = np.column_stack( (r.flatten(), z.flatten()) )
            self.br_rz[...,k] = griddata(points, br.flatten(), (self.R_2D, self.Z_2D),
                                         method='linear', fill_value=0.0)
            self.bz_rz[...,k] = griddata(points, bz.flatten(), (self.R_2D, self.Z_2D),
                                         method='linear', fill_value=0.0)
            self.bphi_rz[...,k] = griddata(points, bphi.flatten(), (self.R_2D, self.Z_2D),
                                           method='linear', fill_value=1.0)

        # Now we have a regular grid in (R,Z,phi) (as zeta==phi), so
        # we can get an interpolation function in 3D
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

        def Byfunc(x, z, phi):
            phi = np.mod(phi, 2.*np.pi)
            return self.bphi_interp((x, z, phi))

        self.adjust_grid(grid)

        self.magnetic_field = MagneticField(grid, Bxfunc, Bzfunc, Byfunc=Byfunc)

class FieldInterpolator(object):
    """Given Bx, By, Bz on a grid, interpolate them for zoidberg

    """
    def __init__(self, Bx, By, Bz):
        pass
