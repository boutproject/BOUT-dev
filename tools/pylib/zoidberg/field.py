from builtins import object

import numpy as np
import sys


try:
    from . import boundary
except ImportError:
    import boundary


if sys.version_info >= (3, 0):
    pickle_read_mode = "rb"
    pickle_write_mode = "wb"
else:
    pickle_read_mode = "r"
    pickle_write_mode = "w"


class MagneticField(object):
    """Represents a magnetic field in either Cartesian or cylindrical
    geometry

    This is the base class, you probably don't want to instantiate one
    of these directly. Instead, create an instance of one of the
    subclasses.

    Functions which can be overridden

    - Bxfunc = Function for magnetic field in x
    - Bzfunc = Function for magnetic field in z
    - Byfunc = Function for magnetic field in y (default = 1.)
    - Rfunc = Function for major radius. If None, y is in meters

    Attributes
    ----------
    boundary
        An object with an "outside" function. See :py:obj:`zoidberg.boundary`

    attributes : A dictionary of string -> function(x,z,phi)
         Contains attributes to be written to the output

    See Also
    --------
    Slab : A straight field in normal Cartesian coordinates
    CurvedSlab : A field in curvilinear coordinates
    StraightStellarator : A rotating ellipse stellarator without curvature
    RotatingEllipse : A rotating ellipse stellarator with curvature
    VMEC : A numerical field from a VMEC equilibrium file
    GEQDSK : A numerical field from an EFIT g-file

    """

    boundary = boundary.NoBoundary() # An optional Boundary object
    attributes = {}

    def Bxfunc(self, x, z, phi):
        """Magnetic field in x direction at given coordinates

        Parameters
        ----------
        x, z, phi : array_like
            X, Z, and toroidal coordinates

        Returns
        -------
        ndarray
            X-component of the magnetic field

        """
        return np.zeros(x.shape)

    def Byfunc(self, x, z, phi):
        """Magnetic field in y direction at given coordinates

        Parameters
        ----------
        x, z, phi : array_like
            X, Z, and toroidal coordinates

        Returns
        -------
        ndarray
            Y-component of the magnetic field

        """
        return np.ones(x.shape)

    def Bzfunc(self, x, z, phi):
        """Magnetic field in z direction at given coordinates

        Parameters
        ----------
        x, z, phi : array_like
            X, Z, and toroidal coordinates

        Returns
        -------
        ndarray
            Z-component of the magnetic field

        """
        return np.zeros(x.shape)

    def Rfunc(self, x, z, phi):
        """Major radius [meters]

        Returns None if in Cartesian coordinates

        Parameters
        ----------
        x, z, phi : array_like
            X, Z, and toroidal coordinates

        Returns
        -------
        ndarray
            The major radius

        """
        return None

    def pressure(self, x, z, phi):
        """Pressure [Pascals]

        Parameters
        ----------
        x, z, phi : array_like
            X, Z, and toroidal coordinates

        Returns
        -------
        ndarray
            The plasma pressure

        """
        return 0.0

    def Bmag(self, x, z, phi):
        """Magnitude of the magnetic field

        .. math ::

           Bmag = \sqrt(B_x^2 + B_y^2 + B_z^2)

        Parameters
        ----------
        x, z, phi : array_like
            X, Z, and toroidal coordinates

        Returns
        -------
        ndarray
            The magnitude of the magnetic field

        """
        return np.sqrt( self.Bxfunc(x,z,phi)**2 + self.Byfunc(x,z,phi)**2 + self.Bzfunc(x,z,phi)**2 )

    def field_direction(self, pos, ycoord, flatten=False):
        """Calculate the direction of the magnetic field
        Returns the change in x with phi and change in z with phi

        Parameters
        ----------
        pos : ndarray
            2-D NumPy array, with the second dimension being [x,z],
            with x and z in meters
        ycoord : float
            Toroidal angle in radians if cylindrical coordinates,
            metres if Cartesian
        flatten : bool, optional
            If True, return a flattened form of the vector
            components. This is useful for passing to
            :py:obj:`~zoidberg.fieldtracer.FieldTracer`

        Returns
        -------
        (dx/dy, dz/dy) : list of floats or ndarray
            - ``= (R*Bx/Bphi, R*Bz/Bphi)`` if cylindrical
            - ``= (Bx/By, Bz/By)`` if Cartesian

        """

        # Input array must have an even number of points
        assert len(pos) % 2 == 0

        if flatten:
            position = pos.reshape((-1, 2))
            x = position[:,0]
            z = position[:,1]
        else:
            x,z = pos

        By = self.Byfunc(x,z,ycoord)
        Rmaj = self.Rfunc(x,z,ycoord) # Major radius. None if Cartesian

        if Rmaj is not None:
            # In cylindrical coordinates

            if np.amin(np.abs(By)) < 1e-8:
                # Very small By
                print(x,z,ycoord, By)
                raise ValueError("Small By ({}) at (x={}, y={}, z={})".format(By, x, ycoord, z))

            R_By = Rmaj / By
            # Rate of change of x location [m] with y angle [radians]
            dxdphi =  R_By * self.Bxfunc(x,z,ycoord)
            # Rate of change of z location [m] with y angle [radians]
            dzdphi =  R_By * self.Bzfunc(x,z,ycoord)
        else:
            # In Cartesian coordinates

            # Rate of change of x location [m] with y angle [radians]
            dxdphi =  self.Bxfunc(x,z,ycoord) / By
            # Rate of change of z location [m] with y angle [radians]
            dzdphi =  self.Bzfunc(x,z,ycoord) / By

        if flatten:
            result = np.column_stack((dxdphi, dzdphi)).flatten()
        else:
            result = [dxdphi, dzdphi]

        return result


class Slab(MagneticField):
    """Represents a magnetic field in an infinite flat slab

    Magnetic field in ``z = Bz + (x - xcentre) * Bzprime``

    Coordinates (x,y,z) assumed to be Cartesian, all in metres

    Parameters
    ----------
    By : float, optional
        Magnetic field in y direction
    Bz : float, optional
        Magnetic field in z at xcentre
    xcentre : float, optional
        Reference x coordinate
    Bzprime : float, optional
        Rate of change of Bz with x

    """
    def __init__(self, By=1.0, Bz=0.1, xcentre=0.0, Bzprime=1.0):

        By = float(By)
        Bz = float(Bz)
        xcentre = float(xcentre)
        Bzprime = float(Bzprime)

        self.By = By
        self.Bz = Bz
        self.xcentre = xcentre
        self.Bzprime = Bzprime

    def Bxfunc(self, x, z, phi):
        return np.zeros(x.shape)

    def Byfunc(self, x, z, phi):
        return np.full(x.shape, self.By)

    def Bzfunc(self, x, z, phi):
        return self.Bz + (x - self.xcentre)*self.Bzprime


class CurvedSlab(MagneticField):
    """
    Represents a magnetic field in a curved slab geometry

    Magnetic field in ``z = Bz + (x - xcentre) * Bzprime``

    x  - Distance in radial direction [m]
    y  - Azimuthal (toroidal) angle
    z  - Height [m]

    Parameters
    ----------
    By : float
        Magnetic field in y direction
    Bz : float
        Magnetic field in z at xcentre (float)
    xcentre : float
        Reference x coordinate
    Bzprime : float
        Rate of change of Bz with x
    Rmaj : float
        Major radius of the slab

    """
    def __init__(self, By=1.0, Bz=0.1, xcentre=0.0, Bzprime=1.0, Rmaj=1.0):

        By = float(By)
        Bz = float(Bz)
        xcentre = float(xcentre)
        Bzprime = float(Bzprime)
        Rmaj = float(Rmaj)

        self.By = By
        self.Bz = Bz
        self.xcentre = xcentre
        self.Bzprime = Bzprime
        self.Rmaj = Rmaj

        # Set poloidal magnetic field
        #Bpx = self.Bp + (self.grid.xarray-self.grid.Lx/2.) * self.Bpprime
        #self.Bpxy = np.resize(Bpx, (self.grid.nz, self.grid.ny, self.grid.nx))
        #self.Bpxy = np.transpose(self.Bpxy, (2,1,0))
        #self.Bxy = np.sqrt(self.Bpxy**2 + self.Bt**2)

    def Bxfunc(self, x, z, phi):
        return np.zeros(x.shape)

    def Byfunc(self, x, z, phi):
        return np.full(x.shape, self.By)

    def Bzfunc(self, x, z, phi):
        return self.Bz + (x - self.xcentre)*self.Bzprime

    def Rfunc(self, x, z, phi):
        return np.full(x.shape, self.Rmaj)

try:
    from sympy import Symbol, atan2, cos, sin, log, pi, sqrt, lambdify, Piecewise, Sum, gamma, And, factorial, diff

    class StraightStellarator(MagneticField):
        """A "rotating ellipse" stellarator without curvature

        Parameters
        ----------
        xcentre : float, optional
            Middle of the domain in x [m]
        zcentre : float, optional
            Middle of the domain in z [m]
        radius : float, optional
            Radius of coils [meters]
        yperiod : float, optional
            The period over which the coils return to their original position
        I_coil : float, optional
            Current in each coil

        """

        def coil(self, xcentre, zcentre, radius, angle, iota, I):
            """Defines a single coil

            Parameters
            ----------
            radius : float
                Radius to coil
            angle : float
                Initial angle of coil
            iota : float
                Rotational transform of coil
            I : float
                Current through coil

            Returns
            -------
            (x, z) - x, z coordinates of coils along phi
            """

            return (xcentre + radius * cos(angle + iota * self.phi),
                    zcentre + radius * sin(angle + iota * self.phi), I)

        def __init__(self, xcentre=0.0, zcentre=0.0, radius=0.8, yperiod=np.pi,
                     I_coil=0.05, smooth=False, smooth_args={}):
            xcentre = float(xcentre)
            zcentre = float(zcentre)
            radius = float(radius)
            yperiod = float(yperiod)

            iota = 2.*np.pi / yperiod

            self.x = Symbol('x')
            self.z = Symbol('z')
            self.y = Symbol('y')
            self.r = Symbol('r')
            self.r = (self.x**2 + self.z**2)**(0.5)
            self.phi = Symbol('phi')

            self.xcentre = xcentre
            self.zcentre = zcentre
            self.radius = radius

            # Four coils equally spaced, alternating direction for current
            self.coil_list = [self.coil(xcentre, zcentre, radius, n*pi, iota, ((-1)**np.mod(i,2))*I_coil)
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

    class RotatingEllipse(MagneticField):
        """A "rotating ellipse" stellarator
        Parameters
        ----------
        xcentre : float, optional
            Middle of the domain in x [m]
        zcentre : float, optional
            Middle of the domain in z [m]
        radius : float, optional
            Radius of coils [meters]
        yperiod : float, optional
            The period over which the coils return to their original position
        I_coil : float, optional
            Current in each coil
        Btor : float, optional
            Toroidal magnetic field strength
        """

        def coil(self, xcentre, zcentre, radius, angle, iota, I):
            """Defines a single coil
            Parameters
            ----------
            radius : float
                Radius to coil
            angle : float
                Initial angle of coil
            iota : float
                Rotational transform of coil
            I : float
                Current through coil
            Returns
            -------
            (x, z) - x, z coordinates of coils along phi
            """

            return (xcentre + radius * cos(angle + iota * self.phi),
                    zcentre + radius * sin(angle + iota * self.phi), I)

        def __init__(self, xcentre=0.0, zcentre=0.0, radius=0.8, yperiod=np.pi,
                     I_coil=0.05, Btor=1.0, smooth=False, smooth_args={}):
            xcentre = float(xcentre)
            zcentre = float(zcentre)
            radius = float(radius)
            yperiod = float(yperiod)
            Btor = float(Btor)

            iota = 2.*np.pi / yperiod

            self.x = Symbol('x')
            self.z = Symbol('z')
            self.y = Symbol('y')
            self.r = Symbol('r')
            self.r = (self.x**2 + self.z**2)**(0.5)
            self.phi = Symbol('phi')

            self.xcentre = xcentre
            self.zcentre = zcentre
            self.radius = radius

            # Four coils equally spaced, alternating direction for current
            self.coil_list = [self.coil(xcentre, zcentre, radius, n*pi, iota, ((-1)**np.mod(i,2))*I_coil)
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

            By = Btor/self.x
            self.Afunc  = lambdify((self.x, self.z, self.phi), A, "numpy")

            self.Bxfunc = lambdify((self.x, self.z, self.phi), Bx, "numpy")
            self.Bzfunc = lambdify((self.x, self.z, self.phi), Bz, "numpy")
            self.Byfunc = lambdify((self.x, self.z, self.phi), By, "numpy")

        def Rfunc(self, x, z, phi):
            return np.full(x.shape, x)
        
    class DommaschkPotentials(MagneticField):
        """A magnetic field generator using the Dommaschk potentials. 
        Parameters
        ----------
        A: Coefficient matrix for the torodial and polidial harmonics. Form: (m,l,(a,b,c,d))
        R_0: major radius [m]
        B_0: magnetic field on axis [T]

        Important Methods
        -----------------
        Bxfunc/Byfunc/Bzfunc(x,z,y): Returns magnetic field in radial/torodial/z-direction
        Sfunc(x,z,y): Returns approximate magnetic surface invariant for Dommaschk potentials. Use this to visualize flux surfaces
        


        """

        
        def __init__(self, A, R_0=1.0, B_0=1.0):

            self.R_0 = R_0
            self.B_0 = B_0
            

            self.R = Symbol('R')
            self.phi = Symbol('phi')
            self.Z = Symbol('Z')


            self.m = Symbol('m')
            self.l = Symbol('l')
            self.n = Symbol('n')
            self.k = Symbol('k')

            self.A = A

            
            self.P = self.U(self.A).doit().subs([(self.R, self.R/self.R_0), (self.Z, self.Z/R_0)])
            self.P_hat = self.U_hat(self.A).doit().subs([(self.R, self.R/self.R_0), (self.Z, self.Z/R_0)])

            S = (0.5*(log(self.R/self.R_0)**2 + (self.Z/R_0)**2) - self.R/self.R_0 * (log(self.R/self.R_0)*self.R_0*diff(self.P_hat, self.R) + self.Z*self.R/self.R_0*diff(self.P_hat, self.Z)))#.subs([(self.R, self.R/self.R_0), (self.Z, self.Z)])


            Bx = R_0*diff(self.P, self.R)
            By = R_0/self.R * diff(self.P, self.phi)
            Bz = diff(self.P, self.Z)

            self.Sf = lambdify((self.R, self.phi, self.Z), S, "numpy")

            self.Bxf = lambdify((self.R, self.phi, self.Z), Bx, "numpy")
            self.Byf = lambdify((self.R, self.phi, self.Z), By, "numpy")
            self.Bzf = lambdify((self.R, self.phi, self.Z), Bz, "numpy")

        def Bxfunc(self, x, z, phi):
            
            return self.Bxf(x, phi, z)/self.Byf(self.R_0, 0, 0)*self.B_0
        
        def Byfunc(self, x, z, phi):


            return self.Byf(x, phi, z)/self.Byf(self.R_0, 0, 0)*self.B_0

        def Bzfunc(self, x, z, phi):

            return self.Bzf(x, phi, z)/self.Byf(self.R_0, 0, 0)*self.B_0

        def Sfunc(self, x, z, y):
            """
            Parameters
            ----------
            x: radial coordinates normalized to R_0
            z: binormal coordinate
            y: torodial angle normalized to 2*pi

            Returns
            -------
            Approximate magnetic surface invariant S at location (x,z,y). This is from the original Dommaschk paper. Use to visualize flux surfaces
            """
            return self.Sf(x,y,z)

        def Rfunc(self, x, z, phi):
            """
            Parameters
            ----------
            x: radial coordinates normalized to R_0
            z: binormal coordinate
            y: torodial angle normalized to 2*pi

            Returns
            -------
            Radial coordinate x
            """

            return x



        def CD(self, m, k):
            """
            Parameters
            ----------
            m: torodial harmonic
            k: summation index in D

            Returns:
            --------
            Sympy function CD_mk (R) (Dirichlet boudary conditions)
            """

            n = Symbol('n')
            b = Symbol('b')
            i = Symbol('i')

            alpha = Piecewise(((((-1)**n)/(gamma(b + n + 1) * gamma(n + 1) * 2**(2 * n + b) )), n >= 0), (0, True))
            alpha_st = Piecewise((alpha * (2*n + b), n >= 0), (0, True))


            beta = Piecewise(((gamma(b  - n))/(gamma(n + 1) * 2**(2*n - b + 1)), And(n >= 0, n < b)), (0, True))
            beta_st = Piecewise((beta * (2*n - b), And(n >= 0, n < b)), (0, True))


            delta = Piecewise((alpha/2 * Sum(1/i + 1/(b + i), (i, 1, n+1)), n > 0), (0, True))
            delta_st = Piecewise((delta * (2 * n + b), n > 0), (0, True))

            j = Symbol('j')

            CD = Sum(-( alpha.subs([(n, j), (b,m)]) * ( alpha_st.subs([(n, k - m - j), (b, m)]) * log(self.R) + delta_st.subs([(n, k - m - j), (b, m)]) - alpha.subs([(n, k - m - j), (b, m)]) ) - delta.subs([(n, j), (b,m)]) * alpha_st.subs([(n, k - m - j), (b, m)]) + alpha.subs([(n, j), (b,m)]) * beta_st.subs([(n, k - j), (b, m)])) * self.R**(2*j + m)\
                + beta.subs([(n, j), (b, m)]) * alpha_st.subs([(n, k - j), (b, m)]) * self.R**(2*j - m), (j, 0, k))

            return CD


        def CN(self, m, k):
            """
            Parameters
            ----------
            m: torodial harmonic
            k: summation index in N

            Returns:
            --------
            Sympy function CN_mk (R) (Neumann boundary conditions)
            """


            n = Symbol('n')
            b = Symbol('b')
            i = Symbol('i')


            alpha = Piecewise(((((-1)**n)/(gamma(b + n + 1) * gamma(n + 1) * 2**(2 * n + b) )), n >= 0), (0, True))
            alpha_st = Piecewise((alpha * (2*n + b), n >= 0), (0, True))


            beta = Piecewise(((gamma(b  - n))/(gamma(n + 1) * 2**(2*n - b + 1)), And(n >= 0, n < b)), (0, True))
            beta_st = Piecewise((beta * (2*n - b), And(n >= 0, n < b)), (0, True))


            delta = Piecewise((alpha/2 * Sum(1/i + 1/(b + i), (i, 1, n+1)), n > 0), (0, True))
            delta_st = Piecewise((delta * (2 * n + b), n > 0), (0, True))

            j = Symbol('j')
            
            CN = Sum(( alpha.subs([(n, j), (b,m)]) * ( alpha.subs( [(n, k - m - j), (b, m)] ) * log(self.R) + delta.subs([(n, k - m - j), (b, m)])) - delta.subs([(n, j), (b,m)]) * alpha.subs( [(n, k - m - j), (b, m)] ) + alpha.subs([(n, j), (b,m)]) * beta.subs([(n, k - j), (b, m)]) ) * self.R**( 2*j + m )\
                 - beta.subs([(n, j), (b,m)]) * alpha.subs([(n, k-j),(b, m)]) * self.R**(2*j - m), (j, 0, k))

            

            return CN


        def D(self, m, n):
            """
            Parameters
            ----------
            m: torodial mode number
            n: summation index in  V

            Returns:
            --------
            Sympy function D_mn (R, Z) (Dirichlet boundary conditions)
            """
            
            i = Symbol('i')
            D = log(1)
            k_arr = np.arange(0,int(n/2)+1,1)
            CD_f = self.CD(m, i)
            

            for k in k_arr:
                D += ((self.Z**(n - 2*k))/factorial(n-2*k)*CD_f.subs(i,k))
                

            return D



        def N(self, m, n):
            """
            Parameters
            ----------
            m: torodial mode number
            n: summation index in V

            Returns:
            --------
            Sympy function N_mn (R, Z) (Neumann boundary conditions)
            """

            i = Symbol('i')
            N = log(1)
            k_arr = np.arange(0, int(n/2)+1, 1)
            CN_f = self.CN(m, i)
            for k in k_arr:
                N += ((self.Z**(n - 2*k))/factorial(n-2*k)*CN_f.subs(i,k))
        
            return N


        def V(self, m, l, a, b, c, d):
            """
            Parameters
            ----------
            m: torodial mode number
            l: polodial mode number
            a,b,c,d: Coefficients for m,l-th Dommaschk potential (elements of matrix A)

            Returns:
            --------
            Sympy function V_ml
            """

            V = ( a*cos(m*self.phi) + b*sin(m*self.phi) ) * self.D(m, l) + ( c*cos(m*self.phi) + d*sin(m*self.phi) ) * self.N(m, l-1) 

            return V



        def U(self, A):
            """
            Parameters
            ----------
            A: Coefficient matrix for the torodial and polidial harmonics. Form: (m,l,(a,b,c,d))
            
            Returns
            -----------------
            U: Superposition of all modes given in A
            
            """
            U = self.phi
            for i in range(A.shape[0]):
                for j  in range(A.shape[1]):
                    if A[i,j,0] or A[i,j,1] or A[i,j,2] or A[i,j,3] != 0:

                        U += self.V(i, j, A[i,j,0], A[i,j,1], A[i,j,2], A[i,j,3])

            return U

        def V_hat(self, m, l, a, b, c, d):
            """
            Parameters
            ----------
            m: torodial mode number
            l: polodial mode number
            a,b,c,d: Coefficients for m,l-th Dommaschk potential (elements of matrix A)

            Returns:
            --------
            Sympy function V_hat_ml; Similar to V; needed for calculation of magnetic surface invariant S
            """

            V = ( a*cos(m*(self.phi-np.pi/(2*m))) + b*sin(m*(self.phi-np.pi/(2*m)) ) ) * self.D(m, l) + ( c*cos(m*(self.phi-np.pi/(2*m))) + d*sin(m*(self.phi-np.pi/(2*m)) ) ) * self.N(m, l-1) 
            

            return V

        def U_hat(self, A):
            """
            Parameters
            ----------
            A: Coefficient matrix for the torodial and polidial harmonics. Form: (m,l,(a,b,c,d))

            Returns
            -----------------
            U: Superposition of all modes given in A

            """
            
            U = log(1)
            for i in range(A.shape[0]):
                for j  in range(A.shape[1]):
                    if A[i,j,0] or A[i,j,1] or A[i,j,2] or A[i,j,3] != 0:
                        U += self.V_hat(i, j, A[i,j,0], A[i,j,1], A[i,j,2], A[i,j,3]) * Piecewise((self.phi, i==0), (1/i, i>0))

            return U




    class Screwpinch(MagneticField):
        
        def __init__(self, xcentre=1.5, zcentre=0.0, shear = 0, yperiod=2*np.pi, Btor=1.0):
            self.x = Symbol('x')
            self.z = Symbol('z')
            self.y = Symbol('y')
            self.r = Symbol('r')
            self.r = ((self.x-xcentre)**2 + (self.z-zcentre)**2)**(0.5)
            
            self.phi = Symbol('phi')
            
            alpha = shear
            self.theta = atan2(self.z - zcentre, self.x - xcentre)
            A = alpha*self.r**2
            Bx = -alpha*self.r*self.r*sin(self.theta)
            Bz = alpha*self.r*self.r*cos(self.theta)
            By = Btor/self.x
            
            self.Afunc  = lambdify((self.x, self.z, self.phi), A, "numpy")
            self.Bxfunc = lambdify((self.x, self.z, self.phi), Bx, "numpy")
            self.Bzfunc = lambdify((self.x, self.z, self.phi), Bz, "numpy")
            self.Byfunc = lambdify((self.x, self.z, self.phi), By, "numpy")
        def Rfunc(self, x, z, phi):
            return np.full(x.shape, x)
               
except ImportError:
    class StraightStellarator(MagneticField):
        """
        Invalid StraightStellarator, since no Sympy module.

        Rather than printing an error on startup, which may
        be missed or ignored, this raises
        an exception if StraightStellarator is ever used.
        """
        def __init__(self, *args, **kwargs):
            raise ImportError("No Sympy module: Can't generate StraightStellarator fields")

    class RotatingEllipse(MagneticField):
        """
        Invalid RotatingEllipse, since no Sympy module.
        Rather than printing an error on startup, which may
        be missed or ignored, this raises
        an exception if StraightStellarator is ever used.
        """
        def __init__(self, *args, **kwargs):
            raise ImportError("No Sympy module: Can't generate RotatingEllipse fields")

    class Screwpinch(MagneticField):
        """
        Invalid screwpinch, since no Sympy module.
        Rather than printing an error on startup, which may
        be missed or ignored, this raises
        an exception if StraightStellarator is ever used.
        """
        def __init__(self, *args, **kwargs):
            raise ImportError("No Sympy module: Can't generate screwpinch fields")
        
class VMEC(MagneticField):
    """A numerical magnetic field from a VMEC equilibrium file

    Parameters
    ----------
    vmec_file : str
        Name of the VMEC file to read
    ntheta : int, optional
        Number of theta points to use (default: use 'mpol' from VMEC file)n
    nzeta : int, optional
        Number of zeta points to use (default: use 'ntor' from VMEC file)
    nr : int
        Number of R points to use
    nz : int
        Number of Z points to use

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

    def __init__(self, vmec_file, ntheta=None, nzeta=None, nr=32, nz=32):
        # Only needed here
        from scipy.interpolate import griddata, RegularGridInterpolator

        self.read_vmec_file(vmec_file, ntheta, nzeta)

        self.nr = nr
        self.nz = nz

        # Make a new rectangular grid in (R,Z)
        self.r_1D = np.linspace(self.r_stz.min(), self.r_stz.max(), nr)
        self.z_1D = np.linspace(self.z_stz.min(), self.z_stz.max(), nz)
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

    def Bxfunc(self, x, z, phi):
        phi = np.mod(phi, 2.*np.pi)
        return self.br_interp((x, z, phi))

    def Bzfunc(self, x, z, phi):
        phi = np.mod(phi, 2.*np.pi)
        return self.bz_interp((x, z, phi))

    def Byfunc(self, x, z, phi):
        phi = np.mod(phi, 2.*np.pi)
        return self.bphi_interp((x, z, phi))

    def Rfunc(self, x, z, phi):
        """
        Major radius
        """
        return x


class SmoothedMagneticField(MagneticField):
    """Represents a magnetic field which is smoothed so it never leaves
    the boundaries of a given grid.

    Parameters
    ----------
    field : :py:obj:`zoidberg.field.MagneticField`
        A MagneticField object
    grid : :py:obj:`zoidberg.grid.Grid`
        A Grid object
    xboundary : int, optional
        Number of grid points in x over which the magnetic field is smoothed
    zboundary : int, optional
        Number of grid points in x over which the magnetic field is smoothed

    """
    def __init__(self, field, grid, xboundary=None, zboundary=None):
        """
        """

        self.field = field
        self.grid = grid
        self.xboundary = xboundary
        self.zboundary = zboundary

    def Bxfunc(self, x, z, phi):
        # Get closest y (phi) grid index
        ind = np.argmin(np.abs( phi - self.grid.ycoords))
        if phi < self.grid.ycoords[ind]:
            ind -= 1
        # phi now between ind and ind+1
        grid_d, y_d = self.grid.getPoloidalGrid(ind)
        grid_u, y_u = self.grid.getPoloidalGrid(ind+1)

        # Get x,z indices from poloidal grids
        x_d, z_d = grid_d.findIndex(x, z)
        x_u, z_u = grid_d.findIndex(x, z)

    def Byfunc(self, x, z, phi):
        """
        Not modified by smoothing
        """
        return self.field.Byfunc(x, z, phi)

    def Bxfunc(self, x, z, phi):
        pass

    def Rfunc(self, x, z, phi):
        return self.field.Rfunc(x, z, phi)

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


class GEQDSK(MagneticField):
    """Read a EFIT G-Eqdsk file for a toroidal equilibrium

    This generates a grid in cylindrical geometry

    Parameters
    ----------
    gfile : str
        Name of the file to open

    """

    def __init__(self, gfile):

        # Import utility to read G-Eqdsk files
        from boututils import geqdsk

        g = geqdsk.Geqdsk()
        g.openFile(gfile)

        # Get the range of major radius
        self.rmin = g.get('rleft')
        self.rmax = g.get('rdim') + self.rmin

        # Range of height
        self.zmin = g.get('zmid') - 0.5*g.get('zdim')
        self.zmax = g.get('zmid') + 0.5*g.get('zdim')

        print("Major radius: {0} -> {1} m".format(self.rmin, self.rmax))
        print("Height: {0} -> {1} m".format(self.zmin, self.zmax))

        # Poloidal flux
        self.psi = np.transpose(g.get('psirz'))
        nr, nz = self.psi.shape

        # Normalising factors: psi on axis and boundary
        self.psi_axis = g.get('simag')
        self.psi_bndry = g.get('sibry')

        # Current flux function f = R * Bt
        self.fpol = g.get('fpol')

        # Pressure [Pascals]
        self.p = g.get('pres')

        self.r = np.linspace(self.rmin, self.rmax, nr)
        self.z = np.linspace(self.zmin, self.zmax, nz)

        # Create a 2D spline interpolation for psi
        from scipy import interpolate
        self.psi_func = interpolate.RectBivariateSpline(self.r, self.z, self.psi)

        # Add to the attributes so that it can be written to file
        self.attributes["psi"] = lambda x,z,phi : self.psi_func(x,z,grid=False)

        # Create a normalised psi array

        self.psinorm = (self.psi - self.psi_axis) / (self.psi_bndry - self.psi_axis)

        # Need to mark areas outside the core as psinorm = 1
        # eg. around coils or in the private flux region
        # Create a boundary

        rb = g.get('rbbbs')
        zb = g.get('zbbbs')
        core_bndry = boundary.PolygonBoundaryXZ(rb, zb)

        # Get the points outside the boundary
        rxz, zxz = np.meshgrid(self.r, self.z, indexing='ij')
        outside = core_bndry.outside(rxz, 0.0, zxz)
        self.psinorm[outside] = 1.0

        self.psinorm_func = interpolate.RectBivariateSpline(self.r, self.z, self.psinorm)

        # Spline for interpolation of f = R*Bt
        psinorm = np.linspace(0.0, 1.0, nr)
        self.f_spl = interpolate.InterpolatedUnivariateSpline(psinorm, self.fpol, ext=3)
        # ext=3 specifies that boundary values are used outside range

        # Spline for interpolation of pressure
        self.p_spl = interpolate.InterpolatedUnivariateSpline(psinorm, self.p, ext=3)

        # Set boundary
        rlim = g.get('rlim')
        zlim = g.get('zlim')
        if len(rlim) > 0:
            # Create a boundary in X-Z with a polygon representation
            self.boundary = boundary.PolygonBoundaryXZ(rlim,zlim)

    def Bxfunc(self, x, z, phi):
        return -self.psi_func(x,z,dy=1,grid=False)/x

    def Bzfunc(self, x, z, phi):
        return self.psi_func(x,z,dx=1,grid=False)/x

    def Byfunc(self, x, z, phi):
        # Interpolate to get flux surface normalised psi
        psinorm = self.psinorm_func(x,z, grid=False)

        # Interpolate fpol array at values of normalised psi
        if hasattr(psinorm, "shape"):
            return np.reshape(self.f_spl(np.ravel(psinorm)),psinorm.shape) / x

        return self.f_spl(psinorm) / x  # f = R*Bt

    def Rfunc(self, x, z, phi):
        return x

    def pressure(self, x, z, phi):
        # Interpolate to get flux surface normalised psi
        psinorm = self.psinorm_func(x,z, grid=False)

        if hasattr(psinorm, "shape"):
            return np.reshape(self.p_spl(np.ravel(psinorm)),psinorm.shape)

        return self.p_spl(psinorm)

class W7X_vacuum(MagneticField):
    def __init__(self,nx=128,ny=32,nz=128,x_range=(4.05,6.55),z_range=(-1.35,1,35), phimax=2.*np.pi, configuration=0, plot_poincare=False,include_plasma_field=False, wout_file='wout_w7x.0972_0926_0880_0852_+0000_+0000.01.00jh.nc'):
        from scipy.interpolate import RegularGridInterpolator
        import numpy as np
        ## create 1D arrays of cylindrical coordinates
        r = np.linspace(x_range[0],x_range[-1],nx)
        phi = np.linspace(0,phimax,ny)
        z = np.linspace(z_range[0],z_range[-1],nz)

        ## make those 1D arrays 3D
        rarray,yarray,zarray = np.meshgrid(r,phi,z,indexing='ij')
    
        ## call vacuum field values
        b_vac  = W7X_vacuum.field_values(rarray,yarray,zarray,configuration,plot_poincare)
        Bx_vac = b_vac[0]
        By_vac = b_vac[1]
        Bz_vac = b_vac[2]

        if(include_plasma_field):
            b_plasma = W7X_vacuum.plasma_field(rarray,yarray,zarray,wout_file=wout_file)
            Bx_plasma = b_plasma[0]
            By_plasma = b_plasma[1]
            Bz_plasma = b_plasma[2] 
        else:
            Bx_plasma = 0
            By_plasma = 0
            Bz_plasma = 0
        
        Bx = Bx_vac + Bx_plasma
        By = By_vac + By_plasma
        Bz = Bz_vac + Bz_plasma
    
        # Now we have a field and regular grid in (R,Z,phi) so
        # we can get an interpolation function in 3D
        points = (r,phi,z)
            
        self.br_interp   = RegularGridInterpolator(points, Bx, bounds_error=False, fill_value=0.0)
        self.bz_interp   = RegularGridInterpolator(points, Bz, bounds_error=False, fill_value=0.0)
        self.bphi_interp = RegularGridInterpolator(points, By, bounds_error=False, fill_value=1.0)

        # if you want non-interpolated, 3D arrays, make this your return function:
        # return Bx,By,Bz
            
        # return points, br_interp, bphi_interp, bz_interp

    ################## Vacuum field service ########################
    # This uses the webservices field line tracer to get           #
    # the vacuum magnetic field given 3d arrrays for R, phi, and Z #
    # http://webservices.ipp-hgw.mpg.de/docs/fieldlinetracer.html  #
    # Only works on IPP network                                    #
    # Contact brendan.shanahan@ipp.mpg.de for questions            #
    ################################################################

    def field_values(r, phi, z, configuration=0, plot_poincare=False):
        from osa import Client
        import os.path
        import pickle
        import matplotlib.pyplot as plt

        tracer = Client('http://esb.ipp-hgw.mpg.de:8280/services/FieldLineProxy?wsdl')

        nx = r.shape[0]
        ny = phi.shape[1]
        nz = z.shape[2]
        
        # create (standardized) file name for saving/loading magnetic field.
        fname = "B.w7x.{}.{}.{}.{:.2f}-{:.2f}.{:.2f}-{:.2f}.{:.2f}-{:.2f}.dat".format(
            nx,
            ny,
            nz,
            r[0, 0, 0],
            r[-1, 0, 0],
            phi[0, 0, 0],
            phi[0, -1, 0],
            z[0, 0, 0],
            z[0, 0, -1],
        )

        if(os.path.isfile(fname)):
            print ("Saved field found, loading from: ", fname)
            with open(fname, pickle_read_mode) as f:
                Br, Bphi, Bz = pickle.load(f)

        else:
            print ("No saved field found -- (re)calculating (must be on IPP network for this to work...)")
            print ("Calculating field for Wendelstein 7-X; nx = ",nx," ny = ",ny," nz = ", nz )

            ## Create configuration objects
            config = tracer.types.MagneticConfig()
            config.configIds = configuration

            Br = np.zeros((nx,ny,nz))
            Bphi = np.ones((nx,ny,nz))
            Bz = np.zeros((nx,ny,nz))
            pos = tracer.types.Points3D()

            pos.x1 = np.ndarray.flatten(np.ones((nx,ny,nz))*r*np.cos(phi)) #x in Cartesian (real-space) 
            pos.x2 = np.ndarray.flatten(np.ones((nx,ny,nz))*r*np.sin(phi)) #y in Cartesian (real-space) 
            pos.x3 = np.ndarray.flatten(z)                                 #z in Cartesian (real-space) 

            ## Call tracer service
            res = tracer.service.magneticField(pos, config)

            ## Reshape to 3d array
            Bx = np.ndarray.reshape(np.asarray(res.field.x1),(nx,ny,nz))
            By = np.ndarray.reshape(np.asarray(res.field.x2),(nx,ny,nz))
            Bz = np.ndarray.reshape(np.asarray(res.field.x3), (nx,ny,nz))

            ## Convert to cylindrical coordinates
            Br = Bx*np.cos(phi) + By*np.sin(phi)
            Bphi = -Bx*np.sin(phi) + By*np.cos(phi)
        
            ## Save so we don't have to do this every time.
            with open(fname, pickle_write_mode) as f:
                pickle.dump([Br,Bphi,Bz],f)
            
        if(plot_poincare):
            ## Poincare plot as done on the web services
            ## Independent of the previously-made field.
        
            print( "Making poincare plot (only works on IPP network)...")
            ## Create configuration objects
            config = tracer.types.MagneticConfig()
            config.configIds = configuration
            pos = tracer.types.Points3D()
            
            pos.x1 = np.linspace(5.6, 6.2, 80)
            pos.x2 = np.zeros(80)
            pos.x3 = np.zeros(80)
            
            poincare = tracer.types.PoincareInPhiPlane()
            poincare.numPoints = 200
            poincare.phi0 = [0.] ## This is where the poincare plane is (bean=0, triangle = pi/5.)

            task = tracer.types.Task()
            task.step = 0.2
            task.poincare = poincare

            res = tracer.service.trace(pos, config, task, None, None)

            for i in range(0, len(res.surfs)):
                plt.scatter(res.surfs[i].points.x1, res.surfs[i].points.x3, color="black", s=0.1)
                plt.show()
            
        return  Br, Bphi, Bz


    ################## Plasma field service ########################
    # This uses EXTENDER via the webservices to get                #
    # the magnetic field from the plasma given 3d arrrays          #
    # for R, phi, and Z                                            #
    # http://webservices.ipp-hgw.mpg.de/docs/extender.html         #
    # Only works on IPP network                                    #
    # Contact brendan.shanahan@ipp.mpg.de for questions            #
    ################################################################

    def plasma_field(r, phi, z, wout_file='wout.nc'):
        from osa import Client
        import os.path
        import pickle
        
        cl = Client('http://esb.ipp-hgw.mpg.de:8280/services/Extender?wsdl')
        vmecURL = "http://svvmec1.ipp-hgw.mpg.de:8080/vmecrest/v1/w7x_ref_1/wout.nc"

        nx = r.shape[0]
        ny = phi.shape[1]
        nz = z.shape[2]

        # create (standardized) file name for saving/loading magnetic field.
        fname = "B.w7x_plasma_field.{}.{}.{}.{:.2f}-{:.2f}.{:.2f}-{:.2f}.{:.2f}-{:.2f}.dat".format(
            nx,
            ny,
            nz,
            r[0, 0, 0],
            r[-1, 0, 0],
            phi[0, 0, 0],
            phi[0, -1, 0],
            z[0, 0, 0],
            z[0, 0, -1],
        )

        if(os.path.isfile(fname)):
            print ("Saved field found, loading from: ", fname)
            with open(fname, pickle_read_mode) as f:
                Br, Bphi, Bz = pickle.load(f)
        else:
            print ("No saved plasma field found -- (re)calculating (must be on IPP network for this to work...)")
            print ("Calculating plasma field for Wendelstein 7-X; nx = ",nx," ny = ",ny," nz = ", nz )
            print ("This part takes AGES... estimate: ",nx*ny*nz/52380., " minutes.")

            points = cl.types.Points3D()

            ## Extender uses cylindrical coordinates, no need to convert, just flatten.
            points.x1 = np.ndarray.flatten(r) #x in Cylindrical 
            points.x2 = np.ndarray.flatten(phi) #y in Cylindrical 
            points.x3 = np.ndarray.flatten(z)  #z in Cylindrical

            ## call EXTENDER on web services
            # if not (os.path.isfile(wout_file)):
            plasmafield = cl.service.getPlasmaField(None, vmecURL, points, None)
            # else:
            # plasmafield = cl.service.getPlasmaField(wout, None, points, None)
                
            ## Reshape to 3d array
            Br = np.ndarray.reshape(np.asarray(plasmafield.x1),(nx,ny,nz))
            Bphi = np.ndarray.reshape(np.asarray(plasmafield.x2),(nx,ny,nz))
            Bz = np.ndarray.reshape(np.asarray(plasmafield.x3), (nx,ny,nz))
        
            ## Save so we don't have to do this every time.
            with open(fname, pickle_write_mode) as f:
                pickle.dump([Br,Bphi,Bz],f)
            
        return  Br, Bphi, Bz
    
    def magnetic_axis(self, phi_axis=0,configuration=0):
        from osa import Client
        
        tracer = Client('http://esb.ipp-hgw.mpg.de:8280/services/FieldLineProxy?wsdl')
        
        config = tracer.types.MagneticConfig()
        config.configIds = configuration
        settings = tracer.types.AxisSettings()
        res = tracer.service.findAxis(0.05, config, settings)
    
        magnetic_axis_x = np.asarray(res.axis.vertices.x1)# (m) 
        magnetic_axis_y = np.asarray(res.axis.vertices.x2)# (m) -- REAL SPACE from an arbitrary start point
        magnetic_axis_z = np.asarray(res.axis.vertices.x3)# (m)
        magnetic_axis_rmaj = np.sqrt(magnetic_axis_x**2 + magnetic_axis_y**2 + magnetic_axis_z**2)

        magnetic_axis_r = np.sqrt(np.asarray(magnetic_axis_x)**2 + np.asarray(magnetic_axis_y**2))
        magnetic_axis_phi = np.arctan(magnetic_axis_y/magnetic_axis_x)

        index = np.where((magnetic_axis_phi >= 0.97*phi_axis ) & (magnetic_axis_phi <= 1.03*phi_axis) )
        index = index[0]

        return np.asarray([magnetic_axis_r[index],magnetic_axis_z[index]])[:,0]

    def Bxfunc(self, x, z, phi):
        phi = np.mod(phi, 2.*np.pi)
        return self.br_interp((x,phi,z))

    def Bzfunc(self, x, z, phi):
        phi = np.mod(phi, 2.*np.pi)
        return self.bz_interp((x,phi,z))

    def Byfunc(self, x, z, phi):
        phi = np.mod(phi, 2.*np.pi)
        # Interpolate to get flux surface normalised psi
        return self.bphi_interp((x,phi,z))
    
    def Rfunc(self, x, z, phi):
        phi = np.mod(phi, 2.*np.pi)
        return x

class W7X_VMEC(MagneticField):
    def __init__(self,nx=512,ny=32,nz=512,x_range=(4.05,6.55),z_range=(-1.35,1,35), phi_range=(0,2*np.pi), vmec_id='w7x_ref_171'):
        from scipy.interpolate import RegularGridInterpolator
        self.nx = nx
        self.ny = ny
        self.nz = nz
        ## create 1D arrays of cylindrical coordinates
        r = np.linspace(x_range[0], x_range[-1], nx)
        phi = np.linspace(phi_range[0], phi_range[-1],ny)
        z = np.linspace(z_range[0], z_range[-1], nz)

        ## make those 1D arrays 3D
        rarray,yarray,zarray = np.meshgrid(r,phi,z,indexing='ij')
    
        ## call vacuum field values
        b_vmec  = self.field_values(rarray,yarray,zarray,vmec_id)
        Bx_vmec = b_vmec[0]
        By_vmec = b_vmec[1]
        Bz_vmec = b_vmec[2]

        # Now we have a field and regular grid in (R,Z,phi) so
        # we can get an interpolation function in 3D
        points = (r,phi,z)
            
        self.br_interp   = RegularGridInterpolator(points, Bx_vmec, bounds_error=False, fill_value=0.0)
        self.bz_interp   = RegularGridInterpolator(points, Bz_vmec, bounds_error=False, fill_value=0.0)
        self.bphi_interp = RegularGridInterpolator(points, By_vmec, bounds_error=False, fill_value=1.0)

    def field_values(self, r,phi,z, vmec_id='w7x_ref_171'):
        from osa import Client
        vmec = Client('http://esb:8280/services/vmec_v5?wsdl')

        pos = vmec.types.Points3D()

        pos.x1 = np.ndarray.flatten(np.ones((self.nx,self.ny,self.nz))*r*np.cos(phi)) #x in Cartesian (real-space) 
        pos.x2 = np.ndarray.flatten(np.ones((self.nx,self.ny,self.nz))*r*np.sin(phi)) #y in Cartesian (real-space) 
        pos.x3 = np.ndarray.flatten(z)                                 #z in Cartesian (real-space) 
        b = vmec.service.magneticField(str(vmec_id), pos)

        ## Reshape to 3d array
        Bx = np.ndarray.reshape(np.asarray(b.field.x1),(self.nx,self.ny,self.nz))
        By = np.ndarray.reshape(np.asarray(b.field.x2),(self.nx,self.ny,self.nz))
        Bz = np.ndarray.reshape(np.asarray(b.field.x3), (self.nx,self.ny,self.nz))
        
        ## Convert to cylindrical coordinates
        Br = Bx*np.cos(phi) + By*np.sin(phi)
        Bphi = -Bx*np.sin(phi) + By*np.cos(phi)

        return Br, Bphi, Bz

    def Bxfunc(self, x, z, phi):
        phi = np.mod(phi, 2.*np.pi)
        return self.br_interp((x,phi,z))

    def Bzfunc(self, x, z, phi):
        phi = np.mod(phi, 2.*np.pi)
        return self.bz_interp((x,phi,z))

    def Byfunc(self, x, z, phi):
        phi = np.mod(phi, 2.*np.pi)
        # Interpolate to get flux surface normalised psi
        return self.bphi_interp((x,phi,z))
    
    def Rfunc(self, x, z, phi):
        phi = np.mod(phi, 2.*np.pi)
        return x
