""" Functions for calculating sources for the
    Method of Manufactured Solutions (MMS)

"""
from __future__ import print_function
from __future__ import division

from sympy import symbols, cos, sin, diff, sqrt, pi, simplify, trigsimp, Wild

from numpy import arange, zeros

# Constants
qe = 1.602e-19
Mp = 1.67262158e-27
mu0 = 4.e-7*3.141592653589793

# Define symbols

x = symbols('x')
y = symbols('y')
z = symbols('z')
t = symbols('t')

class Metric(object):
    def __init__(self):
        # Create an identity metric
        self.x = symbols('x\'')
        self.y = symbols('y\'')
        self.z = symbols('z\'')

        self.g11 = self.g22 = self.g33 = 1.0
        self.g12 = self.g23 = self.g13 = 0.0

        self.g_11 = self.g_22 = self.g_33 = 1.0
        self.g_12 = self.g_23 = self.g_13 = 0.0

        self.J = 1.0
        self.B = 1.0

identity = Metric()

# Basic differencing
def ddt(f):
    """Time derivative"""
    return diff(f, t)


def DDX(f, metric = identity):
    return diff(f, metric.x)

def DDY(f, metric = identity):
    return diff(f, metric.y)

def DDZ(f, metric = identity):
    return diff(f, metric.z)


def D2DX2(f, metric = identity):
    return diff(f, metric.x, 2)

def D2DY2(f, metric = identity):
    return diff(f, metric.y, 2)

def D2DZ2(f, metric = identity):
    return diff(f, metric.z, 2)


def D2DXDY(f, metric = identity):
    message = "* WARNING: D2DXDY is currently not set in BOUT++."+\
              " Check src/sys/derivs.cxx if situation has changed. *"
    print("\n"*3)
    print("*"*len(message))
    print(message)
    print("*"*len(message))
    print("\n"*3)
    return DDX(DDY(f, metric), metric)

def D2DXDZ(f, metric = identity):
    return DDX(DDZ(f, metric), metric)

def D2DYDZ(f, metric = identity):
    return DDY(DDZ(f, metric), metric)

# Operators

def bracket(f, g, metric = identity):
    """
    Calculates [f,g] symbolically
    """

    dfdx = diff(f, metric.x)
    dfdz = diff(f, metric.z)

    dgdx = diff(g, metric.x)
    dgdz = diff(g, metric.z)

    return dfdz * dgdx - dfdx * dgdz

def b0xGrad_dot_Grad(phi, A, metric = identity):
    """
    Perpendicular advection operator, including
    derivatives in y

    Note: If y derivatives are neglected, then this reduces
    to bracket(f, g, metric) * metric.B
    (in a Clebsch coordinate system)

    """
    dpdx = DDX(phi, metric)
    dpdy = DDY(phi, metric)
    dpdz = DDZ(phi, metric)

    vx = metric.g_22*dpdz - metric.g_23*dpdy;
    vy = metric.g_23*dpdx - metric.g_12*dpdz;
    vz = metric.g_12*dpdy - metric.g_22*dpdx;

    return (+ vx*DDX(A, metric)
            + vy*DDY(A, metric)
            + vz*DDZ(A, metric) ) / (metric.J*sqrt(metric.g_22))

def Delp2(f, metric = identity, all_terms=True):
    """ Laplacian in X-Z

    If all_terms is false then first derivative terms are neglected.
    By default all_terms is true, but can be disabled
    in the BOUT.inp file (laplace section)

    """
    d2fdx2 = diff(f, metric.x, 2)
    d2fdz2 = diff(f, metric.z, 2)
    d2fdxdz = diff(f, metric.x, metric.z)

    result = metric.g11*d2fdx2 + metric.g33*d2fdz2 + 2.*metric.g13*d2fdxdz

    if all_terms:
        G1 = (DDX(metric.J*metric.g11, metric) + DDY(metric.J*metric.g12, metric) + DDZ(metric.J*metric.g13, metric)) / metric.J
        G3 = (DDX(metric.J*metric.g13, metric) + DDY(metric.J*metric.g23, metric) + DDZ(metric.J*metric.g33, metric)) / metric.J
        result += G1 * diff(f, metric.x) + G3 * diff(f, metric.z)

    return result

def Delp4(f, metric = identity):
    d4fdx4 = diff(f, metric.x, 4)
    d4fdz4 = diff(f, metric.z, 4)

    return d4fdx4 + d4fdz4

def Grad_par(f, metric = identity):
    """The parallel gradient"""
    return diff(f, metric.y) / sqrt(metric.g_22)

def Vpar_Grad_par(v, f, metric = identity):
    """Parallel advection operator v*grad_||(f)"""
    return v * Grad_par(f, metric=metric)

def Div_par(f, metric=identity):
    '''
    Divergence of magnetic field aligned vector v = \bhat f
    \nabla \cdot (\bhat f) = 1/J \partial_y (f/B)
    = B Grad_par(f/B)
    '''
    return metric.B*Grad_par(f/metric.B, metric)

def Laplace(f, metric=identity):
    """The full Laplace operator"""
    G1 = (DDX(metric.J*metric.g11, metric) + DDY(metric.J*metric.g12, metric) + DDZ(metric.J*metric.g13, metric)) / metric.J
    G2 = (DDX(metric.J*metric.g12, metric) + DDY(metric.J*metric.g22, metric) + DDZ(metric.J*metric.g23, metric)) / metric.J
    G3 = (DDX(metric.J*metric.g13, metric) + DDY(metric.J*metric.g23, metric) + DDZ(metric.J*metric.g33, metric)) / metric.J

    result  = G1*DDX(f, metric) + G2*DDY(f, metric) + G3*DDZ(f, metric)\
	      + metric.g11*D2DX2(f, metric) + metric.g22*D2DY2(f, metric) + metric.g33*D2DZ2(f, metric)\
	      + 2.0*(metric.g12*D2DXDY(f, metric) + metric.g13*D2DXDZ(f, metric) + metric.g23*D2DYDZ(f, metric))

    return result

def Laplace_par(f, metric=identity):
    """
    Div( b (b.Grad(f) ) ) = (1/J) d/dy ( J/g_22 * df/dy )
    """
    return diff( (metric.J/metric.g_22)*diff(f, metric.y), metric.y)/ metric.J

def Laplace_perp(f, metric=identity):
    """
    The perpendicular Laplace operator

    Laplace_perp = Laplace - Laplace_par
    """
    return Laplace(f, metric) - Laplace_par(f, metric)

# Convert expression to string

def trySimplify(expr):
    """
    Tries to simplify an expression
    """
    try:
        return simplify(expr)
    except ValueError:
        return expr

def exprToStr(expr):
    """ Convert a sympy expression to a string for BOUT++ input
    """

    s = str(expr).replace("**", "^") # Replace exponent operator

    # Try to remove lots of 1.0*...
    s = s.replace("(1.0*", "(")
    s = s.replace(" 1.0*", " ")

    return s

def exprMag(expr):
    """
    Estimate the magnitude of an expression

    """

    # Replace all sin, cos with 1
    any = Wild('a') # Wildcard
    expr = expr.replace(sin(any), 1.0)
    expr = expr.replace(cos(any), 1.0)

    # Pick maximum values of x,y,z
    expr = expr.subs(x, 1.0)
    expr = expr.subs(y, 2.*pi)
    expr = expr.subs(z, 2.*pi)

    return expr.evalf()

##################################

class SimpleTokamak(object):
    """
    Simple tokamak

    NOTE: This is NOT an equilibrium calculation. The input
    is intended solely for testing with MMS
    """
    def __init__(self, R = 2, Bt = 1.0, eps = 0.1, dr=0.02, q = lambda x:2+x**2):
        """
        R    - Major radius [m]

        Bt   - Toroidal field [T]

        eps  - Inverse aspect ratio

        dr   - Width of the radial region [m]

        q(x) - A function which returns the safety factor
               as a function of x in range [0,1]


        Coordinates:
        x - Radial, [0,1]
        y - Poloidal, [0,2pi]. Origin is at inboard midplane.


        """
        # X has a range [0,1], and y [0,2pi]
        #x, y = symbols("x y")

        self.x = x
        self.y = y

        self.R = R

        self.dr = dr

        # Minor radius
        self.r = R * eps

        # Get safety factor
        self.q = q(x)

        # Toroidal angle of a field-line as function
        # of poloidal angle y
        self.zShift = self.q*(y + eps * sin(y))

        # Field-line pitch
        self.nu = self.q*(1 + eps*cos(y)) #diff(self.zShift, y)

        # Coordinates of grid points
        self.Rxy = R - self.r * cos(y)
        self.Zxy = self.r * sin(y)

        # Poloidal arc length
        self.hthe = self.r + 0.*x

        # Toroidal magnetic field
        self.Btxy = Bt * R / self.Rxy

        # Poloidal magnetic field
        self.Bpxy = self.Btxy * self.hthe / (self.nu * self.Rxy)

        # Total magnetic field
        self.Bxy = sqrt(self.Btxy**2 + self.Bpxy**2)

        # Approximate poloidal field for radial width calculation
        Bp0 = Bt * self.r / (q(0.5) * R)
        print("Bp0 = %e" % Bp0)

        # dx = Bp * R * dr  -- width of the box in psi space
        self.psiwidth = Bp0 * R * dr
        print("psi width = %e" % self.psiwidth)

        # Integrated shear
        self.sinty = diff(self.zShift, x) / self.psiwidth

        # Extra expressions to add to grid file
        self._extra = {}

    def add(self, expr, name):
        """
        Add an additional expression to be written to the grid files

        """
        self._extra[name] = expr


    def write(self, nx, ny, output, MXG=2):
        """
        Outputs a tokamak shape to a grid file

        nx - Number of radial grid points, not including guard cells
        ny - Number of poloidal (parallel) grid points
        output - boututils.datafile object, e.g., an open netCDF file
        MXG, Number of guard cells in the x-direction
        """       

        ngx = nx + 2*MXG
        ngy = ny

        # Create an x and y grid to evaluate expressions on
        xarr = (arange(nx + 2*MXG) - MXG + 0.5) / nx
        yarr = 2.*pi*arange(ny)/ny

        output.write("nx", ngx)
        output.write("ny", ngy)

        dx = self.psiwidth / nx  + 0.*self.x
        dy = 2.*pi / ny + 0.*self.x

        for name, var in [ ("dx", dx),
                           ("dy", dy),
                           ("Rxy", self.Rxy),
                           ("Zxy", self.Zxy),
                           ("Btxy", self.Btxy),
                           ("Bpxy", self.Bpxy),
                           ("Bxy", self.Bxy),
                           ("hthe", self.hthe),
                           ("sinty", self.sinty),
                           ("zShift", self.zShift)]:

            # Note: This is slow, and could be improved using something like lambdify
            values = zeros([ngx, ngy])
            for i, x in enumerate(xarr):
                for j, y in enumerate(yarr):
                    values[i,j] = var.evalf(subs={self.x:x, self.y:y})

            output.write(name, values)

        for name, var in list(self._extra.items()):
            values = zeros([ngx, ngy])
            for i, x in enumerate(xarr):
                for j, y in enumerate(yarr):
                    values[i,j] = var.evalf(subs={self.x:x, self.y:y})

            output.write(name, values)

        shiftAngle = zeros(ngx)
        for i, x in enumerate(xarr):
            shiftAngle[i] = 2.*pi*self.q.evalf(subs={self.x:x})

        output.write("ShiftAngle", shiftAngle)

    def metric(self):
        """
        Returns an analytic metric tensor
        """
        m = Metric()

        # Set symbols for x and y directions
        m.x = self.x
        m.y = self.y

        # Calculate metric tensor

        m.g11 = (self.Rxy * self.Bpxy)**2
        m.g22 = 1./self.hthe**2
        m.g33 = self.sinty**2*m.g11 + self.Bxy**2/m.g11
        m.g12 = 0.0*x
        m.g13 = -self.sinty*m.g11
        m.g23 = -self.Btxy / (self.hthe * self.Bpxy * self.R)

        m.g_11 = 1./m.g11 + (self.sinty*self.Rxy)**2
        m.g_22 = (self.Bxy * self.hthe / self.Bpxy)**2
        m.g_33 = self.Rxy**2
        m.g_12 = self.Btxy*self.hthe*self.sinty*self.Rxy / self.Bpxy
        m.g_13 = self.sinty*self.Rxy**2
        m.g_23 = self.Btxy*self.hthe*self.Rxy / self.Bpxy

        m.J = self.hthe / self.Bpxy
        m.B = self.Bxy

        # Convert all "x" symbols from [0,1] into flux
        m.Lx = self.psiwidth
        xsub = m.x / self.psiwidth

        m.g11 = m.g11.subs(x, xsub)
        m.g22 = m.g22.subs(x, xsub)
        m.g33 = m.g33.subs(x, xsub)
        m.g12 = m.g12.subs(x, xsub)
        m.g13 = m.g13.subs(x, xsub)
        m.g23 = m.g23.subs(x, xsub)

        m.g_11 = m.g_11.subs(x, xsub)
        m.g_22 = m.g_22.subs(x, xsub)
        m.g_33 = m.g_33.subs(x, xsub)
        m.g_12 = m.g_12.subs(x, xsub)
        m.g_13 = m.g_13.subs(x, xsub)
        m.g_23 = m.g_23.subs(x, xsub)

        m.J = m.J.subs(x, xsub)
        m.B = m.B.subs(x, xsub)

        return m

##########################
# Shaped tokamak

class ShapedTokamak(object):
    def __init__(self, Rmaj=6.0, rmin=2.0, dr=0.1, kappa=1.0, delta=0.0, b=0.0, ss=0.0, Bt0=1.0, Bp0 = 0.2):
        """
        Rmaj  - Major radius [m]
        rmin  - Minor radius [m]
        dr    - Radial width of region [m]

        kappa - Ellipticity, 1 for a circle
        delta - Triangularity, 0 for circle
        b     - Indentation ("bean" shape), 0 for circle

        ss    - Shafranov shift [m]

        Bt0   - Toroidal magnetic field on axis [T]. Varies as 1/R
        Bp0   - Poloidal field at outboard midplane [T]

        Outputs
        -------

        Assigns member variables

        x, y    - Symbols for x and y coordinates

        R (x,y)
        Z (x,y)

        """

        # X has a range [0,1], and y [0,2pi]
        x, y = symbols("x y")

        # Minor radius as function of x
        rminx = rmin + (x-0.5)*dr

        # Analytical expression for R and Z coordinates as function of x and y
        Rxy = Rmaj - b + (rminx + b*cos(y))*cos(y + delta*sin(y)) + ss*(0.5-x)*(dr/rmin)
        Zxy = kappa * rminx * sin(y)

        # Toroidal magnetic field
        Btxy = Bt0 * Rmaj / Rxy

        # Poloidal field. dx constant, so set poloidal field
        # at outboard midplane (y = 0)
        # NOTE: Approximate calculation

        # Distance between flux surface relative to outboard midplane.
        expansion = (1 - (old_div(ss,rmin))*cos(y))/(1 - (ss/rmin))

        Bpxy = Bp0 * ((Rmaj + rmin) / Rxy) / expansion

        # Calculate hthe
        hthe = sqrt(diff(Rxy, y)**2 + diff(Zxy, y)**2)
        try:
            hthe = trigsimp(hthe)
        except ValueError:
            pass

        # Field-line pitch
        nu = Btxy * hthe / (Bpxy * Rxy)

        # Shift angle
        # NOTE: Since x has a range [0,1] this could be done better
        # than ignoring convergence conditions
        self.zShift = integrate(nu, y, conds='none')

        # Safety factor
        self.shiftAngle = self.zShift.subs(y, 2*pi) - self.zShift.subs(y, 0)

        # Integrated shear
        self.I = diff(self.zShift, x)

        self.x = x
        self.y = y

        self.R = Rxy
        self.Z = Zxy

        self.Bt = Btxy
        self.Bp = Bpxy
        self.B = sqrt(Btxy**2 + Bpxy**2)

        self.hthe = hthe

    def write(self, nx, ny, filename, MXG=2):
        """
        Outputs a tokamak shape to a grid file

        nx - Number of radial grid points, not including guard cells
        ny - Number of poloidal (parallel) grid points
        output - boututils.datafile object, e.g., an open netCDF file
        MXG, Number of guard cells in the x-direction
        """ 

        ngx = nx + 2*MXG
        ngy = ny

        # Create an x and y grid to evaluate expressions on
        xarr = (arange(nx + 2*MXG) - MXG + 0.5) / nx
        yarr = 2.*pi*arange(ny)/ny

        Rxy = zeros([ngx, ngy])
        Zxy = zeros([ngx, ngy])

        Btxy = zeros([ngx, ngy])
        Bpxy = zeros([ngx, ngy])

        hthe = zeros([ngx, ngy])


        I = zeros([ngx, ngy])

        # Note: This is slow, and could be improved using something like lambdify
        for i, x in enumerate(xarr):
            for j, y in enumerate(yarr):
                Rxy[i,j] = self.R.evalf(subs={self.x:x, self.y:y})
                Zxy[i,j] = self.Z.evalf(subs={self.x:x, self.y:y})

                Btxy[i,j] = self.Bt.evalf(subs={self.x:x, self.y:y})
                Bpxy[i,j] = self.Bp.evalf(subs={self.x:x, self.y:y})

                hthe[i,j] = self.hthe.evalf(subs={self.x:x, self.y:y})


            plt.plot(Rxy[i,:], Zxy[i,:])
        plt.show()

        Bxy = sqrt(Btxy**2 + Bpxy**2)

    def metric(self):
        """
        Returns an analytic metric tensor
        """
        m = Metric()

        # Set symbols for x and y directions
        m.x = self.x
        m.y = self.y

        # Calculate metric tensor

        m.g11 = (self.R * self.Bp)**2
        m.g22 = 1./self.hthe**2
        m.g33 = self.I**2*m.g11 + self.B**2 / m.g11
        m.g12 = 0.0
        m.g13 = -self.I*m.g11
        m.g23 = -self.Bt / (self.hthe * self.Bp * self.R)

        m.g_11 = 1./m.g11 + (self.I*self.R)**2
        m.g_22 = (self.B * self.hthe / self.Bpxy)**2
        m.g_33 = self.R**2
        m.g_12 = self.Bt*self.hthe*self.I*self.R / self.Bp
        m.g_13 = self.I*self.R**2
        m.g_23 = self.Bt*self.hthe*self.R / self.Bp

        m.J = self.hthe / self.Bp
        m.B = self.B

        return m


