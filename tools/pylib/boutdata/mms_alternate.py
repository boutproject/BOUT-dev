""" Functions for calculating sources for the
    Method of Manufactured Solutions (MMS)

"""
from __future__ import print_function
from __future__ import division
#from builtins import str
#from builtins import object

from sympy import symbols, cos, sin, diff, sqrt, pi, simplify, trigsimp, Wild, integrate

from numpy import arange, newaxis
#from numpy import arange, zeros

global metric

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
        self.x = x
        self.y = y
        self.z = z

        self.g11 = self.g22 = self.g33 = 1.0
        self.g12 = self.g23 = self.g13 = 0.0

        self.g_11 = self.g_22 = self.g_33 = 1.0
        self.g_12 = self.g_23 = self.g_13 = 0.0

        self.J = 1.0
        self.B = 1.0

metric = Metric()

# Basic differencing
def ddt(f):
    """Time derivative"""
    return diff(f, t)


def DDX(f):
    return diff(f, metric.x)

def DDY(f):
    return diff(f, metric.y)

def DDZ(f):
    return diff(f, metric.z)


def D2DX2(f):
    return diff(f, metric.x, 2)

def D2DY2(f):
    return diff(f, metric.y, 2)

def D2DZ2(f):
    return diff(f, metric.z, 2)


def D2DXDY(f):
    message = "* WARNING: D2DXDY is currently not set in BOUT++."+\
              " Check src/sys/derivs.cxx if situation has changed. *"
    print("\n"*3)
    print("*"*len(message))
    print(message)
    print("*"*len(message))
    print("\n"*3)
    return DDX(DDY(f))

def D2DXDZ(f):
    return DDX(DDZ(f))

def D2DYDZ(f):
    return DDY(DDZ(f))

# Operators

def bracket(f, g):
    """
    Calculates [f,g] symbolically
    """

    dfdx = diff(f, metric.x)
    dfdz = diff(f, metric.z)

    dgdx = diff(g, metric.x)
    dgdz = diff(g, metric.z)

    return dfdz * dgdx - dfdx * dgdz

def b0xGrad_dot_Grad(phi, A):
    """
    Perpendicular advection operator, including
    derivatives in y

    Note: If y derivatives are neglected, then this reduces
    to bracket(f, g) * metric.B
    (in a Clebsch coordinate system)

    """
    dpdx = DDX(phi)
    dpdy = DDY(phi)
    dpdz = DDZ(phi)

    vx = metric.g_22*dpdz - metric.g_23*dpdy;
    vy = metric.g_23*dpdx - metric.g_12*dpdz;
    vz = metric.g_12*dpdy - metric.g_22*dpdx;

    return (+ vx*DDX(A)
            + vy*DDY(A)
            + vz*DDZ(A) ) / (metric.J*sqrt(metric.g_22))

def Delp2(f, all_terms=True):
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
        G1 = (DDX(metric.J*metric.g11) + DDY(metric.J*metric.g12) + DDZ(metric.J*metric.g13)) / metric.J
        G3 = (DDX(metric.J*metric.g13) + DDY(metric.J*metric.g23) + DDZ(metric.J*metric.g33)) / metric.J
        result += G1 * diff(f, metric.x) + G3 * diff(f, metric.z)

    return result

def Delp4(f):
    d4fdx4 = diff(f, metric.x, 4)
    d4fdz4 = diff(f, metric.z, 4)

    return d4fdx4 + d4fdz4

def Grad_par(f):
    """The parallel gradient"""
    return diff(f, metric.y) / sqrt(metric.g_22)

def Vpar_Grad_par(v, f):
    """Parallel advection operator v*grad_||(f)"""
    return v * Grad_par(f)

def Div_par(f):
    '''
    Divergence of magnetic field aligned vector v = \bhat f
    \nabla \cdot (\bhat f) = 1/J \partial_y (f/B)
    = B Grad_par(f/B)
    '''
    return metric.B*Grad_par(f/metric, metric)

def Laplace(f):
    """The full Laplace operator"""
    G1 = (DDX(metric.J*metric.g11) + DDY(metric.J*metric.g12) + DDZ(metric.J*metric.g13)) / metric.J
    G2 = (DDX(metric.J*metric.g12) + DDY(metric.J*metric.g22) + DDZ(metric.J*metric.g23)) / metric.J
    G3 = (DDX(metric.J*metric.g13) + DDY(metric.J*metric.g23) + DDZ(metric.J*metric.g33)) / metric.J

    result  = G1*DDX(f) + G2*DDY(f) + G3*DDZ(f)\
	      + metric.g11*D2DX2(f) + metric.g22*D2DY2(f) + metric.g33*D2DZ2(f)\
	      + 2.0*(metric.g12*D2DXDY(f) + metric.g13*D2DXDZ(f) + metric.g23*D2DYDZ(f))

    return result

def Laplace_par(f):
    """
    Div( b (b.Grad(f) ) ) = (1/J) d/dy ( J/g_22 * df/dy )
    """
    return diff( (metric.J/metric.g_22)*diff(f, metric.y), metric.y)/ metric.J

def Laplace_perp(f):
    """
    The perpendicular Laplace operator

    Laplace_perp = Laplace - Laplace_par
    """
    return Laplace(f) - Laplace_par(f)

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

class BaseTokamak(object):
    """
    Virtual base class defining various useful functions for child *Tokamak classes

    Child class __init__ methods must define the member fields:
        self.x
        self.y

        # Major radius of axis
        self.R

        self.dr

        # Minor radius
        self.r

        # Safety factor
        self.q

        # Toroidal angle of a field-line as function
        # of poloidal angle y
        self.zShift

        # Field-line pitch
        self.nu

        # Coordinates of grid points
        self.Rxy
        self.Zxy

        # Poloidal arc length
        self.hthe

        # Toroidal magnetic field
        self.Btxy

        # Poloidal magnetic field
        self.Bpxy

        # Total magnetic field
        self.Bxy

        # dx = Bp * R * dr  -- width of the box in psi space
        self.psiwidth

        # Integrated shear
        self.sinty

        # Extra expressions to add to grid file
        self._extra
    """

    def __init__(self):
        raise ValueError("Error: An instance of BaseTokamak should never be created, use a derived class like SimpleTokamak or ShapedTokamak instead")

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
        xarr = xarr[:,newaxis]
        yarr = 2.*pi*arange(ny)/ny
        yarr = yarr[newaxis,:]

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

            varfunc = lambdify([x, y], var)
            values = varfunc(xarr, yarr)
            ## Note: This is slow, and could be improved using something like lambdify
            #values = zeros([ngx, ngy])
            #for i, x in enumerate(xarr):
            #    for j, y in enumerate(yarr):
            #        values[i,j] = var.evalf(subs={self.x:x, self.y:y})

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
        Calculates an analytic metric tensor
        """

        # Set symbols for x and y directions
        metric.x = self.x
        metric.y = self.y

        # Calculate metric tensor
        sbp = 1. # sign of Bp
        if (self.Bpxy.evalf(subs={x:0., y:0.}) < 0.):
            sbp = -1.

        metric.g11 = (self.Rxy * self.Bpxy)**2
        metric.g22 = 1/self.hthe**2
        metric.g33 = self.sinty**2*metric.g11 + self.Bxy**2/metric.g11
        metric.g12 = 0*x
        metric.g13 = -self.sinty*metric.g11
        metric.g23 = -sbp*self.Btxy / (self.hthe * self.Bpxy * self.Rxy)

        metric.g_11 = 1./metric.g11 + (self.sinty*self.Rxy)**2
        metric.g_22 = (self.Bxy * self.hthe / self.Bpxy)**2
        metric.g_33 = self.Rxy**2
        metric.g_12 = sbp*self.Btxy*self.hthe*self.sinty*self.Rxy / self.Bpxy
        metric.g_13 = self.sinty*self.Rxy**2
        metric.g_23 = sbp*self.Btxy*self.hthe*self.Rxy / self.Bpxy

        metric.J = self.hthe / self.Bpxy
        metric.B = self.Bxy

        # Convert all "x" symbols from [0,1] into flux
        metric.Lx = self.psiwidth
        xsub = metric.x / self.psiwidth

        #metric.g11 = metric.g11.subs(x, xsub)
        #metric.g22 = metric.g22.subs(x, xsub)
        #metric.g33 = metric.g33.subs(x, xsub)
        #metric.g12 = metric.g12.subs(x, xsub)
        #metric.g13 = metric.g13.subs(x, xsub)
        #metric.g23 = metric.g23.subs(x, xsub)

        #metric.g_11 = metric.g_11.subs(x, xsub)
        #metric.g_22 = metric.g_22.subs(x, xsub)
        #metric.g_33 = metric.g_33.subs(x, xsub)
        #metric.g_12 = metric.g_12.subs(x, xsub)
        #metric.g_13 = metric.g_13.subs(x, xsub)
        #metric.g_23 = metric.g_23.subs(x, xsub)

        #metric.J = metric.J.subs(x, xsub)
        #metric.B = metric.B.subs(x, xsub)

        self.metric_is_set = True

    def print_mesh(self):
        """
        Prints the metrics to stdout to be copied to a BOUT.inp file
        """
        if not self.metric_is_set:
            raise ValueError("Error: metric has not been calculated yet, so cannot print")
        print("g11 = "+exprToStr(metric.g11))
        print("g22 = "+exprToStr(metric.g22))
        print("g33 = "+exprToStr(metric.g33))
        print("g12 = "+exprToStr(metric.g12))
        print("g13 = "+exprToStr(metric.g13))
        print("g23 = "+exprToStr(metric.g23))
        print("g_11 = "+exprToStr(metric.g_11))
        print("g_22 = "+exprToStr(metric.g_22))
        print("g_33 = "+exprToStr(metric.g_33))
        print("g_12 = "+exprToStr(metric.g_12))
        print("g_13 = "+exprToStr(metric.g_13))
        print("g_23 = "+exprToStr(metric.g_23))
        print("J = "+exprToStr(metric.J))
        print("Bxy = "+exprToStr(metric.B))

##################################

class SimpleTokamak(BaseTokamak):
    """
    Simple tokamak

    NOTE: This is NOT an equilibrium calculation. The input
    is intended solely for testing with MMS
    """
    def __init__(self, R = 2, Bt = 1.0, eps = 0.1, dr=0.02, q = lambda x:2+x**2):
        """
        R    - Major radius [metric]

        Bt   - Toroidal field [T]

        eps  - Inverse aspect ratio

        dr   - Width of the radial region [metric]

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
        #print("Bp0 = %e" % Bp0)

        # dx = Bp * R * dr  -- width of the box in psi space
        self.psiwidth = Bp0 * R * dr
        #print("psi width = %e" % self.psiwidth)

        # Integrated shear
        self.sinty = diff(self.zShift, x) / self.psiwidth

        # Extra expressions to add to grid file
        self._extra = {}

        # Calculate metric terms
        self.metric()


##########################
# Shaped tokamak

class ShapedTokamak(object):
    def __init__(self, Rmaj=6.0, rmin=2.0, dr=0.1, kappa=1.0, delta=0.0, b=0.0, ss=0.0, Bt0=1.0, Bp0 = 0.2):
        """
        Rmaj  - Major radius [metric]
        rmin  - Minor radius [metric]
        dr    - Radial width of region [metric]

        kappa - Ellipticity, 1 for a circle
        delta - Triangularity, 0 for circle
        b     - Indentation ("bean" shape), 0 for circle

        ss    - Shafranov shift [metric]

        Bt0   - Toroidal magnetic field on axis [T]. Varies as 1/R
        Bp0   - Poloidal field at outboard midplane [T]

        Outputs
        -------

        Assigns member variables

        x, y    - Symbols for x and y coordinates

        R (x,y)
        Z (x,y)

        """

        # Have we calculated metric components yet?
        self.metric_is_set = False

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
        expansion = (1 - ss/rmin)*cos(y)/(1 - (ss/rmin))

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
        self.zShift = integrate(nu, (y,0,y), conds='none')

        # Safety factor
        self.shiftAngle = self.zShift.subs(y, 2*pi) - self.zShift.subs(y, 0)

        # Integrated shear
        self.I = diff(self.zShift, x)

        # X has a range [0,1], and y [0,2pi]
        self.x = x
        self.y = y

        self.R = Rxy
        self.Z = Zxy

        self.Bt = Btxy
        self.Bp = Bpxy
        self.B = sqrt(Btxy**2 + Bpxy**2)

        self.hthe = hthe
