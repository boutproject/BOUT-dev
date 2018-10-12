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

        # coordinate transformations for 'mesh refinement'
        # expressions for these must average to 1.
        self.scalex = 1.0
        self.scaley = 1.0

metric = Metric()

# Basic differencing
def ddt(f):
    """Time derivative"""
    return diff(f, t)


def DDX(f):
    # psiwidth = dx/dx_in
    return diff(f, metric.x)/metric.psiwidth/metric.scalex

def DDY(f):
    return diff(f, metric.y)/metric.scaley

def DDZ(f):
    return diff(f, metric.z)*metric.zperiod


def D2DX2(f):
    return DDX(DDX(f))

def D2DY2(f):
    return DDY(DDY(f))

def D2DZ2(f):
    return DDZ(DDZ(f))


# don't include derivatives of scalex/scaley, to match BOUT++ implementation
# where D4D*4 are used just for numerical 'hyperdiffusion'
def D4DX4(f):
    return diff(f, metric.x, 4)/metric.psiwidth**4/metric.scalex**4

def D4DY4(f):
    return diff(f, metric.y, 4)/metric.scaley**4

def D4DZ4(f):
    return diff(f, metric.z, 4)*metric.zperiod**4


def D2DXDY(f):
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

    dfdx = DDX(f)
    dfdz = DDZ(f)

    dgdx = DDX(g)
    dgdz = DDZ(g)

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
    d2fdx2 = D2DX2(f)
    d2fdz2 = D2DZ2(f)
    d2fdxdz = D2DXDZ(f)

    result = metric.g11*d2fdx2 + metric.g33*d2fdz2 + 2.*metric.g13*d2fdxdz

    if all_terms:
        result += metric.G1 * DDX(f) + metric.G3 * DDZ(f)

    return result

def Delp4(f):
    d4fdx4 = D2DX2(D2DX2(f))
    d4fdz4 = D2DZ2(D2DZ2(f))

    return d4fdx4 + d4fdz4

def Grad_par(f):
    """The parallel gradient"""
    return DDY(f) / sqrt(metric.g_22)

def Grad2_par2(f):
    """The parallel 2nd derivative"""
    return Grad_par(Grad_par(f))

def Vpar_Grad_par(v, f):
    """Parallel advection operator v*grad_||(f)"""
    return v * Grad_par(f)

def Div_par_K_Grad_par(K, f):
    """Parallel diffusion operator Div(b K b.Grad(f))"""
    return Div_par(K*Grad_par(f))

def Div_par(f):
    '''
    Divergence of magnetic field aligned vector v = \bhat f
    \nabla \cdot (\bhat f) = 1/J \partial_y (f/B)
    = B Grad_par(f/B)
    '''
    return metric.B*Grad_par(f/metric.B)

def Laplace(f):
    """The full Laplace operator"""
    result  = metric.G1*DDX(f) + metric.G2*DDY(f) + metric.G3*DDZ(f)\
	      + metric.g11*D2DX2(f) + metric.g22*D2DY2(f) + metric.g33*D2DZ2(f)\
	      + 2.0*(metric.g12*D2DXDY(f) + metric.g13*D2DXDZ(f) + metric.g23*D2DYDZ(f))

    return result

def Laplace_par(f):
    """
    Div( b (b.Grad(f) ) ) = (1/J) d/dy ( J/g_22 * df/dy )
    """
    return DDY( metric.J/metric.g_22 * DDY(f) )/ metric.J

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
    """
    Convert a sympy expression to a string for BOUT++ input
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

        # If true, set integrated shear to zero in metric (for shifted-metric schemes)
        self.shifted

        # Extra expressions to add to grid file
        self._extra
    """

    psiwidth = None
    r0 = None

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

        dx = self.psiwidth * metric.scalex / nx  + 0.*self.x
        dy = 2.*pi * metric.scaley / ny + 0.*self.x

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
            shiftAngle[i] = self.shiftAngle.evalf(subs={self.x:x})

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

        metric.zShift = self.zShift
        metric.shiftAngle = self.shiftAngle

        metric.J = self.hthe / self.Bpxy
        metric.B = self.Bxy

        metric.psiwidth = self.psiwidth
        metric.zperiod = self.zperiod

        # normalize lengths
        if self.r0 is not None:
            metric.g11 *= self.r0**2
            metric.g22 *= self.r0**2
            metric.g33 *= self.r0**2
            metric.g12 *= self.r0**2
            metric.g13 *= self.r0**2
            metric.g23 *= self.r0**2
            metric.g_11 /= self.r0**2
            metric.g_22 /= self.r0**2
            metric.g_33 /= self.r0**2
            metric.g_12 /= self.r0**2
            metric.g_13 /= self.r0**2
            metric.g_23 /= self.r0**2

        self.metric_is_set = True

        # Christoffel symbols
        metric.G1_11 = 0.5 * metric.g11 * DDX(metric.g_11) + metric.g12 * (DDX(metric.g_12) - 0.5 * DDY(metric.g_11)) + metric.g13 * (DDX(metric.g_13) - 0.5 * DDZ(metric.g_11))
        metric.G1_22 = metric.g11 * (DDY(metric.g_12) - 0.5 * DDX(metric.g_22)) + 0.5 * metric.g12 * DDY(metric.g_22) + metric.g13 * (DDY(metric.g_23) - 0.5 * DDZ(metric.g_22))
        metric.G1_33 = metric.g11 * (DDZ(metric.g_13) - 0.5 * DDX(metric.g_33)) + metric.g12 * (DDZ(metric.g_23) - 0.5 * DDY(metric.g_33)) + 0.5 * metric.g13 * DDZ(metric.g_33)
        metric.G1_12 = 0.5 * metric.g11 * DDY(metric.g_11) + 0.5 * metric.g12 * DDX(metric.g_22) + 0.5 * metric.g13 * (DDY(metric.g_13) + DDX(metric.g_23) - DDZ(metric.g_12))
        metric.G1_13 = 0.5 * metric.g11 * DDZ(metric.g_11) + 0.5 * metric.g12 * (DDZ(metric.g_12) + DDX(metric.g_23) - DDY(metric.g_13)) + 0.5 * metric.g13 * DDX(metric.g_33)
        metric.G1_23 = 0.5 * metric.g11 * (DDZ(metric.g_12) + DDY(metric.g_13) - DDX(metric.g_23)) + 0.5 * metric.g12 * (DDZ(metric.g_22) + DDY(metric.g_23) - DDY(metric.g_23)) + 0.5 * metric.g13 * DDY(metric.g_33)

        metric.G2_11 = 0.5 * metric.g12 * DDX(metric.g_11) + metric.g22 * (DDX(metric.g_12) - 0.5 * DDY(metric.g_11)) + metric.g23 * (DDX(metric.g_13) - 0.5 * DDZ(metric.g_11))
        metric.G2_22 = metric.g12 * (DDY(metric.g_12) - 0.5 * DDX(metric.g_22)) + 0.5 * metric.g22 * DDY(metric.g_22) + metric.g23 * (DDY(metric.g23) - 0.5 * DDZ(metric.g_22))
        metric.G2_33 = metric.g12 * (DDZ(metric.g_13) - 0.5 * DDX(metric.g_33)) + metric.g22 * (DDZ(metric.g_23) - 0.5 * DDY(metric.g_33)) + 0.5 * metric.g23 * DDZ(metric.g_33)
        metric.G2_12 = 0.5 * metric.g12 * DDY(metric.g_11) + 0.5 * metric.g22 * DDX(metric.g_22) + 0.5 * metric.g23 * (DDY(metric.g_13) + DDX(metric.g_23) - DDZ(metric.g_12))
        metric.G2_13 = 0.5 * metric.g12 * (DDZ(metric.g_11) + DDX(metric.g_13) - DDX(metric.g_13)) + 0.5 * metric.g22 * (DDZ(metric.g_12) + DDX(metric.g_23) - DDY(metric.g_13)) + 0.5 * metric.g23 * DDX(metric.g_33)
        metric.G2_23 = 0.5 * metric.g12 * (DDZ(metric.g_12) + DDY(metric.g_13) - DDX(metric.g_23)) + 0.5 * metric.g22 * DDZ(metric.g_22) + 0.5 * metric.g23 * DDY(metric.g_33)

        metric.G3_11 = 0.5 * metric.g13 * DDX(metric.g_11) + metric.g23 * (DDX(metric.g_12) - 0.5 * DDY(metric.g_11)) + metric.g33 * (DDX(metric.g_13) - 0.5 * DDZ(metric.g_11))
        metric.G3_22 = metric.g13 * (DDY(metric.g_12) - 0.5 * DDX(metric.g_22)) + 0.5 * metric.g23 * DDY(metric.g_22) + metric.g33 * (DDY(metric.g_23) - 0.5 * DDZ(metric.g_22))
        metric.G3_33 = metric.g13 * (DDZ(metric.g_13) - 0.5 * DDX(metric.g_33)) + metric.g23 * (DDZ(metric.g_23) - 0.5 * DDY(metric.g_33)) + 0.5 * metric.g33 * DDZ(metric.g_33)
        metric.G3_12 = 0.5 * metric.g13 * DDY(metric.g_11) + 0.5 * metric.g23 * DDX(metric.g_22) + 0.5 * metric.g33 * (DDY(metric.g_13) + DDX(metric.g_23) - DDZ(metric.g_12))
        metric.G3_13 = 0.5 * metric.g13 * DDZ(metric.g_11) + 0.5 * metric.g23 * (DDZ(metric.g_12) + DDX(metric.g_23) - DDY(metric.g_13)) + 0.5 * metric.g33 * DDX(metric.g_33)
        metric.G3_23 = 0.5 * metric.g13 * (DDZ(metric.g_12) + DDY(metric.g_13) - DDX(metric.g_23)) + 0.5 * metric.g23 * DDZ(metric.g_22) + 0.5 * metric.g33 * DDY(metric.g_33)

        metric.G1 = (DDX(metric.J * metric.g11) + DDY(metric.J * metric.g12) + DDZ(metric.J * metric.g13)) / metric.J
        metric.G2 = (DDX(metric.J * metric.g12) + DDY(metric.J * metric.g22) + DDZ(metric.J * metric.g23)) / metric.J
        metric.G3 = (DDX(metric.J * metric.g13) + DDY(metric.J * metric.g23) + DDZ(metric.J * metric.g33)) / metric.J

    def print_mesh(self):
        """
        Prints the metrics to stdout to be copied to a BOUT.inp file
        """
        if not self.metric_is_set:
            raise ValueError("Error: metric has not been calculated yet, so cannot print")

        print("dx = "+exprToStr(metric.psiwidth*metric.scalex)+"/(nx-2*mxg)")
        print("dy = 2.*pi*"+exprToStr(metric.scaley)+"/ny")
        print("dz = 2.*pi/nz/"+exprToStr(metric.zperiod))
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
        print("G1_11 = "+exprToStr(metric.G1_11))
        print("G1_22 = "+exprToStr(metric.G1_22))
        print("G1_33 = "+exprToStr(metric.G1_33))
        print("G1_12 = "+exprToStr(metric.G1_12))
        print("G1_13 = "+exprToStr(metric.G1_13))
        print("G1_23 = "+exprToStr(metric.G1_23))
        print("G2_11 = "+exprToStr(metric.G2_11))
        print("G2_22 = "+exprToStr(metric.G2_22))
        print("G2_33 = "+exprToStr(metric.G2_33))
        print("G2_12 = "+exprToStr(metric.G2_12))
        print("G2_13 = "+exprToStr(metric.G2_13))
        print("G2_23 = "+exprToStr(metric.G2_23))
        print("G3_11 = "+exprToStr(metric.G3_11))
        print("G3_22 = "+exprToStr(metric.G3_22))
        print("G3_33 = "+exprToStr(metric.G3_33))
        print("G3_12 = "+exprToStr(metric.G3_12))
        print("G3_13 = "+exprToStr(metric.G3_13))
        print("G3_23 = "+exprToStr(metric.G3_23))
        print("G1 = "+exprToStr(metric.G1))
        print("G2 = "+exprToStr(metric.G2))
        print("G3 = "+exprToStr(metric.G3))

        print("Lx = "+exprToStr(metric.psiwidth))
        print("Lz = "+exprToStr((2.*pi*self.R/metric.zperiod).evalf()))

    def set_scalex(self, expr):
        # check average of expr is 1.
        average = integrate(expr, (x, 0, 1), (y, 0, 2*pi))/2/pi
        if not average == 1:
            raise ValueError("scalex must average to 1. Got "+str(average))
        else:
            metric.scalex = expr

    def set_scaley(self, expr):
        # check average of expr is 1.
        average = integrate(expr, (x, 0, 1), (y, 0, 2*pi))/2/pi
        if not average == 1:
            raise ValueError("scaley must average to 1. Got "+str(average))
        else:
            metric.scaley = expr

##################################

class SimpleTokamak(BaseTokamak):
    """
    Simple tokamak

    NOTE: This is NOT an equilibrium calculation. The input
    is intended solely for testing with MMS
    """
    def __init__(self, R = 2, Bt = 1.0, eps = 0.1, dr=0.02, psiN0=0.5, zperiod=1, r0=1., q = lambda psiN:2+psiN**2, lower_legs_fraction=0., upper_legs_fraction=0., shifted=False):
        """
        R    - Major radius [metric]

        Bt   - Toroidal field [T]

        eps  - Inverse aspect ratio

        dr   - Width of the radial region [metric]

        q(psiN) - A function which returns the safety factor
               as a function of psiN in range [0,1]

        psiN0- Normalized poloidal flux of inner edge of grid

        r0   - Length-scale to normalize metric components (default 1m)

        lower_legs_fraction - what fraction of the poloidal grid is in the lower divertor legs (per leg, so lower_legs_fraction+upper_legs_fraction<=0.5)

        upper_legs_fraction - what fraction of the poloidal grid is in the upper divertor legs (per leg, so lower_legs_fraction+upper_legs_fraction<=0.5)

        Coordinates:
        x - Radial, [0,psiwidth], x=psi-psi_inner so dx=dpsi but x=0 at the inner edge of the grid
        y - Poloidal, [0,2pi]. Origin is at inboard midplane.


        """

        # Have we calculated metric components yet?
        self.metric_is_set = False

        self.x = x
        self.y = y
        self.zperiod = zperiod
        self.r0 = r0
        self.R = R
        self.dr = dr
        self.psiN0 = psiN0
        self.lower_legs_fraction = lower_legs_fraction
        self.upper_legs_fraction = upper_legs_fraction
        self.shifted = shifted

        # Minor radius
        self.r = R * eps

        # Approximate poloidal field for radial width calculation
        Bp0 = Bt * self.r*self.psiN0 / (q(self.psiN0) * self.R)

        # dpsi = Bp * R * dr  -- width of the box in psi space
        self.psiwidth = Bp0 * self.R * self.dr
        self.psi0 = Bp0 * R * self.r # value of psi at 'separatrix' taken to be at r, psi=0 at magnetic axis

        # Get safety factor
        self.q = q((x + self.psiN0*self.psi0)/self.psi0)

        # Toroidal angle of a field-line as function
        # of poloidal angle y
        self.zShift = self.q*(self.y-pi + eps * sin(y-pi))

        # Toroidal shift of field line at branch cut where poloidal angle goes 2pi->0
        self.shiftAngle = 2*pi*self.q

        # Field-line pitch
        self.nu = self.q*(1 + eps*cos(y-pi)) #diff(self.zShift, y)

        # Coordinates of grid points
        self.Rxy = R - self.r*self.psiN0 * cos(y-pi)
        self.Zxy = self.r*self.psiN0 * sin(y-pi)

        # Poloidal arc length
        self.hthe = self.r*self.psiN0 + 0.*y

        # Toroidal magnetic field
        self.Btxy = Bt * R / self.Rxy

        # Poloidal magnetic field
        self.Bpxy = self.Btxy * self.hthe / (self.nu * self.Rxy)

        # Total magnetic field
        self.Bxy = sqrt(self.Btxy**2 + self.Bpxy**2)

        # Integrated shear
        if shifted:
            self.sinty = 0*y
        else:
            self.sinty = diff(self.zShift, x)

        # Convert all "x" symbols from flux to [0,1]
        xsub = metric.x * self.psiwidth

        self.q = self.q.subs(x, xsub)
        self.zShift = self.zShift.subs(x, xsub)
        self.shiftAngle = self.shiftAngle.subs(x, xsub)
        self.nu = self.nu.subs(x, xsub)
        self.Rxy = self.Rxy.subs(x, xsub)
        self.Zxy = self.Zxy.subs(x, xsub)
        self.hthe = self.hthe.subs(x, xsub)
        self.Btxy = self.Btxy.subs(x, xsub)
        self.Bpxy = self.Bpxy.subs(x, xsub)
        self.Bxy = self.Bxy.subs(x, xsub)
        self.sinty = self.sinty.subs(x, xsub)

        # calculating grid spacing needs grid szie: initialize as None
        self.dx = None
        self.dy = None
        self.dz = None

        # calculating topology parameters needs grid size: initialize as None
        self.ixseps1 = None
        self.ixseps2 = None
        self.jyseps1_1 = None
        self.jyseps1_2 = None
        self.jyseps2_1 = None
        self.jyseps2_2 = None

        # Extra expressions to add to grid file
        self._extra = {}

        # Calculate metric terms
        self.metric()

    def setNs(self, nx=None, ny=None, nz=None, mxg=2):
        """
        Pass in the grid sizes and hence calculate grid spacing and topology parameters
        """
        self.nx = nx
        self.ny = ny
        self.nz = nz

        core_fraction = 1.-2.*self.lower_legs_fraction-2.*self.upper_legs_fraction

        # calculate grid spacings
        self.dx = metric.psiwidth*metric.scalex/self.nx
        self.dy = 2.*pi * metric.scaley / (core_fraction * self.ny)
        self.dz = 2.*pi/self.nz

        # find position of separatrix (uniform grid in x), i.e. x-index for which psiN=1
        self.ixseps1 = int(self.nx * (1. - self.psiN0)*self.psi0/self.psiwidth)-1 + mxg # numbering for ixseps* includes x-guard cells
        self.ixseps2 = self.nx+2*mxg # assume connected double-null if double-null: set ixseps2 to BoutMesh default

        if self.ixseps1 < 0 or self.ixseps1 >= self.nx:
            # no separatrix, so no branch cuts: make everywhere 'core', like BoutMesh defaults
            self.jyseps1_1 = -1
            self.jyseps1_2 = self.nx/2
            self.jyseps2_1 = self.jyseps1_2
            self.jyseps2_2 = self.ny-1
        else:
            self.jyseps1_1 = int(self.lower_legs_fraction*self.ny)-1
            self.jyseps2_2 = int((1.-self.lower_legs_fraction)*self.ny)-1
            self.jyseps1_2 = int(self.ny//2 - self.upper_legs_fraction*self.ny)-1
            self.jyseps2_1 = int(self.ny//2 + self.upper_legs_fraction*self.ny)-1
