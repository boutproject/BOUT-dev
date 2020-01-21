""" Functions for calculating sources for the
    Method of Manufactured Solutions (MMS)

"""

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
    # psiwidth = dx/dx_in
    return diff(f, metric.x)/metric.psiwidth

def DDY(f):
    return diff(f, metric.y)

def DDZ(f):
    return diff(f, metric.z)*metric.zperiod


def D2DX2(f):
    return diff(f, metric.x, 2)/metric.psiwidth**2

def D2DY2(f):
    return diff(f, metric.y, 2)

def D2DZ2(f):
    return diff(f, metric.z, 2)*metric.zperiod**2


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
    result  = metric.G1*DDX(f) + metric.G2*DDY(f) + metric.G3*DDZ(f)\
	      + metric.g11*D2DX2(f) + metric.g22*D2DY2(f) + metric.g33*D2DZ2(f)\
	      + 2.0*(metric.g12*D2DXDY(f) + metric.g13*D2DXDZ(f) + metric.g23*D2DYDZ(f))

    return result

def Laplace_par(f):
    """
    Div( b (b.Grad(f) ) ) = (1/J) d/dy ( J/g_22 * df/dy )
    """
    return DDY(metric.J/metric.g_22)*DDY(f)/ metric.J

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

        metric.zShift = self.zShift

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

        print("dx = "+exprToStr(metric.psiwidth)+"/(nx-2*mxg)")
        print("dy = 2.*pi/ny")
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

##################################

class SimpleTokamak(BaseTokamak):
    """
    Simple tokamak

    NOTE: This is NOT an equilibrium calculation. The input
    is intended solely for testing with MMS
    """
    def __init__(self, R = 2, Bt = 1.0, eps = 0.1, dr=0.02, psiN0=0.5, zperiod=1, r0=1., q = lambda x:2+x**2):
        """
        R    - Major radius [metric]

        Bt   - Toroidal field [T]

        eps  - Inverse aspect ratio

        dr   - Width of the radial region [metric]

        q(x) - A function which returns the safety factor
               as a function of x in range [0,1]
        
        psiN0- Normalized poloidal flux of inner edge of grid

        r0   - Length-scale to normalize metric components (default 1m)

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

        # Minor radius
        self.r = R * eps

        # Approximate poloidal field for radial width calculation
        Bp0 = Bt * self.r*psiN0 / (q(psiN0) * self.R)

        # dpsi = Bp * R * dr  -- width of the box in psi space
        self.psiwidth = Bp0 * self.R * self.dr
        self.psi0 = Bp0 * R * self.r # value of psi at 'separatrix' taken to be at r, psi=0 at magnetic axis

        # Get safety factor
        self.q = q((x + psiN0*self.psi0)/self.psi0)

        # Toroidal angle of a field-line as function
        # of poloidal angle y
        self.zShift = self.q*(self.y-pi + eps * sin(y-pi))

        # Field-line pitch
        self.nu = self.q*(1 + eps*cos(y-pi)) #diff(self.zShift, y)

        # Coordinates of grid points
        self.Rxy = R - self.r*psiN0 * cos(y-pi)
        self.Zxy = self.r*psiN0 * sin(y-pi)

        # Poloidal arc length
        self.hthe = self.r*psiN0 + 0.*y

        # Toroidal magnetic field
        self.Btxy = Bt * R / self.Rxy

        # Poloidal magnetic field
        self.Bpxy = self.Btxy * self.hthe / (self.nu * self.Rxy)

        # Total magnetic field
        self.Bxy = sqrt(self.Btxy**2 + self.Bpxy**2)

        # Integrated shear
        self.sinty = diff(self.zShift, x)

        # Convert all "x" symbols from flux to [0,1]
        xsub = metric.x * self.psiwidth

        self.q = self.q.subs(x, xsub)
        self.zShift = self.zShift.subs(x, xsub)
        self.nu = self.nu.subs(x, xsub)
        self.Rxy = self.Rxy.subs(x, xsub)
        self.Zxy = self.Zxy.subs(x, xsub)
        self.hthe = self.hthe.subs(x, xsub)
        self.Btxy = self.Btxy.subs(x, xsub)
        self.Bpxy = self.Bpxy.subs(x, xsub)
        self.Bxy = self.Bxy.subs(x, xsub)
        self.sinty = self.sinty.subs(x, xsub)

        # Extra expressions to add to grid file
        self._extra = {}

        # Calculate metric terms
        self.metric()


###########################
## Shaped tokamak
###########################
#
#class ShapedTokamak(object):
#    def __init__(self, Rmaj=6.0, rmin=2.0, dr=0.1, kappa=1.0, delta=0.0, b=0.0, ss=0.0, Bt0=1.0, Bp0 = 0.2, zperiod=1):
#        """
#        Rmaj  - Major radius [metric]
#        rmin  - Minor radius [metric]
#        dr    - Radial width of region [metric]
#
#        kappa - Ellipticity, 1 for a circle
#        delta - Triangularity, 0 for circle
#        b     - Indentation ("bean" shape), 0 for circle
#
#        ss    - Shafranov shift [metric]
#
#        Bt0   - Toroidal magnetic field on axis [T]. Varies as 1/R
#        Bp0   - Poloidal field at outboard midplane [T]
#
#        Outputs
#        -------
#
#        Assigns member variables
#
#        x, y    - Symbols for x and y coordinates
#
#        R (x,y)
#        Z (x,y)
#
#        """
#
#        xin = symbols("xin") # xin=x/psiwidth
#        self.zperiod = zperiod
#
#        # Have we calculated metric components yet?
#        self.metric_is_set = False
#
#        # Minor radius as function of xin
#        rminx = rmin + (xin-0.5)*dr
#
#        # Analytical expression for R and Z coordinates as function of x and y
#        self.Rxy = Rmaj - b + (rminx + b*cos(y))*cos(y + delta*sin(y)) + ss*(0.5-xin)*(dr/rmin)
#        self.Zxy = kappa * rminx * sin(y)
#
#        # Toroidal magnetic field
#        self.Btxy = Bt0 * Rmaj / self.Rxy
#
#        # Poloidal field. dx constant, so set poloidal field
#        # at outboard midplane (y = 0)
#        # NOTE: Approximate calculation
#
#        # Distance between flux surface relative to outboard midplane.
#        expansion = (1 - ss/rmin)*cos(y)/(1 - (ss/rmin))
#
#        self.Bpxy = Bp0 * ((Rmaj + rmin) / self.Rxy) / expansion
#        self.B = sqrt(self.Btxy**2 + self.Bpxy**2)
#
#        # Calculate hthe
#        self.hthe = sqrt(diff(self.Rxy, y)**2 + diff(self.Zxy, y)**2)
#        try:
#            self.hthe = trigsimp(self.hthe)
#        except ValueError:
#            pass
#
#        # calculate width in psi
#        drdxin = diff(self.Rxy, xin).subs(y, 0)
#        dpsidr = (self.Bpxy * self.Rxy).subs(y, 0)
#        self.psiwidth = integrate(dpsidr * drdxin, (xin, 0, 1))
#
#        # Field-line pitch
#        nu = self.Btxy * self.hthe / (self.Bpxy * self.Rxy)
#
#        # Shift angle
#        # NOTE: Since x has a range [0,1] this could be done better
#        # than ignoring convergence conditions
#        print(nu)
#        self.zShift = integrate(nu, (y,0,y), conds='none')
#
#        # Safety factor
#        self.shiftAngle = self.zShift.subs(y, 2*pi) - self.zShift.subs(y, 0)
#        self.q = self.shiftAngle/2/pi
#
#        # Integrated shear
#        self.I = diff(self.zShift.subs(xin, x/self.psiwidth), x)
#
#        # X has a range [0,psiwidth], and y [0,2pi]
#        self.x = x
#        self.y = y
#
#        # Convert all "x" symbols from flux to [0,1]
#        # Then convert "xin" (which is already [0,1]) to "x"
#        xsub = metric.x * self.psiwidth
#
#        self.q = self.q.subs(x, xsub).subs(xin, x)
#        self.zShift = self.zShift.subs(x, xsub).subs(xin, x)
#        self.nu = self.nu.subs(x, xsub).subs(xin, x)
#        self.Rxy = self.Rxy.subs(x, xsub).subs(xin, x)
#        self.Zxy = self.Zxy.subs(x, xsub).subs(xin, x)
#        self.hthe = self.hthe.subs(x, xsub).subs(xin, x)
#        self.Btxy = self.Btxy.subs(x, xsub).subs(xin, x)
#        self.Bpxy = self.Bpxy.subs(x, xsub).subs(xin, x)
#        self.Bxy = self.Bxy.subs(x, xsub).subs(xin, x)
#        self.sinty = self.sinty.subs(x, xsub).subs(xin, x)
#
#        # Extra expressions to add to grid file
#        self._extra = {}
#
#        # Calculate metric terms
#        self.metric()
