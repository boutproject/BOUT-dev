""" Functions for calculating sources for the
    Method of Manufactured Solutions (MMS)

"""

from sympy import symbols, cos, sin, diff, sqrt, pi

# Constants
qe = 1.602e-19
Mp = 1.67262158e-27
mu0 = 4.e-7*3.141592653589793

# Define symbols

x = symbols('x')
y = symbols('y')
z = symbols('z')
t = symbols('t')

class Metric:
    def __init__(self):
        # Create an identity metric
        self.x = symbols('x\'')
        self.y = symbols('y\'')
        self.z = symbols('z\'')

        self.g11 = self.g22 = self.g33 = 1.0
        self.g12 = self.g23 = 0.0

        self.g_11 = self.g_22 = self.g_33 = 1.0
        self.g_12 = self.g_23 = 0.0

        self.J = 1.0
        self.B = 1.0

identity = Metric()

# Basic differencing

def DDZ(f, metric = identity):
    return diff(f, metric.z)

def DDX(f, metric = identity):
    return diff(f, metric.x)

def D2DX2(f, metric = identity):
    return diff(f, metric.x, 2)

def D2DY2(f, metric = identity):
    return diff(f, metric.y, 2)

def D2DZ2(f, metric = identity):
    return diff(f, metric.z, 2)

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

def Delp2(f, metric = identity):
    """ Laplacian in X-Z
    """
    d2fdx2 = diff(f, metric.x, 2)
    d2fdz2 = diff(f, metric.z, 2)

    return metric.g11*d2fdx2 + metric.g33*d2fdz2

def Delp4(f, metric = identity):
    d4fdx4 = diff(f, metric.x, 4)
    d4fdz4 = diff(f, metric.z, 4)

    return d4fdx4 + d4fdz4

def Grad_par(f, metric = identity):
    return diff(f, metric.y) / sqrt(metric.g_22)

def Vpar_Grad_par(v, f, metric = identity):
    return v * Grad_par(f, metric=metric)


# Convert expression to string

def exprToStr(expr):
    """ Convert a sympy expression to a string for BOUT++ input
    """
    return str(expr).replace("**", "^") # Replace exponent operator


