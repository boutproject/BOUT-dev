from __future__ import print_function
from builtins import str

#
# Generate the test case using SymPy
#
# Equations are:
#
# d/dt(n) = Dx * D2DX2(n) + Dy * D2DY2(n) + Dz * D2DZ2(n)
#

from sympy import symbols, sin, diff

####


def bracket(f, g):
    """
    Calculates [f,g] symbolically
    """

    dfdx = diff(f, x)
    dfdz = diff(f, z)

    dgdx = diff(g, x)
    dgdz = diff(g, z)

    return dfdz * dgdx - dfdx * dgdz


def DDZ(f):
    return diff(f, z)


def Delp2(f):
    """Laplacian in X-Z"""
    d2fdx2 = diff(f, x, 2)
    d2fdz2 = diff(f, z, 2)

    return d2fdx2 + d2fdz2


def Delp4(f):
    d4fdx4 = diff(f, x, 4)
    d4fdz4 = diff(f, z, 4)

    return d4fdx4 + d4fdz4


def D2DX2(f):
    return diff(f, x, 2)


def D2DY2(f):
    return diff(f, y, 2)


def D2DZ2(f):
    return diff(f, z, 2)


def exprToStr(expr):
    """Convert a sympy expression to a string for BOUT++ input"""
    return str(expr).replace("**", "^")  # Replace exponent operator


####

# Parameters
Dx = 1.0
Dy = 0
Dz = 0

# Define symbols

x = symbols("x")
y = symbols("y")
z = symbols("z")
t = symbols("t")
pi = symbols("pi")

# Define the manufactured solution

# n = 0.9 + 0.9*x + 0.2*sin(5.*x**2 - 2*z) + cos(y)
n = 0.9 + 0.9 * x + 0.2 * sin(5.0 * x ** 2)


# Calculate gradients for boundaries

dndx = diff(n, x)
dndy = diff(n, y)

# Calculate RHS function

dndt = Dx * D2DX2(n) + Dy * D2DY2(n) + Dz * D2DZ2(n)

# Calculate sources

Sn = diff(n, t) - dndt

#######################

print("[n]")
print("solution = " + exprToStr(n))
print("\nddx     = " + exprToStr(dndx))
print("\nddy     = " + exprToStr(dndy))
print("\nsource   = " + exprToStr(Sn))
