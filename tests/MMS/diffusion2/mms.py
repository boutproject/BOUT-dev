from __future__ import print_function
from builtins import str

#
# Generate the test case using SymPy
#
# Equations are:
#
# d/dt(n) = Dx * D2DX2(n) + Dy * D2DY2(n) + Dz * D2DZ2(n)
#

from sympy import symbols, cos, sin, diff

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
    tmp = str(expr).replace("**", "^")  # Replace exponent operator
    tmp = tmp.replace("xl", "N:xl")
    return tmp.replace("yl", "N:yl")


####

# Parameters
Dx = 1.0
Dy = 1.0
Dz = 1.0

# Define symbols

x = symbols("x")
y = symbols("y")
xl = symbols("xl")
yl = symbols("yl")
z = symbols("z")
t = symbols("t")
pi = symbols("pi")

# Define the manufactured solution

n = (
    0.9
    + 0.9 * x
    + 0.2 * sin(5.0 * x ** 2 - 2 * z)
    + cos(10 * z) * cos(y) * sin(y * 7 + 1.234)
)


# Calculate gradients for boundaries

dndx = diff(n, x)
dndy = diff(n, y)

# Calculate RHS function

dndt = Dx * D2DX2(n) + Dy * D2DY2(n) + Dz * D2DZ2(n)

# Calculate sources

Sn = diff(n, t) - dndt

# x and y-domains are Lx and Ly long, respectively
# change to corresponding coordinates. Standard BOUT++ coordinates
# have Lx = 1, Ly = 2pi.
# Scaling takes place in BOUT.inp
scale_coordinates = [(x, xl), (y, yl)]
n = n.subs(scale_coordinates)
dndx = dndx.subs(scale_coordinates)
dndy = dndy.subs(scale_coordinates)
Sn = Sn.subs(scale_coordinates)
#######################

print("[n]")
print("solution = " + exprToStr(n))
print("\nddx     = " + exprToStr(dndx))
print("\nddy     = " + exprToStr(dndy))
print("\nsource   = " + exprToStr(Sn))
