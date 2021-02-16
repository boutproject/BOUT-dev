from __future__ import print_function
from builtins import str

#
# Generate the test case using SymPy
#
# Equations are:
#
# d/dt(n) = -[phi,n] + alpha*(phi - n) - kappa*DDZ(phi) + Sn
# d/dt(vort) = -[phi,vort] + alpha*(phi - n)  + Svort
#
# Delp2(phi) = vort + Sphi

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


def exprToStr(expr):
    """Convert a sympy expression to a string for BOUT++ input"""
    return str(expr).replace("**", "^")  # Replace exponent operator


####

# Parameters
alpha = 1.0
kappa = 0.5
Dn = 1.0
Dvort = 1.0

# Define symbols

x = symbols("x")
z = symbols("z")
t = symbols("t")
pi = symbols("pi")

# Define the manufactured solution

n = 0.9 + 0.9 * x + 0.2 * cos(10 * t) * sin(5.0 * x ** 2 - 2 * z)
vort = 0.9 + 0.7 * x + 0.2 * cos(7 * t) * sin(2.0 * x ** 2 - 3 * z)
phi = sin(pi * x) * (
    0.5 * x - cos(7 * t) * sin(3.0 * x ** 2 - 3 * z)
)  # Must satisfy Dirichlet BCs for now

# Calculate gradients in x for boundaries

dndx = diff(n, x)
dvortdx = diff(vort, x)
dphidx = diff(phi, x)

# Calculate RHS function

dndt = -bracket(phi, n) + alpha * (phi - n) - kappa * DDZ(phi) + Dn * Delp2(n)

dvortdt = -bracket(phi, vort) + alpha * (phi - n) + Dvort * Delp2(vort)

# Calculate sources

Sn = diff(n, t) - dndt
Svort = diff(vort, t) - dvortdt

Sphi = Delp2(phi) - vort

#######################

print("[n]")
print("solution = " + exprToStr(n))
print("\nddx     = " + exprToStr(dndx))
print("\nsource   = " + exprToStr(Sn))

print("\n[vort]")
print("solution = " + exprToStr(vort))
print("\nddx     = " + exprToStr(dvortdx))
print("\nsource   = " + exprToStr(Svort))

print("\n[phi]")
print("solution = " + exprToStr(phi))
print("\nddx     = " + exprToStr(dphidx))
print("\nsource   = " + exprToStr(Sphi))
