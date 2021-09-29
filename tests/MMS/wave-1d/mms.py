from __future__ import print_function
from builtins import str

#
# Generate the test case using SymPy
#
# Equations are:
#
# df/dt = DDX(g)
# dg/dt = DDX(f)
#
#

from sympy import symbols, cos, sin, diff

# Define symbols

x = symbols("x")
t = symbols("t")

# Define the manufactured solution

f = 0.9 + 0.9 * x + 0.2 * cos(10 * t) * sin(5.0 * x ** 2)
g = 0.9 + 0.7 * x + 0.2 * cos(7 * t) * sin(2.0 * x ** 2)

# Calculate gradients in x for boundaries

dfdx = diff(f, x)
dgdx = diff(g, x)

# Calculate sources

Sf = diff(f, t) - dgdx
Sg = diff(g, t) - dfdx

#######################

print("F:")
print("solution = " + str(f))
print("d/dx     = " + str(dfdx))
print("source   = " + str(Sf))

print("\nG:")
print("solution = " + str(g))
print("d/dx     = " + str(dgdx))
print("source   = " + str(Sg))
