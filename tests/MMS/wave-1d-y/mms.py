from __future__ import print_function

#
# Generate the test case using SymPy
#
# Equations are:
#
# df/dt = DDY(g)
# dg/dt = DDY(f)
#
#

from sympy import cos, sin, diff

from boutdata.mms import exprToStr, y, t

# Define the manufactured solution

f = y + cos(y) - sin(t) * cos(0.5 * y)
g = y ** 2 + sin(y) + cos(t) * cos(0.1 * y * y)

# Calculate gradients in x for boundaries

dfdy = diff(f, y)
dgdy = diff(g, y)

# Calculate sources

Sf = diff(f, t) - dgdy
Sg = diff(g, t) - dfdy

#######################

print("\n[f]")
print("solution = " + exprToStr(f))
print("ddy      = " + exprToStr(dfdy))
print("source   = " + exprToStr(Sf))

print("\n[g]")
print("solution = " + exprToStr(g))
print("ddy      = " + exprToStr(dgdy))
print("source   = " + exprToStr(Sg))
