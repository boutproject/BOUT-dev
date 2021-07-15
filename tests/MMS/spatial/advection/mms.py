from __future__ import division
from __future__ import print_function

from boutdata.mms import *

from math import pi

Lx = 2 * pi
Lz = 2 * pi

ZMAX = Lz / (2 * pi)

metric = Metric()  # Identity metric

# Define solution in terms of input x,y,z

g = sin(6 * x ** 2 - z)  # Constant drive for advection

f = cos(4 * x ** 2 + z)

# Turn solution into real x and z coordinates
replace = [(x, metric.x / Lx), (z, metric.z / ZMAX)]

f = f.subs(replace)
g = g.subs(replace)

# Calculate time derivatives

solution = bracket(g, f, metric)

# Substitute back to get input x and z coordinates
replace = [(metric.x, x * Lx), (metric.z, z * ZMAX)]

g = g.subs(replace)
f = f.subs(replace)
solution = solution.subs(replace)

print("g = " + exprToStr(g))
print("f = " + exprToStr(f))
print("solution = " + exprToStr(solution))
