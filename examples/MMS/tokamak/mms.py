
# Generate MMS solutions for a tokamak geometry

from boutdata.mms import *

from sympy import sin, cos

shape = SimpleTokamak()
metric = shape.metric()  # Get the metric tensor

###
# Define solution in normalised x,y coordinates

g = sin(6*x**2 - z)   # Constant drive for advection

f = cos(4*x**2 + z)
#f = cos(4*x**2 + z) + sin(t)*sin(3*x + 2*z)

ZMAX = 0.1

# Turn solution into real x and z coordinates
replace = [ (x, metric.x/metric.Lx), (z, metric.z / ZMAX) ]

f = f.subs(replace)
g = g.subs(replace)

# Calculate time derivatives

dfdt = 1e-5*Delp2(f, metric) #-1e-3*bracket(g, f, metric)

# Calculate source
S = diff(f, t) - dfdt

# Differentials for boundary conditions
dfdx = diff(f, metric.x)

# Change back to normalised coordinates
replace = [ (metric.x, x * metric.Lx), (metric.z, z * ZMAX) ]

print exprMag( (metric.g11*diff(f, metric.x, 2)).subs(replace) )
print exprMag( (metric.g33*diff(f, metric.z, 2)).subs(replace) )

g    = g.subs(replace)
f    = f.subs(replace)
dfdt = dfdt.subs(replace)
S    = S.subs(replace)
dfdx = dfdx.subs(replace)

print("[g]")
print("solution = "+exprToStr(g))

print("\n[f]")
print("solution = "+exprToStr(f))
print("\nddx = "+exprToStr(dfdx))
print("\nsource = "+exprToStr(S))

