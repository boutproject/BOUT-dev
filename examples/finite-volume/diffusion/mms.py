from __future__ import division
from __future__ import print_function

from boutdata.mms import Metric, sin, Div_par, Grad_par, exprToStr, diff, y, t
from math import pi

# Length of the y domain
Ly = 10.

# metric tensor
metric = Metric()  # Identity

# Define solution in terms of input x,y,z

f = 1 + 0.1*sin(2*y - t)
k = 1 + 0.1*sin(y)

# Turn solution into real x and z coordinates
replace = [ (y, metric.y*2*pi/Ly) ]

f = f.subs(replace)
k = k.subs(replace) 

##############################
# Calculate time derivatives

dfdt = Div_par( k * Grad_par(f) )

#############################
# Calculate sources

Sf = diff(f, t) - dfdt

# Substitute back to get input y coordinates
replace = [ (metric.y, y*Ly/(2*pi) ) ]

k = k.subs(replace)
f = f.subs(replace)

Sf = Sf.subs(replace)

print("[mesh]")
print("k = " + exprToStr(k))

print("\n[f]")
print("solution = " + exprToStr(f))
print("\nsource = " + exprToStr(Sf))
