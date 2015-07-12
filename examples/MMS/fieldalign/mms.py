from __future__ import print_function

# Test
from sympy import symbols, cos, sin, diff
from numpy import pi
from boutdata.mms import exprToStr

# Define symbols for orthogonal (theta, phi) coordinates

theta = symbols("theta")
phi   = symbols("phi")
t     = symbols("t")

# Define solution (theta, phi)

f = cos(theta + phi - t)**2

# Calculate solution in (theta, phi)

dfdt = diff(f, theta) + diff(f, phi)

# Source

Sf = diff(f, t) - dfdt

# Define new coordinate system
# 
# y = theta
# z = phi - zShift
#
# where 
# 
# zShift = int_{theta_0}^{theta} nu(theta) dtheta
#
# and nu(theta) is the field-line pitch
#

zShift = 0.0*theta #0.1*theta
nu = diff(zShift, theta)

shiftangle = zShift.subs( [(theta, 2*pi)] ) - zShift.subs( [(theta, 0.0)] )

y = symbols("y")
z = symbols("z")

replace = [ (phi, z + zShift),    # z = phi - zShift
            (theta, y) ]          # y = theta

# Replace and print sources, solutions etc.

print("[mesh]")
print("zShift = "+exprToStr( zShift.subs(replace) ))
print("nu = "+exprToStr( nu.subs(replace) ))
print("shiftangle = "+exprToStr(shiftangle))

print("\n[f]")
print("solution = "+exprToStr( f.subs(replace) ))
print("\nsource = "+exprToStr( Sf.subs(replace) ))

