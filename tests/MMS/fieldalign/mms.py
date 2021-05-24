from __future__ import print_function

# Test
from sympy import symbols, cos, diff
from numpy import pi as pin
from boutdata.mms import exprToStr

# Define symbols for orthogonal (theta, phi) coordinates

theta = symbols("theta")
phi = symbols("phi")
psi = symbols("psi")
t = symbols("t")
pi = symbols("pi")

# Define solution (theta, phi)

# f = cos(theta + phi - t)**2
f = cos(psi + theta + phi - t) ** 2

# Calculate solution in (theta, phi)

dfdt = diff(f, theta) + diff(f, phi) + diff(f, psi)
dfdpsi = diff(f, psi)
dfdtheta = diff(f, theta)

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

x = symbols("x")
y = symbols("y")
z = symbols("z")

zShift = 0.1 * theta  # *psi
nu = diff(zShift, theta)
I = diff(zShift, psi)
H = diff(zShift, theta)
shiftangle = zShift.subs([(theta, 2 * pin)]) - zShift.subs([(theta, 0.0)])

yShift = 0.2 * psi  # 0.05*(0.5-psi)*sin(theta) #0.2*psi #
eta = diff(yShift, psi)
G = 1 - diff(yShift, theta)

replace = [(phi, z + zShift), (theta, y + yShift), (psi, x)]  # z = phi - zShift
# replace = [ (phi, z + zShift), (psi, x) ]    # z = phi - zShift

# Replace and print sources, solutions etc.

print("[mesh]")
print("zShift = " + exprToStr(zShift.subs(replace)))
print("nu = " + exprToStr(nu.subs(replace)))
print("shiftangle = " + exprToStr(shiftangle.subs(replace)))
print("I = " + exprToStr(I.subs(replace)))
print("HH = " + exprToStr(H.subs(replace)))
print("")
print("eta = " + exprToStr(eta.subs(replace)))
print("yShift = " + exprToStr(yShift.subs(replace)))
print("G = " + exprToStr(G.subs(replace)))


print("\n[f]")
print("solution = " + exprToStr(f.subs(replace)))
print("\nsource = " + exprToStr(Sf.subs(replace)))
print("dfdpsi = " + exprToStr(dfdpsi.subs(replace)))
print("dfdtheta = " + exprToStr(dfdtheta.subs(replace)))
