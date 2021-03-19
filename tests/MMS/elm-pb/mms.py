# Generate MMS solutions for tokamak geometry

from __future__ import division
from __future__ import print_function

from boutdata.mms import *

from sympy import sin, cos

from math import pi

shape = SimpleTokamak()
metric = shape.metric()  # Get the metric tensor

###
# Define solution in normalised x,y coordinates
# NOTE: These are orthogonal tokamak coordinates
#       so y is poloidal angle, not parallel coordinate

phi = (sin(z - x + t) + 0.001 * cos(y - z)) * sin(
    2.0 * pi * x
)  # Must satisfy Dirichlet BCs for now

Psi = 1e-2 * cos(4 * x ** 2 + z - y)  # + sin(t)*sin(3*x + 2*z - y))

U = 2.0 * cos(2 * t) * cos(x - z + 4 * y)

P = 1 + 0.5 * cos(t) * cos(3 * x ** 2 - z + y) + 0.005 * sin(y - z) * sin(t)

P0 = 2 + cos(x * pi)  # Pressure pedestal
J0 = 1 - x - sin(x * pi) ** 2 * cos(y)  # Parallel current
bxcvz = -((1.0 / shape.Rxy) ** 2) * cos(y)  # Curvature

eta = 1e-1  # core_resist =  1 / core_lund
hyperresist = -1e-6  # negative -> none

viscos_par = 1.0

ZMAX = 1

nonlinear = True
diamag = True

# Turn solution into real x and z coordinates
# NOTE: Z is shifted, so y is now the parallel coordinate
zShift = shape.zShift.subs(x, metric.x / metric.Lx)
sinty = shape.sinty.subs(x, metric.x / metric.Lx)
replace = [(x, metric.x / metric.Lx), (z, metric.z / ZMAX + zShift)]

phi = phi.subs(replace)
Psi = Psi.subs(replace)
U = U.subs(replace)
P = P.subs(replace)

P0 = P0.subs(replace)
J0 = J0.subs(replace)

# Calculate time derivatives

B0 = metric.B


def Grad_parP(f):
    result = Grad_par(f, metric)
    if nonlinear:
        result -= B0 * bracket(Psi, f, metric)

    return result


##########################################
# Normalise

MU0 = 4.0e-7 * pi
Mi = 2.0 * 1.6726e-27  # Ion mass [kg]

Bbar = 1  # Tesla
Lbar = 1  # meters
# Note: With this choice, all geometrical normalisation factors go to 1.

Jpar = Delp2(Psi, metric)

##########################################
# Parallel electric field

dPsidt = -Grad_parP(phi) + eta * Jpar

if hyperresist > 0.0:
    dPsidt -= eta * hyperresist * Delp2(Jpar, metric)

# a = Grad_par(B0*phi, metric)/B0
# b = bracket(Psi, B0*phi, metric)
# c = eta * Jpar
# d = eta*hyperresist * Delp2(Jpar, metric)

##########################################
# Vorticity

dUdt = (
    B0 ** 2 * b0xGrad_dot_Grad(Psi, J0, metric)
    + bxcvz * diff(P, metric.z)
    - (B0 ** 2) * Grad_parP(Jpar)
)
if nonlinear:
    # Bracket method '0' (STD) goes to b0xGrad_dot_Grad
    dUdt -= b0xGrad_dot_Grad(phi, U, metric)

# a = B0**2 * b0xGrad_dot_Grad(Psi, J0, metric)
# b = bxcvz*diff(P, metric.z)
# c = (B0**2) * Grad_parP(Jpar)

##########################################
# Pressure

dPdt = -b0xGrad_dot_Grad(phi, P0, metric)

if nonlinear:
    # Bracket method '0' (STD) goes to b0xGrad_dot_Grad
    dPdt -= b0xGrad_dot_Grad(phi, P, metric)

vars = [(Psi, dPsidt, "Psi"), (U, dUdt, "U"), (P, dPdt, "P")]

# Change back to normalised coordinates
# NOTE: Z remains shifted, so y is still parallel coordinate
#       Need to scale the toroidal angle part of z, not zShift part
replace = [
    (metric.x, x * metric.Lx),
    (metric.z, (z - shape.zShift) * ZMAX + shape.zShift),
]

# For applying boundary conditions to shifted fields, remove zShift
replace_shiftbc = [(metric.x, x * metric.Lx), (metric.z, (z - shape.zShift) * ZMAX)]

# print "MAG: ", exprMag(a.subs(replace)), exprMag(b.subs(replace)), exprMag(c.subs(replace))#, exprMag(d.subs(replace))


# Potential
if diamag:
    # Delp2(phi + 0.5*P/B0) = U + Sphi
    Sphi = Delp2(phi + 0.5 * P / B0, metric) - U
else:
    # Delp2(phi) = U + Sphi
    Sphi = Delp2(phi, metric) - U
phi = phi.subs(replace)
Sphi = Sphi.subs(replace)
print("[phi]")
print("solution = " + exprToStr(phi))
print("\nsource = " + exprToStr(Sphi))

Jpar = Jpar.subs(replace)
print("\n[J]")
print("solution = " + exprToStr(Jpar))


# Loop over variables and print solution, source etc.
for f, dfdt, name in vars:
    # Calculate source
    S = diff(f, t) - dfdt

    # Differentials for boundary conditions
    dfdx = diff(f, metric.x)
    dfdy = diff(f, metric.y)

    # Substitute back to get in terms of x,y,z

    fbc = f.subs(replace_shiftbc)
    f = f.subs(replace)
    dfdt = dfdt.subs(replace)
    S = S.subs(replace)
    dfdx = dfdx.subs(replace)
    dfdy = dfdy.subs(replace)

    print("\n[" + name + "]")
    print("solution = " + exprToStr(f))
    print("solution_zshift = " + exprToStr(fbc))
    print("\nddx = " + exprToStr(dfdx))
    print("\nddy = " + exprToStr(dfdy))
    print("\nsource = " + exprToStr(S))

    # print("\nSource magnitude: %e" % exprMag(S))
