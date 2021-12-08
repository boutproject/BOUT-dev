from __future__ import print_function
from __future__ import division
from past.utils import old_div

# Generate MMS solutions for a tokamak geometry

from boutdata.mms import *

from sympy import sin, cos

shape = SimpleTokamak()
metric = shape.metric()  # Get the metric tensor

###
# Define solution in normalised x,y coordinates
# NOTE: These are orthogonal tokamak coordinates
#       so y is poloidal angle, not parallel coordinate

drive = sin(6 * x ** 2 - z + y)  # Constant drive for advection

advect = cos(4 * x ** 2 + z - y) + sin(t) * sin(3 * x + 2 * z - y)
delp2 = cos(4 * x ** 2 + z - y) + sin(t) * sin(3 * x + 2 * z - y)
laplacepar = cos(4 * x ** 2 + z - y) + sin(t) * sin(3 * x + 2 * z - y)

ZMAX = 1

# Turn solution into real x and z coordinates
# NOTE: Z is shifted, so y is now the parallel coordinate
zShift = shape.zShift.subs(x, old_div(metric.x, metric.Lx))
sinty = shape.sinty.subs(x, old_div(metric.x, metric.Lx))
replace = [(x, old_div(metric.x, metric.Lx)), (z, old_div(metric.z, ZMAX) + zShift)]

drive = drive.subs(replace)
advect = advect.subs(replace)
delp2 = delp2.subs(replace)
laplacepar = laplacepar.subs(replace)

# Calculate time derivatives

dadt = -1e-3 * bracket(drive, advect, metric)

dddt = 1e-5 * Delp2(delp2, metric)
dgdt = Laplace_par(laplacepar, metric)

vars = [
    (advect, dadt, "advect"),
    (delp2, dddt, "delp2"),
    (laplacepar, dgdt, "laplacepar"),
]


# Change back to normalised coordinates
# NOTE: Z remains shifted, so y is still parallel coordinate
#       Need to scale the toroidal angle part of z, not zShift part
replace = [
    (metric.x, x * metric.Lx),
    (metric.z, (z - shape.zShift) * ZMAX + shape.zShift),
]

drive = drive.subs(replace)
print("[drive]")
print("solution = " + exprToStr(drive))

# Loop over variables and print solution, source etc.
for f, dfdt, name in vars:
    # Calculate source
    S = diff(f, t) - dfdt

    # Differentials for boundary conditions
    dfdx = diff(f, metric.x)
    dfdy = diff(f, metric.y)

    # Substitute back to get in terms of x,y,z

    f = f.subs(replace)
    dfdt = dfdt.subs(replace)
    S = S.subs(replace)
    dfdx = dfdx.subs(replace)
    dfdy = dfdy.subs(replace)

    print("\n[" + name + "]")
    print("solution = " + exprToStr(f))
    print("\nddx = " + exprToStr(dfdx))
    print("\nddy = " + exprToStr(dfdy))
    print("\nsource = " + exprToStr(S))

    # print("\nSource magnitude: %e" % exprMag(S))
