from __future__ import print_function
from __future__ import division
from past.utils import old_div

#
# Generate the test case using SymPy
#
# Equations are:
#
#
# Delp2(phi) = vort + Sphi

from boutdata.mms import *

from sympy import log, Wild

####


def C(f):
    """Curvature"""
    return bxcvz * diff(f, metric.z)


#######################################
# Inputs, essentially as in BOUT.inp

# Switches
ionvis = False  # Ion Viscosity
Ti = 10  # Ion temperature for viscosity calculation
elecvis = False  # Electron viscosity
resistivity = True

parallel = True  # Parallel dynamics

# Parameters
Tnorm = 3  # Electron Temperature (eV)
Nnorm = 1e19  # Background plasma density (m^-3)
Bnorm = 0.1  # Magnetic field [T]
AA = 0.1  # Ion atomic mass

#######################################
#

mi_me = AA * 1.67262158e-27 / 9.109e-31
beta_e = qe * Tnorm * Nnorm / (old_div(Bnorm ** 2, mu0))

# Normalisation parameters
Cs0 = sqrt(qe * Tnorm / (AA * Mp))
Omega_ci = qe * Bnorm / (AA * Mp)
rho_s0 = old_div(Cs0, Omega_ci)
Coulomb = 6.6 - 0.5 * log(Nnorm * 1e-20) + 1.5 * log(Tnorm)
tau_e0 = old_div(
    1.0, (2.91e-6 * (old_div(Nnorm, 1e6)) * Coulomb * Tnorm ** (old_div(-3.0, 2)))
)
# Input settings

ZMAX = 1e-3

MXG = 2
MYG = 2

bxcvz = 100  # Curvature

Rxy = 1.5  # Major radius
Bpxy = 0.35  # Poloidal magnetic field
Bxy = 0.35  # Total magnetic field
Btxy = 0.0  # Toroidal magnetic field
hthe = 0.1  # Poloidal arc length

sinty = 0.0
sbp = 1.0

dx = 2e-5
dy = 1e-3

nx = 132
ny = 64

estatic = True

# Normalise

Rxy /= rho_s0
Bpxy /= Bnorm
Btxy /= Bnorm
Bxy /= Bnorm
hthe /= rho_s0

dx /= rho_s0 ** 2 * Bnorm

bxcvz *= rho_s0 ** 2

# Define a metric
metric = Metric()  # Identity

Lx = dx * (nx - 2.0 * MXG)  # Size of the X domain
Ly = dy * ny  # Size of the Y domain

metric.g11 = (Rxy * Bpxy) ** 2
metric.g22 = old_div(1.0, (hthe ** 2))
metric.g33 = (sinty ** 2) * metric.g11 + old_div((Bxy ** 2), metric.g11)
metric.g12 = 0.0
metric.g13 = -sinty * metric.g11
metric.g23 = -sbp * Btxy / (hthe * Bpxy * Rxy)

metric.J = old_div(hthe, Bpxy)
B = metric.B = Bxy

metric.g_11 = old_div(1.0, metric.g11) + ((sinty * Rxy) ** 2)
metric.g_22 = (Bxy * hthe / Bpxy) ** 2
metric.g_33 = Rxy * Rxy
metric.g_12 = sbp * Btxy * hthe * sinty * Rxy / Bpxy
metric.g_13 = sinty * Rxy * Rxy
metric.g_23 = sbp * Btxy * hthe * Rxy / Bpxy


print("dx = %e, dy = %e" % (dx, dy))
print("g11 = %e, g22 = %e, g33 = %e" % (metric.g11, metric.g22, metric.g33))
print("g12 = %e, g23 = %e" % (metric.g12, metric.g23))
print("g_11 = %e, g_22 = %e, g_33 = %e" % (metric.g_11, metric.g_22, metric.g_33))
print("g_12 = %e, g_23 = %e" % (metric.g_12, metric.g_23))

# Define the manufactured solution in terms of input x,y and z

Ne = 0.9 + 0.9 * x + 0.5 * cos(t) * sin(5.0 * x ** 2 - z) + 0.01 * sin(y - z)
Te = 1 + 0.5 * cos(t) * cos(3 * x ** 2 - 2 * z) + 0.005 * sin(y - z) * sin(t)
Vort = 2.0 * sin(2 * t) * cos(x - z + 4 * y)
VePsi = cos(1.5 * t) * (2.0 * sin((x - 0.5) ** 2 + z) + 0.05 * cos(3 * x ** 2 + y - z))
Vi = -0.01 * cos(7 * t) * cos(3 * x ** 2 + 2 * y - 2 * z)

# phi = sin(z - x + t)*sin(2.*pi*x)
phi = (sin(z - x + t) + 0.001 * cos(y - z)) * sin(
    2.0 * pi * x
)  # Must satisfy Dirichlet BCs for now
psi = sin(pi * x) * (
    0.5 * x - 0.1 * cos(7 * t) * sin(3.0 * x ** 2 + y - z)
)  # Must satisfy Dirichlet BCs for now

# Substitute to get in terms of actual x,y,z coordinates

replace = [
    (x, old_div(metric.x, Lx)),
    (y, 2 * pi * metric.y / Ly),
    (z, old_div(metric.z, ZMAX)),
]

Ne = Ne.subs(replace)
Te = Te.subs(replace)
Vort = Vort.subs(replace)
VePsi = VePsi.subs(replace)
Vi = Vi.subs(replace)
phi = phi.subs(replace)
psi = psi.subs(replace)

# Calculate RHS function

Pe = Ne * Te
if estatic:
    Ve = VePsi
else:
    Ve = VePsi - 0.5 * mi_me * beta_e * psi

print("mi_me = ", mi_me)
print("beta_e = ", beta_e)
print("mi_me*beta_e = ", mi_me * beta_e)

Gi = 0.0
Ge = 0.0

tau_e = Omega_ci * tau_e0 * (Te ** 1.5) / Ne
# Normalised collision time

nu = old_div(1.0, (1.96 * Ne * tau_e * mi_me))

### Electron density

dNedt = -bracket(phi, Ne, metric) + (old_div(2, B)) * (C(Pe) - Ne * C(phi))

if parallel:
    dNedt -= Ne * Grad_par(Ve, metric) + Vpar_Grad_par(Ve, Ne, metric)

### Electron temperature

dTedt = -bracket(phi, Te, metric) + (old_div(4.0, 3)) * (old_div(Te, B)) * (
    (old_div(7.0, 2)) * C(Te) + (old_div(Te, Ne)) * C(Ne) - C(phi)
)

if parallel:
    dTedt -= Vpar_Grad_par(Ve, Te, metric)
    dTedt += (
        (old_div(2.0, 3.0))
        * Te
        * (
            0.71 * Grad_par(Vi, metric)
            - 1.71 * Grad_par(Ve, metric)
            + 0.71 * (Vi - Ve) * Grad_par(Ne, metric) / Ne
        )
    )

### Vorticity

dVortdt = -bracket(phi, Vort, metric) + 2.0 * B * C(Pe) / Ne + B * C(Gi) / (3.0 * Ne)

if parallel:
    dVortdt += -Vpar_Grad_par(Vi, Vort, metric) + B ** 2 * (
        Grad_par(Vi - Ve, metric) + (Vi - Ve) * Grad_par(Ne, metric) / Ne
    )

### Parallel Ohm's law

if parallel:
    dVePsidt = (
        -bracket(phi, Ve, metric)
        - Vpar_Grad_par(Ve, Ve, metric)
        - mi_me * (old_div(2.0, 3.0)) * Grad_par(Ge, metric)
        - mi_me * nu * (Ve - Vi)
        + mi_me * Grad_par(phi, metric)
        - mi_me * (Te * Grad_par(Ne, metric) / Ne + 1.71 * Grad_par(Te, metric))
    )
else:
    dVePsidt = 0.0

### Parallel ion velocity

if parallel:
    dVidt = (
        -bracket(phi, Vi, metric)
        - Vpar_Grad_par(Vi, Vi, metric)
        - (old_div(2.0, 3.0)) * Grad_par(Gi, metric)
        - (Grad_par(Te, metric) + Te * Grad_par(Ne, metric) / Ne)
    )
else:
    dVidt = 0.0

#######################

# Create list of all evolving variables
vars = [
    (Ne, dNedt, "Ne"),
    (Te, dTedt, "Te"),
    (Vort, dVortdt, "Vort"),
    (VePsi, dVePsidt, "VePsi"),
    (Vi, dVidt, "Vi"),
]

replace = [(metric.x, x * Lx), (metric.y, y * Ly / (2 * pi)), (metric.z, z * ZMAX)]

# Loop over variables and print solution, source etc.
for v, dvdt, name in vars:
    # Calculate source
    S = diff(v, t) - dvdt

    dvdx = diff(v, metric.x)
    dvdy = diff(v, metric.y)

    # Substitute back to get in terms of x,y,z
    v = v.subs(replace)
    dvdx = dvdx.subs(replace)
    dvdy = dvdy.subs(replace)
    S = S.subs(replace)

    print("[" + name + "]")
    print("solution = " + exprToStr(trySimplify(v)))
    print("\nddx     = " + exprToStr(trySimplify(dvdx)))
    print("\nddy     = " + exprToStr(trySimplify(dvdy)))
    print("\nsource   = " + exprToStr(S))

# Potential
Sphi = Delp2(phi, metric) - Vort

print("\n[phi]")
print("\nsolution = " + exprToStr(phi.subs(replace)))
print("\nsource   = " + exprToStr(Sphi.subs(replace)))


# print "\n\nDelp2 phi = ", Delp2(phi, metric).subs(replace)

if not estatic:

    Spsi = Delp2(psi, metric) - 0.5 * Ne * mi_me * beta_e * psi - Ne * (Vi - VePsi)
    print("\n[psi]")
    print("\nsolution = " + exprToStr(psi.subs(replace)))
    print("\nsource = " + exprToStr(Spsi.subs(replace)))


##########################################
# Check magnitudes of individual terms


def exprmag(expr):
    """
    Estimate the magnitude of an expression

    """
    # First put into the BOUT.inp x,y,z
    expr = expr.subs(replace)

    # Replace all sin, cos with 1
    any = Wild("a")  # Wildcard
    expr = expr.replace(sin(any), 1.0)
    expr = expr.replace(cos(any), 1.0)

    # Pick maximum values of x,y,z
    expr = expr.subs(x, 1.0)
    expr = expr.subs(y, 2.0 * pi)
    expr = expr.subs(z, 2.0 * pi)

    return expr.evalf()


print("\n\n########################################")
print("\n\nDensity terms:")
print("  bracket : ", exprmag(-bracket(phi, Ne, metric)), "\n")
print("  (2/B)*C(Pe) : ", exprmag((old_div(2, B)) * C(Pe)), "\n")
print("  (2/B)*Ne*C(phi) : ", exprmag((old_div(2, B)) * Ne * C(phi)), "\n")
print("  Ne*Grad_par(Ve) :", exprmag(Ne * Grad_par(Ve, metric)), "\n")
print("  Vpar_Grad_par(Ve, Ne)", exprmag(Vpar_Grad_par(Ve, Ne, metric)))

print("\n\n########################################")
print("\n\nTemperature terms:")
print("  bracket : ", exprmag(-bracket(phi, Te, metric)), "\n")
print(
    "  C(Te)   : ",
    exprmag((old_div(4.0, 3)) * (old_div(Te, B)) * (old_div(7.0, 2)) * C(Te)),
    "\n",
)
print(
    "  (Te/Ne)*C(Ne) : ",
    exprmag((old_div(4.0, 3)) * (old_div(Te, B)) * (old_div(Te, Ne)) * C(Ne)),
    "\n",
)
print("  C(phi)  : ", exprmag(-(old_div(4.0, 3)) * (old_div(Te, B)) * C(phi)), "\n")
print("  Vpar_Grad_par(Ve, Te)", exprmag(Vpar_Grad_par(Ve, Te, metric)), "\n")
print(
    "  (2./3.)*Te*( 0.71*Grad_par(Vi)",
    exprmag((old_div(2.0, 3.0)) * Te * 0.71 * Grad_par(Vi, metric)),
    "\n",
)
print(
    "  (2./3.)*Te*1.71*Grad_par(Ve)",
    exprmag(-(old_div(2.0, 3.0)) * Te * 1.71 * Grad_par(Ve, metric)),
    "\n",
)
print(
    "  (2./3.)*Te*0.71*(Vi-Ve)*Grad_par(log(Ne))",
    exprmag((old_div(2.0, 3.0)) * Te * 0.71 * (Vi - Ve) * Grad_par(log(Ne), metric)),
    "\n",
)


print("\n\n########################################")
print("\n\nVorticity terms:")
print("  bracket : ", exprmag(-bracket(phi, Vort, metric)), "\n")
print("  C(Pe)   : ", exprmag(2.0 * B * C(Pe) / Ne), "\n")
print("  C(Gi)   : ", exprmag(B * C(Gi) / (3.0 * Ne)), "\n")
print("  Vpar_Grad_par(Vi, Vort) : ", exprmag(-Vpar_Grad_par(Vi, Vort, metric)), "\n")
print("  B**2*Grad_par(Vi - Ve) :", exprmag(B ** 2 * Grad_par(Vi - Ve, metric)), "\n")
print(
    "  B**2*(Vi - Ve)*Grad_par(log(Ne)) :",
    exprmag(B ** 2 * (Vi - Ve) * Grad_par(log(Ne), metric)),
    "\n",
)

print("\n\n########################################")
print("\n\nOhm's law terms:")

print("  bracket : ", exprmag(-bracket(phi, Ve, metric)), "\n")
print("  Vpar_Grad_par(Ve, Ve) : ", exprmag(-Vpar_Grad_par(Ve, Ve, metric)), "\n")
print(
    "  mi_me*(2./3.)*Grad_par(Ge)",
    exprmag(-mi_me * (old_div(2.0, 3.0)) * Grad_par(Ge, metric)),
    "\n",
)
print("  mi_me*nu*(Ve - Vi)", exprmag(-mi_me * nu * (Ve - Vi)), "\n")
print("  mi_me*Grad_par(phi)", exprmag(mi_me * Grad_par(phi, metric)), "\n")
print(
    "  mi_me*Te*Grad_par(log(Ne))",
    exprmag(-mi_me * Te * Grad_par(log(Ne), metric)),
    "\n",
)
print(
    "  mi_me*1.71*Grad_par(Te) : ", exprmag(-mi_me * 1.71 * Grad_par(Te, metric)), "\n"
)


print("\nVe mag: ", exprmag(Ve), exprmag(VePsi), exprmag(0.5 * mi_me * beta_e * psi))

print("\n\n########################################")
print("\n\nVi terms:")

print("  bracket : ", exprmag(-bracket(phi, Vi, metric)), "\n")
print("  Vpar_Grad_par(Vi, Vi)", exprmag(-Vpar_Grad_par(Vi, Vi, metric)), "\n")
print(
    "  (2./3.)*Grad_par(Gi)", exprmag((old_div(2.0, 3.0)) * Grad_par(Gi, metric)), "\n"
)
print("  Grad_par(Te", exprmag(-Grad_par(Te, metric)), "\n")
print("  Te*Grad_par(log(Ne))", exprmag(-Te * Grad_par(log(Ne), metric)), "\n")
