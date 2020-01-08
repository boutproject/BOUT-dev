#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division

from builtins import str, range

from numpy import *
from scipy.integrate import quad

from boututils.datafile import DataFile

######################################################

nx = 68 # Number of radial grid points
ny = 32 # Number of poloidal (parallel) grid points

varyBp = False

output = "cyclone_"+str(nx)+"x"+str(ny)+".nc"

######################################################

Ni = 1.   # Ion density in 10^20 m^-3
Ti = 1000 # Temperature in eV (Te = Ti)
Rmaj = 4  # Major radius [meters]
q = 1.4   # Safety factor q = r*Bt/(R*Bp)
s = 0.776 # Magnetic shear s = (r/q) dq/dr
eta_i = 3.114   # Ratio of density to temp. length scales eta = L_n / L_T
epsilon = 0.18  # Inverse aspect ratio epsilon = r / R
Rnorm = 6.92    # Ratio of major radius to L_T  Rnorm  = R / L_T
rho_norm = 0.01 # Normalised ion gyro-radius rho_norm = rho_i / L_T
r_wid = 100     # Radial extent, normalised to gyro-radius r_wid = dr / rho_i

Mi = 2.*1.67262158e-27   # Ion mass [kg]. Deuterium

######################################################

def eps_integral(eps, theta=None):
    if theta == None:
        theta = 2.*pi
    return (quad(lambda t: 1./((1. - eps*cos(t))**2), 0., theta))[0]

rminor = Rmaj * epsilon  # Minor radius [m]
L_T = Rmaj / Rnorm       # Temp. length scale [m]
L_n = eta_i * L_T        # Density length scale [m]
rho_i = rho_norm * L_T   # Ion Larmor radius [m]
Bt0 = sqrt(2.*Ti*Mi / 1.602e-19)/rho_i # Toroidal field from rho_i [T]
Bp = rminor * Bt0 * eps_integral(epsilon)/ (2.*pi * q * Rmaj) # Poloidal field [T]

dr = r_wid * rho_i       # Width of domain [m]

theta = 2.*pi * arange(0,float(ny)) / float(ny)

Rxy = zeros([nx, ny])
Zxy = Rxy.copy()
for i in range(ny):
    Rxy[:,i] = Rmaj - rminor*cos(theta[i])
    Zxy[:,i] = rminor * sin(theta[i])

dy = zeros([nx,ny]) + 2.*pi / float(ny)
hthe = zeros([nx,ny]) + rminor

Btxy = Bt0 * Rmaj / Rxy

print("Toroidal field varies from "+str(Bt0*Rmaj/(Rmaj + rminor)) + \
    " to "+str(Bt0*Rmaj/(Rmaj - rminor)))

# Minor radius offset
drprof = dr*((arange(nx)/float(nx-1)) - 0.5)

# q profile
qprof = q + (s*q/rminor) * drprof

print("q varies from "+str(min(qprof))+" to "+str(max(qprof)))

ShiftAngle = qprof * 2.*pi

Bpxy = zeros([nx,ny])
if varyBp:
    # Vary Bp to get shear
    for y in range(ny):
        Bpxy[:,y] = Bp * q / qprof
    print("Poloidal field varies from "+str(amin(Bpxy))+" to "+str(amax(Bpxy)))
else:
    # Constant Bp, but shift angle varies
    Bpxy += Bp

dx = Bp * (dr/float(nx-1)) * Rxy

Bxy = sqrt(Btxy**2 + Bpxy**2)

zShift = zeros([nx, ny])
qint = eps_integral(epsilon)

for i in range(1,ny):
    zShift[:,i] = ShiftAngle * eps_integral(epsilon, theta=theta[i]) / qint

# Make zShift = 0 on outboard midplane (for plotting mainly)
y0 = int(ny/2)
zs0 = zShift[:,y0]
for i in range(ny):
    zShift[:,i] -= zs0

Ni0 = zeros([nx, ny])
Ti0 = Ni0
for i in range(ny):
    Ni0[:,i] = Ni * exp(-drprof / L_n)
    Ti0[:,i] = Ti * exp(-drprof / L_T)
Te0 = Ti0

pressure = Ni0 * (Ti0 + Te0) * 1.602e-19*1.0e20 # In Pascals

Jpar0 = zeros([nx, ny])

# Shape       : Rxy, Zxy
# Differencing: hthe, dx, dy
# Profiles    : Ni0, Ti0, Te0, pressure, Jpar0
# B field     : Btxy, Bpxy, Bxy
# q profile   : qprof

######################################################
# Curvature

# Bxy is constant in x, so need to supply logB too

logB = zeros([nx, ny])

for x in range(nx):
    for y in range(ny):
        rpos = (float(x)/float(nx-1) - 0.5) * dr
        R = Rmaj - (rminor + rpos)*cos(theta[y])
        Bt = Bt0 * Rmaj / R
        logB[x,y] = log(sqrt(Bt**2 + Bp**2))

######################################################
# Topology: Just in the core

ixseps1 = nx
ixseps2 = nx
jyseps1_1 = -1
jyseps1_2 = int(ny/2)
jyseps2_1 = jyseps1_2
jyseps2_2 = ny-1
ny_inner = jyseps1_2

# Only one region
yup_xsplit   = [nx]
ydown_xsplit = [nx]
yup_xin = [0]
yup_xout = [-1]
ydown_xin = [0]
ydown_xout = [-1]
nrad = [nx]
npol = [ny]

######################################################

print("Writing grid to file "+output)

of = DataFile()
of.open(output, create=True)

of.write("nx", nx)
of.write("ny", ny)

# Topology for original scheme
of.write("ixseps1", ixseps1)
of.write("ixseps2", ixseps2)
of.write("jyseps1_1", jyseps1_1)
of.write("jyseps1_2", jyseps1_2)
of.write("jyseps2_1", jyseps2_1)
of.write("jyseps2_2", jyseps2_2)
of.write("ny_inner", ny_inner)

# Grid spacing
of.write("dx", dx)
of.write("dy", dy)

of.write("ShiftAngle", ShiftAngle)
of.write("zShift", zShift)

of.write("Rxy", Rxy)
of.write("Zxy", Zxy)
of.write("Bpxy", Bpxy)
of.write("Btxy", Btxy)
of.write("Bxy", Bxy)
of.write("hthe", hthe)

# Topology for general configurations
of.write("yup_xsplit", yup_xsplit)
of.write("ydown_xsplit", ydown_xsplit)
of.write("yup_xin", yup_xin)
of.write("ydown_xin", ydown_xin)
of.write("ydown_xout", ydown_xout)
of.write("nrad", nrad)
of.write("npol", npol)

# plasma profiles
of.write("pressure", pressure)
of.write("Jpar0", Jpar0)
of.write("Ni0", Ni0)
of.write("Te0", Te0)
of.write("Ti0", Ti0)
of.write("Ni_x", Ni)
of.write("Te_x", Ti)
of.write("Ti_x", Ti)
of.write("bmag", Bt0)
of.write("rmag", Rmaj)

# Curvature
of.write("logB", logB)

of.close()

print("Done")
