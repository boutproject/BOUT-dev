#!/usr/bin/env python3

#
# Generate input mesh for slab shear Alfven wave case
# 
############################################################
# Slab inputs


# Size of the domain

Lx = 0.1     # Radial domain size [m]
Ltor = 10.   # "Toroidal" length [m]
Lpol = 1.    # "Poloidal" length [m]

Bt  = 1.0   # Magnetic field [T]
Bp  = 0.1   # Poloidal field at the middle of the domain [T]
Bpprime = 1.0  # Bp gradient [T/m] Bp(x) = Bp + Bpprime * x 

# Number of grid points in each direction

nx = 34
ny = 64
nz = 64   # Note: Never used here

ymid = 32  # Index of the midplane where X-Z mesh is orthogonal

############################################################

from math import pi

# Effective major radius
R = Ltor / (2.*pi)

# Effective minor radius
hthe = Lpol / (2.*pi)

# Safety factor
q = Bt * hthe / (Bp * R)
print("Safety factor: %e" % q)

# Magnetic shear
s = (-q/Bp)*Bpprime * (hthe/q)
print("Magnetic shear: %e" % s)

# Normalised x coordinate from 0 to 1

from numpy import linspace, ndarray, sqrt, arange, amin, amax
x = linspace(0,1, nx)

# Set poloidal magnetic field

Bpx = Bp + (x-0.5) * Lx * Bpprime

Bpxy = ndarray([nx, ny])
for x in range(nx):
    Bpxy[x,:] = Bpx[x]

Bxy = sqrt(Bpxy**2 + Bt**2)

# Calculate change in poloidal flux
dr = Lx / nx        # Constant mesh spacing in radius
dx = Bpxy * R * dr  # Slightly non-uniform mesh in x, since Bp*R is not constant

# Poloidal angle
dy = 2.*pi / ny

# Calculate zShift and magnetic shear.
# Since metrics vary only in x, integration is trivial

zShift = ndarray([nx,ny])
sinty = ndarray([nx,ny])
TwistShift = ndarray(nx)

for x in range(nx):
    qsafe = Bt * hthe / (Bpxy[x,0] * R)  # safety factor
    zShift[x,:] = qsafe * (arange(ny) - ymid) * dy
    TwistShift[x] = qsafe * 2.*pi  # Angle to use for twist-shift condition
    
    # dq/dx = dq/dr * dr/dx 
    dqdx = (-qsafe * Bpprime / Bpxy[x,0]) * (dr / dx[x,0])
    
    # Integrate to get shear
    sinty[x,:] = dqdx * (arange(ny) - ymid) * dy

print("Safety factor varies between %e and %e" % (amin(TwistShift)/(2.*pi), amax(TwistShift)/(2.*pi)))
############################################################
# Write grid file

from boututils.datafile import DataFile

with DataFile("slab.grd.nc", create=True) as d:
    
    d.write("nx", nx)
    d.write("ny", ny)
    
    d.write("dx", dx)
    d.write("dy", dy)

    d.write("Rxy", R)
    d.write("hthe", hthe)
    d.write("Bpxy", Bpxy)
    d.write("Btxy", Bt)
    d.write("Bxy", Bxy)
    
    d.write("sinty", sinty)
    d.write("zShift", zShift)
    d.write("ShiftAngle", TwistShift)
    
