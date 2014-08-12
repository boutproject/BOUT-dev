
nxlist = [4, 8, 16, 32, 64, 128]


from boututils import DataFile

from os.path import isfile

from boutdata.mms import SimpleTokamak, x, y
from sympy import sin, cos

from math import pi

shape = SimpleTokamak()

## Add equilibrium profiles
MU0 = 4.0e-7*pi

J0 = 1 - x - sin(x*pi)**2 * cos(y)   # Parallel current
P0 = 2 + cos(x*pi)   # Pressure pedestal
bxcvz = -(1./shape.Rxy)**2*cos(y)  # Curvature

# Normalisation
Lbar = 1.
Bbar = 1.
J0 = - J0 * shape.Bxy / (MU0 * Lbar)  # Turn into A/m^2
P0 = P0 * Bbar**2 / (2.0*MU0)  # Pascals

shape.add(P0, "pressure")
shape.add(J0, "Jpar0")
shape.add(bxcvz, "bxcvz")

for nx in nxlist:
    # Generate a new mesh file
    
    filename = "grid%d.nc" % nx
    
    if isfile(filename):
        print("Grid file '%s' already exists" % filename)
    else:
        print("Creating grid file '%s'" % filename)
        f = DataFile(filename, create=True)
        shape.write(nx,nx, f)
        f.close()
