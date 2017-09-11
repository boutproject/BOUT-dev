#
# Generate manufactured solution and sources for FCI test
#

from __future__ import division
from __future__ import print_function

from boutdata.mms import *

from sympy import sin, cos, sqrt

from math import pi

f = sin(y - z) + cos(t)*sin(y - 2*z)

g = cos(y - z) - cos(t)*sin(y - 2*z)

Lx = 0.1
Ly = 10.
Lz = 1.

Bt = 1.0
Bp = 0.05
Bpprime = 0.1

Bpx = Bp + (x-0.5)*Lx * Bpprime  # Note: x in range [0,1]
B = sqrt(Bpx**2 + Bt**2)

def FCI_Grad_par(f):
    return ( Bt * diff(f, y)*2.*pi/Ly + Bpx * diff(f, z)*2.*pi/Lz ) / B

############################################
# Equations solved

dfdt = FCI_Grad_par(g)
dgdt = FCI_Grad_par(f)

# Loop over variables and print solution, source etc.
for v, dvdt, name in [ (f, dfdt, "f"), (g, dgdt, "g") ]:
    # Calculate source
    S = diff(v, t) - dvdt
    
    print("\n["+name+"]")
    print("solution = "+exprToStr(v))
    print("\nsource = "+exprToStr(S))
    print("\nbndry_par_all = parallel_dirichlet("+name+":solution)")
    
