#!/usr/bin/env python3
#
# Generate manufactured solution and sources for FCI test
#

from __future__ import division
from __future__ import print_function

from boutdata.mms import *

from sympy import sin, cos, sqrt

from math import pi

f = sin(y - z) + sin(y - 2 * z)

Lx = 0.1
Ly = 10.0
Lz = 1.0

Bt = 1.0
Bp = 0.05
Bpprime = 0.1

Bpx = Bp + (x - 0.5) * Lx * Bpprime  # Note: x in range [0,1]
B = sqrt(Bpx ** 2 + Bt ** 2)


def FCI_ddy(f):
    return (Bt * diff(f, y) * 2.0 * pi / Ly + Bpx * diff(f, z) * 2.0 * pi / Lz) / B


############################################
# Equations solved

print("input = " + exprToStr(f))
print("solution = " + exprToStr(FCI_ddy(f)))
