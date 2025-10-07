#!/usr/bin/env python3
#
# Generate manufactured solution and sources for FCI test
#

from math import pi
import warnings

from boututils.boutwarnings import AlwaysWarning
from boutdata.data import BoutOptionsFile
from boutdata.mms import diff, exprToStr, x, y, z
from sympy import sin, sqrt, Expr

warnings.simplefilter("ignore", AlwaysWarning)

f = sin(y - z) + sin(y - 2 * z)

Lx = 0.1
Ly = 10.0
Lz = 1.0

Bt = 1.0
Bp = 0.05
Bpprime = 0.1

Bpx = Bp + (x - 0.5) * Lx * Bpprime  # Note: x in range [0,1]
B = sqrt(Bpx**2 + Bt**2)


def FCI_ddy(f: Expr) -> Expr:
    return (Bt * diff(f, y) * 2.0 * pi / Ly + Bpx * diff(f, z) * 2.0 * pi / Lz) / B


def FCI_div_par(f: Expr) -> Expr:
    return Bpx * FCI_ddy(f / Bpx)


############################################
# Equations solved

input_field = exprToStr(f)
grad_par_solution = exprToStr(FCI_ddy(f))
grad2_par2_solution = exprToStr(FCI_ddy(FCI_ddy(f)))
div_par_solution = exprToStr(FCI_div_par(f))

print(f"input_field = {input_field}")
print(f"grad_par_solution = {grad_par_solution}")
print(f"grad2_par2_solution = {grad2_par2_solution}")
print(f"div_par_solution = {div_par_solution}")
print(f"div_par_K_grad_par_solution = {div_par_K_grad_par_solution}")

options = BoutOptionsFile("data/BOUT.inp")
options["input_field"] = input_field
options["grad_par_solution"] = grad_par_solution
options["grad2_par2_solution"] = grad2_par2_solution
options["div_par_solution"] = div_par_solution
options.write("data/BOUT.inp", overwrite=True)
