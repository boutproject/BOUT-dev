#!/usr/bin/env python3
#
# Generate manufactured solution and sources for FCI test
#

from math import pi
import warnings

from boututils.boutwarnings import AlwaysWarning
from boutdata.data import BoutOptionsFile
from boutdata.mms import diff, exprToStr, x, y, z
from sympy import sin, cos, sqrt, Expr

warnings.simplefilter("ignore", AlwaysWarning)

f = sin(y - z) + sin(y - 2 * z)
K = cos(z - y)


Lx = 0.1
Ly = 10.0
Lz = 1.0

Bt = 1.0
Bp = 0.05
Bpprime = 0.1

Bpx = Bp + (x - 0.5) * Lx * Bpprime  # Note: x in range [0,1]
B = sqrt(Bpx**2 + Bt**2)


def FCI_grad_par(f: Expr) -> Expr:
    return (Bt * diff(f, y) * 2.0 * pi / Ly + Bpx * diff(f, z) * 2.0 * pi / Lz) / B


def FCI_grad2_par2(f: Expr) -> Expr:
    return FCI_grad_par(FCI_grad_par(f))


def FCI_div_par(f: Expr) -> Expr:
    return Bpx * FCI_grad_par(f / Bpx)


def FCI_div_par_K_grad_par(f: Expr, K: Expr) -> Expr:
    return (K * FCI_grad2_par2(f)) + (FCI_div_par(K) * FCI_grad_par(f))


def FCI_Laplace_par(f: Expr) -> Expr:
    return FCI_div_par(FCI_grad_par(f))


############################################
# Equations solved

options = BoutOptionsFile("data/BOUT.inp")

for name, expr in (
    ("input_field", f),
    ("K", K),
    ("grad_par_solution", FCI_grad_par(f)),
    ("grad2_par2_solution", FCI_grad2_par2(f)),
    ("div_par_solution", FCI_div_par(f)),
    ("div_par_K_grad_par_solution", FCI_div_par_K_grad_par(f, K)),
    ("laplace_par_solution", FCI_Laplace_par(f)),
    ("FV_div_par_mod_solution", FCI_div_par(f * K)),
    ("FV_div_par_fvv_solution", FCI_div_par(f * K * K)),
):
    expr_str = exprToStr(expr)
    print(f"{name} = {expr_str}")
    options[name] = expr_str

options.write("data/BOUT.inp", overwrite=True)
