#!/usr/bin/env python3
#
# Generate manufactured solution and sources for FCI test
#

from math import pi
import warnings

from boututils.boutwarnings import AlwaysWarning
from boutdata.data import BoutOptionsFile
from boutdata.mms import exprToStr, y, Grad_par, Div_par, Metric
from sympy import sin, cos, Expr

warnings.simplefilter("ignore", AlwaysWarning)

# Length of the y domain
Ly = 10.0

# Identity
metric = Metric()

# Define solution in terms of input x,y,z
f = 1 + 0.1 * sin(2 * y)
fv = 1 + 0.1 * cos(3 * y)


# Turn solution into real x and z coordinates
replace = [(y, metric.y * 2 * pi / Ly)]

f = f.subs(replace)
fv = fv.subs(replace)
v = fv / f

# Substitute back to get input y coordinates
replace = [(metric.y, y * Ly / (2 * pi))]


def Grad2_par2(f: Expr) -> Expr:
    return Grad_par(Grad_par(f))


def Div_par_K_Grad_par(f: Expr, K: Expr) -> Expr:
    return (K * Grad2_par2(f)) + (Div_par(K) * Grad_par(f))


############################################
# Equations solved

options = BoutOptionsFile("data/BOUT.inp")

for name, expr in (
    ("input_field", f),
    ("v", v),
    ("FV_Div_par_solution", Div_par(f * v)),
    ("FV_Div_par_K_Grad_par_solution", Div_par_K_Grad_par(f, v)),
    ("FV_Div_par_K_Grad_par_mod_solution", Div_par_K_Grad_par(f, v)),
    ("FV_Div_par_mod_solution", Div_par(f * v)),
    ("FV_Div_par_fvv_solution", Div_par(f * v * v)),
):
    expr_str = exprToStr(expr.subs(replace))
    print(f"{name} = {expr_str}")
    options[name] = expr_str

options.write("data/BOUT.inp", overwrite=True)
