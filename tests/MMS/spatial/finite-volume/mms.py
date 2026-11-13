#!/usr/bin/env python3
#
# Generate manufactured solution and sources for FCI test
#

import warnings
from math import pi

from boutdata.data import BoutOptionsFile
from boutdata.mms import Div_par, Grad_par, Metric, exprToStr, y
from boututils.boutwarnings import AlwaysWarning
from sympy import Expr, Integer, cos, sin

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

# Velocity profile with stagnation points at y = 0 and y = pi
v_stagnation = 0.5 * sin(y)

f = f.subs(replace)
fv = fv.subs(replace)
v = fv / f
v_stagnation = v_stagnation.subs(replace)

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
    ("v_stagnation", v_stagnation),
    # Constant wave speed for the stagnation tests, larger than |v| and |v_stagnation|
    ("wave_speed", Integer(2)),
    ("FV_Div_par_mod_stagnation_solution", Div_par(f * v_stagnation)),
    ("FV_Div_par_fvv_stagnation_solution", Div_par(f * v_stagnation * v_stagnation)),
):
    expr_str = exprToStr(expr.subs(replace))
    print(f"{name} = {expr_str}")
    options[name] = expr_str

options.write("data/BOUT.inp", overwrite=True)
