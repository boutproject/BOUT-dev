#!/bin/python3
# -*- coding: utf-8

#!##################################################################
#
#       2D blob simulations
#
#  Copyright:   NR Walkden, B Dudson, D Schwörer;  2012, , 2017
#
###################################################################

import boutpp as bc
from numpy import sqrt
from boutpp import bracket, DDZ, Delp2
import sys

bc.init("-d blob".split(" ") + sys.argv[1:])
mesh = bc.Mesh.getGlobal()
n = bc.Field3D.fromMesh(mesh)
omega = bc.Field3D.fromMesh(mesh)
phi = bc.Field3D.fromMesh(mesh)

phiSolver = bc.Laplacian()

options = bc.Options("model")
# Temperature in eV
Te0 = options.get("Te0", 30)
e = options.get("e", 1.602e-19)
m_i = options.get("m_i", 2 * 1.667e-27)
m_e = options.get("m_e", 9.11e-31)

# Background density in cubic m
n0 = options.get("n0", 1e19)
# Viscous diffusion coefficient
D_vort = options.get("D_vort", 0)
# Density diffusion coefficient
D_n = options.get("D_n", 0)

# Radius of curvature [m]
R_c = options.get("R_c",   1.5)
# Parallel connection length [m]
L_par = options.get("L_par", 10)

# Value of magnetic field strength [T]
B0 = options.get("B0", 0.35)


# System option switches

# Include compressible ExB term in density equation
compressible = options.get("compressible", False)
# Use Boussinesq approximation in vorticity
boussinesq = options.get("boussinesq", True)
# Sheath closure
sheath = options.get("sheath", True)


Omega_i = e * B0 / m_i           # Cyclotron Frequency
c_s = sqrt(e * Te0 / m_i)      # Bohm sound speed
rho_s = c_s / Omega_i          # Bohm gyro-radius

print("\n\n\t----------Parameters: ------------ \n\tOmega_i = %e /s,\n\t"
      "c_s = %e m/s,\n\trho_s = %e m\n" % (Omega_i, c_s, rho_s))

# Calculate delta_*, blob size scaling
print("\tdelta_* = rho_s * (dn/n) * %e "
      % (pow(L_par * L_par / (R_c * rho_s), 1. / 5)))

# /************ Create a solver for potential ********/

if boussinesq:
    # BOUT.inp section "phiBoussinesq"
    phiSolver = bc.Laplacian(bc.Options("phiBoussinesq"))
else:
        # BOUT.inp section "phiSolver"
    phiSolver = bc.Laplacian(bc.Options("phiSolver"))

phi = bc.Field3D.fromMesh(mesh)
phi.set(0.0)  # Starting guess for first solve (if iterative)

# /************ Tell BOUT++ what to solve ************/

model = bc.PhysicsModelBase()
model.solve_for(n=n, omega=omega)
# model.save_repeat(phi=phi)
# model.save_once(rho_s=rho_s,c_s=c_s,Omega_i=Omega_i)


def rhs(time):
    global n, omega, phi
    # Run communications
    ######################################
    mesh.communicate(n, omega)

    # Invert div(n grad(phi)) = grad(n) grad(phi) + n Delp_perp^2(phi) = omega
    ######################################
    # Set the time derivative by adding/... to it
    # make sure to never overwrite it
    # ddt_n = bla does NOT set the time derivative
    ddt_n = n.ddt()
    ddt_n.set(0)
    if not boussinesq:
        # Including full density in vorticit inversion
          # Update the 'C' coefficient. See invert_laplace.hxx
        phiSolver.setCoefC(n)
        # Use previous solution as guess
        phi = phiSolver.solve(omega / n, phi)
    else:
        # Background density only (1 in normalised units)
        phi = phiSolver.solve(omega, phi)

    mesh.communicate(phi)

    # Density Evolution
    # /

    ddt_n += -bracket(phi, n, "BRACKET_SIMPLE")    # ExB term
    ddt_n += 2 * DDZ(n) * (rho_s / R_c)               # Curvature term
    ddt_n += D_n * Delp2(n)                       # Diffusion term

    if compressible:
        ddt_n -= 2 * n * DDZ(phi) * (rho_s / R_c)       # ExB Compression term

    if sheath:
        # Sheath closure
        ddt_n += n * phi * (rho_s / L_par)  # - (n - 1)*(rho_s/L_par)

    # Vorticity evolution
    # /
        # ExB term
    ddt_omega = -bracket(phi, omega, "BRACKET_SIMPLE")
    ddt_omega += 2 * DDZ(n) * (rho_s / R_c) / n
    # Viscous diffusion term
    ddt_omega += D_vort * Delp2(omega) / n

    if sheath:
        ddt_omega += phi * (rho_s / L_par)
    # other option to set time derivaitve:
    # create a field and set it in the end
    omega.ddt(ddt_omega)

# Start the simulation
model.setRhs(rhs)
model.solve()
