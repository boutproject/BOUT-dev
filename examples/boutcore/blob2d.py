#!/bin/python3
# -*- coding: utf-8

#!##################################################################
#
#       2D blob simulations
#
#  Copyright:   NR Walkden, B Dudson, D Schwörer;  2012, 2017, 2018
#
###################################################################

import boutcore as bc
from numpy import sqrt
from boutcore import bracket, DDZ, Delp2, PhysicsModel
import sys
import os


class Blob2D(PhysicsModel):
    def init(self, restart):
        self.mesh = bc.Mesh.getGlobal()
        self.n = bc.Field3D.fromMesh(self.mesh)
        self.omega = bc.Field3D.fromMesh(self.mesh)

        self.phiSolver = bc.Laplacian()

        options = bc.Options("model")
        # Temperature in eV
        Te0 = options.get("Te0", 30)
        e = options.get("e", 1.602e-19)
        m_i = options.get("m_i", 2 * 1.667e-27)
        m_e = options.get("m_e", 9.11e-31)

        # Viscous diffusion coefficient
        self.D_vort = options.get("D_vort", 0)
        # Density diffusion coefficient
        self.D_n = options.get("D_n", 0)

        # Radius of curvature [m]
        self.R_c = options.get("R_c", 1.5)
        # Parallel connection length [m]
        self.L_par = options.get("L_par", 10)

        # Value of magnetic field strength [T]
        B0 = options.get("B0", 0.35)

        # System option switches

        # Include compressible ExB term in density equation
        self.compressible = options.get("compressible", False)
        # Use Boussinesq approximation in vorticity
        self.boussinesq = options.get("boussinesq", True)
        # Sheath closure
        self.sheath = options.get("sheath", True)

        Omega_i = e * B0 / m_i  # Cyclotron Frequency
        c_s = sqrt(e * Te0 / m_i)  # Bohm sound speed
        self.rho_s = c_s / Omega_i  # Bohm gyro-radius

        print(
            "\n\n\t----------Parameters: ------------ \n\tOmega_i = %e /s,\n\t"
            "c_s = %e m/s,\n\trho_s = %e m\n" % (Omega_i, c_s, self.rho_s)
        )

        # Calculate delta_*, blob size scaling
        print(
            "\tdelta_* = rho_s * (dn/n) * %e "
            % (pow(self.L_par * self.L_par / (self.R_c * self.rho_s), 1.0 / 5))
        )

        # /************ Create a solver for potential ********/

        if self.boussinesq:
            # BOUT.inp section "phiBoussinesq"
            self.phiSolver = bc.Laplacian(bc.Options("phiBoussinesq"))
        else:
            # BOUT.inp section "phiSolver"
            self.phiSolver = bc.Laplacian(bc.Options("phiSolver"))

        # Starting guess for first solve (if iterative)
        self.phi = bc.create3D("0")

        # /************ Tell BOUT++ what to solve ************/

        self.solve_for(n=self.n, omega=self.omega)
        # model.save_repeat(phi=phi)
        # model.save_once(rho_s=rho_s,c_s=c_s,Omega_i=Omega_i)

    def rhs(self, time):
        # Run communications
        ######################################
        self.mesh.communicate(self.n, self.omega)

        # Invert div(n grad(phi)) = grad(n) grad(phi) + n Delp_perp^2(phi) = omega
        ######################################
        # Set the time derivative by adding/... to it
        # make sure to never overwrite it
        # ddt_n = bla does NOT set the time derivative
        ddt_n = self.n.ddt()
        ddt_n.set(0)
        if not self.boussinesq:
            # Including full density in vorticit inversion
            # Update the 'C' coefficient. See invert_laplace.hxx
            self.phiSolver.setCoefC(n)
            # Use previous solution as guess
            self.phi = self.phiSolver.solve(omega / n, self.phi)
        else:
            # Background density only (1 in normalised units)
            self.phi = self.phiSolver.solve(self.omega, self.phi)

        self.mesh.communicate(self.phi)

        # Density Evolution
        # /

        # ExB term
        ddt_n += -bracket(self.phi, self.n, "BRACKET_SIMPLE")
        # Curvature term
        ddt_n += 2 * DDZ(self.n) * (self.rho_s / self.R_c)
        # Diffusion term
        ddt_n += self.D_n * Delp2(self.n)

        if self.compressible:
            # ExB Compression term
            ddt_n -= 2 * self.n * DDZ(self.phi) * (self.rho_s / self.R_c)

        if self.sheath:
            # Sheath closure
            ddt_n += (
                self.n * self.phi * (self.rho_s / self.L_par)
            )  # - (n - 1)*(rho_s/L_par)

        # Vorticity evolution
        # /
        # ExB term
        ddt_omega = -bracket(self.phi, self.omega, "BRACKET_SIMPLE")
        ddt_omega += 2 * DDZ(self.n) * (self.rho_s / self.R_c) / self.n
        # Viscous diffusion term
        ddt_omega += self.D_vort * Delp2(self.omega) / self.n

        if self.sheath:
            ddt_omega += self.phi * (self.rho_s / self.L_par)
        # other option to set time derivaitve:
        # create a field and set it in the end
        self.omega.ddt(ddt_omega)


# Ensure the blob folder exists
def ensure_blob():
    if not os.path.isdir("blob"):
        print("Setting up folder blob for simulation ...")
        os.mkdir("blob")
    if not os.path.exists("blob/BOUT.inp"):
        with open("blob/BOUT.inp", "w") as f:
            f.write(
                """\
# settings file for BOUT++
#
# Blob simulation in a 2D slab
#
# This case has blob size
#
# delta = 0.3*256 ~ 10 * delta_*


# settings used by the core code

NOUT = 50      # number of time-steps
TIMESTEP = 50  # time between outputs [1/wci]


MXG = 2      # Number of X guard cells
MYG = 0      # No y derivatives, so no guard cells needed in y

[mesh]

nx = 260    # Note: 4 guard cells
ny = 1
nz = 256

dx = 0.3      # Grid spacing [rho_s]
dz = 0.3

##################################################
# derivative methods

[mesh:ddx]

first = C2
second = C2
upwind = W3

[mesh:ddy]

first = C2
second = C2
upwind = W3

[mesh:ddz]

first = FFT
second = FFT
upwind = W3

###################################################
# Time-integration solver

[solver]

ATOL = 1.0e-10  # absolute tolerance
RTOL = 1.0e-5   # relative tolerance
mxstep = 10000  # Maximum internal steps per output

###################################################
# Electrostatic potential solver
# These options are used if boussinesq = false

[phiSolver]
type = petsc  # Needed if Boussinesq = false
pctype = user  # Preconditioning type

fourth_order = true  # 4th order or 2nd order

flags = 0  # inversion flags for phi
             # 0  = Zero value
             # 10 = Zero gradient AC inner & outer
             # 15 = Zero gradient AC and DC
             # 768 = Zero laplace inner & outer

[phiSolver:precon]  # Preconditioner (if pctype=user)
filter     = 0.     # Must not filter solution
flags      = 49152  # set_rhs i.e. identity matrix in boundaries

###################################################
# Electrostatic potential solver (Boussinesq)

[phiBoussinesq]
# By default type is tri (serial) or spt (parallel)
flags = 0

##################################################
# general settings for the model

[model]

Te0 = 5    # Electron Temperature (eV)

n0 = 2e18  # Background plasma density (m^-3)

compressible = false  # Compressibility?

boussinesq = true  # Boussinesq approximation (no perturbed n in vorticity)

D_vort = 1e-6  # Viscosity
D_n = 1e-6    # Diffusion

R_c = 1.5  # Radius of curvature (m)

# settings for individual variables
# The section "All" defines default settings for all variables
# These can be overridden for individual variables in
# a section of that name.

[All]
scale = 0.0 # default size of initial perturbations

bndry_all = neumann # Zero-gradient on all boundaries

[n]  # Density
scale = 1.0 # size of perturbation

height = 0.5
width = 0.05

function = 1 + height * exp(-((x-0.25)/width)^2 - ((z/(2*pi) - 0.5)/width)^2)
"""
            )


if __name__ == "__main__":
    if "--create" in sys.argv:
        sys.argv.remove("--create")
        ensure_blob()
    bc.init("-d blob".split(" ") + sys.argv[1:])

    # Create an instance
    blob2d = Blob2D()
    # Start the simulation
    blob2d.solve()
