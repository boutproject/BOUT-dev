#!/bin/python3
# -*- coding: utf-8

#!##################################################################
#
#       2D blob simulations
#
#  Copyright:   NR Walkden, B Dudson, D Schwörer;  2012, 2017, 2018
#
###################################################################

# from boutcore import *
# init("-d mini -q -q -q")

# class MyModel(PhysicsModel):
#          def init(self,restart):
#                      self.n=create3D("dens:function")
#                              self.solve_for(dens=self.n)
#                                  def rhs(self,time):
#                                              self.n.ddt(DDX(self.n))

#                                              model=MyModel()
#                                              model.solve()

import boutcore as bc
from numpy import sqrt
from boutcore import bracket, DDZ, Delp2, PhysicsModel
import sys
bc.init("-d blob".split(" ") + sys.argv[1:])


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
        self.R_c = options.get("R_c",   1.5)
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

        Omega_i = e * B0 / m_i           # Cyclotron Frequency
        c_s = sqrt(e * Te0 / m_i)      # Bohm sound speed
        self.rho_s = c_s / Omega_i          # Bohm gyro-radius

        print("\n\n\t----------Parameters: ------------ \n\tOmega_i = %e /s,\n\t"
              "c_s = %e m/s,\n\trho_s = %e m\n" % (Omega_i, c_s, self.rho_s))

        # Calculate delta_*, blob size scaling
        print("\tdelta_* = rho_s * (dn/n) * %e "
              % (pow(self.L_par * self.L_par / (self.R_c * self.rho_s), 1. / 5)))

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
            ddt_n += self.n * self.phi * \
                (self.rho_s / self.L_par)  # - (n - 1)*(rho_s/L_par)

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


# Create an instance
blob2d = Blob2D()
# Start the simulation
blob2d.solve()
