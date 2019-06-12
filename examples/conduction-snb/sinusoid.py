#
# SNB model example
# A small sinusoidal temperature perturbation is added to a periodic 1D domain

import numpy as np
import matplotlib.pyplot as plt

from boututils.run_wrapper import shell_safe, launch_safe
from boutdata.collect import collect
from sys import exit

shell_safe("make > make.log")

# Electron temperature in eV
Telist = 10 ** np.linspace(0,3,20)

# Electron density in m^-3
Ne = 1e20

# Length of the domain in m
length = 1.0

c   = 299792458
mu0 = 4.e-7*np.pi
e0  = 1/(c*c*mu0)
qe = 1.602176634e-19
me = 9.10938356e-31

thermal_speed = np.sqrt(2.*qe  * Telist / me)
Y = (qe**2 / (e0 * me))**2 / (4 * np.pi)
coulomb_log = 6.6 - 0.5 * np.log(Ne * 1e-20) + 1.5 * np.log(Telist)
lambda_ee_T = thermal_speed**4 / (Y * Ne  * coulomb_log)

beta_max_list = [5, 10, 20, 40]
colors = ['k','b','g','r']

ngroups_list = [20, 40, 80]
syms = ['x', 'o', 'D']

for beta_max, color in zip(beta_max_list, colors):
    for ngroups, sym in zip(ngroups_list, syms):

        flux_ratio = []
        for Te in Telist:
            cmd = "./conduction-snb \"Te={0}+0.01*sin(y)\" Ne={1} mesh:length={2} snb:beta_max={3} snb:ngroups={4}".format(Te, Ne, length, beta_max, ngroups)

            # Run the case
            s, out = launch_safe(cmd, nproc=1, mthread=1, pipe=True)

            div_q = collect("Div_Q", path="data").ravel()
            div_q_SH = collect("Div_Q_SH", path="data").ravel()

            # Get the index of maximum S-H heat flux
            ind = np.argmax(div_q_SH)

            flux_ratio.append(div_q[ind] / div_q_SH[ind])

        plt.plot(lambda_ee_T / length, flux_ratio, '-'+sym+color, label=r"$\beta_{{max}}={0}, N_g={1}$".format(beta_max,ngroups))

plt.legend()
plt.xlabel(r"$\lambda_{ee,T} / L$")
plt.ylabel(r"$q / q_{SH}$")
plt.xscale("log")

plt.savefig("snb-sinusoidal.png")
plt.show()
