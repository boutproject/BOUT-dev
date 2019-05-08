#!/usr/bin/env python

from __future__ import print_function
from builtins import str
from boututils.run_wrapper import shell, launch
from boutdata.collect import collect

import numpy as np

print("Making 3-field ELM simulation")
shell("make > make.log")

# Run simulation
nproc = 2

s, out = launch("./elm_pb ", nproc=nproc, pipe=True)

with open("run.log", "w") as f:
    f.write(out)

# Get time base
tarr = collect("t_array", path="data")
nt = len(tarr)

# Read pressure
p = collect("P", path="data")

# Calculate RMS in toroidal direction
prms = np.sqrt(np.mean(p**2, axis=3))

growth = np.gradient(np.log(prms[:,42,32]))

# Final growth-rate
gamma = growth[-2]

import matplotlib.pyplot as plt

plt.plot(tarr, prms[:,42,32], label='Outboard midplane')
plt.plot( [tarr[0], tarr[-1]],
	  [prms[-1,42,32]*np.exp(gamma*(tarr[0] - tarr[-1])), prms[-1,42,32]], '--', label=r'$\gamma =$'+str(gamma))

plt.yscale('log')
plt.grid()

plt.xlabel(r"Time [$1/\tau_A$]")
plt.ylabel("Pressure RMS")

plt.legend(loc="upper left")

plt.savefig("growth.pdf")

plt.show()

############### Poloidal plot

plt.figure()

# Take a poloidal slice at fixed toroidal angle
from boutdata.pol_slice import pol_slice
p2d = pol_slice(p[-1,:,:,:], 'cbm18_dens8.grid_nx68ny64.nc', n=15, zangle=0.0)

# Read grid file to get coordinates
from boututils.datafile import DataFile
g = DataFile('cbm18_dens8.grid_nx68ny64.nc')

Rxy = g.read("Rxy") # Major radius [m]
Zxy = g.read("Zxy") # Height [m]

plt.contourf(Rxy, Zxy, p2d, 30)
plt.axis('equal')  # Maintain aspect ratio

plt.colorbar()   # Plot a bar down the side with a color scale

plt.savefig("poloidal_slice.pdf")

plt.show()
