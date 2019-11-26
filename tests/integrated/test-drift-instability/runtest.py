#!/usr/bin/env python3

#
# Run the test, compare results against the benchmark
#

from __future__ import print_function
from __future__ import division

from builtins import str, range

from math import isnan

nproc = 2       # Number of processors to run on

# Relative tolerance in frequency and growth rate
omega_tol = 1e-2
gamma_tol = 1e-2

from boututils.run_wrapper import shell, shell_safe, launch_safe
from boututils.file_import import file_import
from boututils.calculus import deriv
from boututils.linear_regression import linear_regression

from boutdata.collect import collect
import numpy as np
from sys import exit ,argv

nthreads=1


print("Making resistive drift instability test")
shell_safe("make > make.log")

zlist          = [2, 32, 256]  # Just test a few

# Values from revision c4f7ec92786b333a5502c5256b5e602ba867090f
# 26th Oct 2011
omega_orig     = {1:1.00819536681,
                  2:0.975756752742,
                  4:0.888269347989,
                  8:0.741253948437,
                  16:0.574003058717,
                  32:0.428644078148, #0.430745728591, Changed 25th April 2014, revision fd032da
                  64:0.336073734175,
                  128:0.2208085,
                  256:0.155890075514} #0.156853396108} Changed 25th April 2014
gamma_orig     = {1:0.0784576199501,
                  2:0.145174416274, #0.143251663975, Changed 25th April 2014, revision fd032da
                  4:0.226979299615,
                  8:0.28632661549,
                  16:0.295998650267,
                  32:0.271288568509,
                  64:0.222344013151,
                  128:0.1716229,
                  256:0.12957680451} #0.130220286897} Changed 25th April 2014

#zlist = map(lambda x:2**x, range(9))

# Create a directory for the data
shell_safe("mkdir -p data")

# Import the grid file
grid = file_import("uedge.grd_std.cdl")

code = 0 # Return code
for zeff in zlist:
    # Create the input file, setting Zeff
    timestep = 5e3
    if zeff < 128:
        # reduce time-step. At large times these cases produce noise
        timestep = 1e3

    # Delete old output files
    shell("rm data/BOUT.dmp.*.nc")

    print("Running drift instability test, zeff = ", zeff)

    # Run the case
    s, out = launch_safe("./2fluid 2fluid:Zeff={} timestep={}"
                         .format(zeff, timestep),
                         nproc=nproc, mthread=nthreads, pipe=True)
    f = open("run.log."+str(zeff), "w")
    f.write(out)
    f.close()

    # Collect data
    Ni = collect("Ni", path="data", xind=2, yind=20, info=False)
    phi = collect("phi", path="data", xind=2, yind=20, info=False)

    zmax     = collect("ZMAX", path="data", info=False)
    rho_s    = collect("rho_s", path="data", info=False)
    wci      = collect("wci", path="data", info=False)
    t_array  = collect("t_array", path="data", info=False)

    dims = np.shape(Ni)
    nt = dims[0]
    ny = dims[2]
    nz = dims[3]

    ##### Calculate geometric and physical quantities
    lZeta  = 1e2*zmax*2*np.pi*grid['R0']    # toroidal range [cm]
    lbNorm = lZeta*(grid['Bpxy'][0, ny // 2] / grid['Bxy'][0, ny // 2])  # binormal coord range [cm]
    zPerp  = lbNorm*np.arange(nz)/(nz-1)    # binormal coordinate [cm]

    cLight = 3e10                        # speed of light [cm/s]
    vTe    = 4.2e7*np.sqrt(grid['Te_x']) # electron thermal speed [cm/s]
    kperp  = 2*np.pi/lbNorm              # binormal wavenumber, [cm-1]
    wce    = 1.76e7*1e4*grid['bmag']     # electron cyclotron frequency, [rad/s]

    # Ni scale length [cm]
    Ln = np.mean(np.abs(grid['Ni0'][:, ny // 2] /
                        deriv(grid['Rxy'][:, ny // 2]*1e2, grid['Ni0'][:, ny // 2])))

    vPe   = (vTe)**2 / (wce*Ln) # electron diamagnetic drift speed [cm/s]
    wstar = vPe*kperp

    logLam = 24. - np.log(np.sqrt(grid['Ni_x']*1e14 / grid['Te_x']))
    nuei   = zeff*2.91e-6*(grid['Ni_x']*1e14)*logLam/(grid['Te_x'])**1.5
    #wci=9.58e3*(1./d.AA)*1e4*du.Bmag

    lpar = np.sum(((grid['Bxy']/grid['Bpxy']))*grid['dlthe']) / grid['nx'] # [m], average over flux surfaces
    kpar = 2*np.pi/(1e2*lpar) # cm-1
    spar = (kpar / kperp)**2 * wci * wce / (0.51 * nuei) # [1/s]
    sparn = spar / wstar

    wpe   = 5.64e4*np.sqrt(1e14*grid['Ni_x']) # electron plasma frequency, [rad/s]
    mu    = (cLight*kperp/wpe)**2
    sperp = (0.51 * nuei) * mu # [1/s]

    ##### Analyse data

    nt0 = 15  # Skip initial part of the curve (transients)

    # Measure exponential growth-rate using maximum value
    maxVal = np.zeros(nt - nt0)

    # Track the motion of the peak to infer phase velocity
    peak = np.zeros(nt-nt0)

    for t in range(nt0, nt):
        ind = np.argmax(Ni[t,0,0,:])  # Index of maximum value
        maxVal[t-nt0] = Ni[t,0,0,ind] # Store maximum value (for gamma)
        # Refine maximum location by fitting a quadratic through 3 points
        c = Ni[t,0,0,ind]
        m = Ni[t,0,0,(ind-1) % nz] # Python '%' always positive
        p = Ni[t,0,0,(ind+1) % nz]
        # Fit y = c + ax + bx**2
        a = 0.5*(p-m)
        b = p - (c + a)
        peak[t-nt0] = ind - 0.5*a/b # Maximum

    # Check for periodic recurrence
    if peak[1] > peak[0]:
        # Increasing; watch for wrapping around the top
        for i in range(nt-nt0):
            if peak[i] < peak[i-1]:
                peak[i:(nt-nt0)] = peak[i:(nt-nt0)] + nz
    else:
        # Decreasing; watch for wrapping around the bottom
        for i in range(nt-nt0):
            if peak[i] > peak[i-1]:
                peak[i:(nt-nt0)] = peak[i:(nt-nt0)] - nz

    # Fit y = a + gamma*x
    a, gamma = linear_regression(t_array[nt0:nt] / wci, np.log(maxVal))

    # Get phase velocity
    a, Vphase = linear_regression(t_array[nt0:nt] / wci, peak*lbNorm/(nz-1))

    # Calculate normalised quantities
    omega = np.abs(Vphase)*kperp/wstar
    gamma = gamma / wstar

    # Calculate analytic result
    t = 0.5*( np.sqrt(sparn**4 + 16*sparn**2) - sparn**2 )
    wr = 0.5*np.sqrt(t)
    wi = sparn / np.sqrt(t) - 0.5*sparn

    try:
        origr = omega_orig[zeff]
        origi = gamma_orig[zeff]

        omegadiff = abs(omega - origr) / origr
        gammadiff = abs(gamma - origi) / origi
    except:
        origr = None
        origi = None
        omegadiff = None
        gammadiff = None

    print("  Normalised omega = ", omega, " analytic = ", wr, " original = ", origr, " (", 100.*omegadiff,"%)")
    print("  Normalised gamma = ", gamma, " analytic = ", wi, " original = ", origi, " (",100.*gammadiff,"%)")

    if omegadiff != None:
        if isnan(omegadiff) or (omegadiff > omega_tol) or (gammadiff > gamma_tol):
            code = 1 # Failed test
            print("  => FAILED")
        else:
            print("  => PASSED")
    else:
        print("  => No original to compare against")

exit(code)
