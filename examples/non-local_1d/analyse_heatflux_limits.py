#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Runs the conduction example, produces some output
#

from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div

from boutdata.collect import collect
from math import sqrt, pow, pi
from matplotlib import pyplot

nproc = 1  # Number of processors to use

def sign(x):
  if x>=0.:
    return 1.
  else:
    return -1.

massunit = 1.602176565e-19
electron_mass = old_div(9.10938291e-31,massunit)
electron_charge = -1.602176565e-19
ion_mass =  old_div(3.34358348e-27,massunit)
ion_charge = 1.602176565e-19
epsilon_0 = old_div(8.85418781762039e-12,pow(massunit,-1))
logLambda = 16.0

stagger = True
t_interval = 10

# Collect the data
Te = collect("T_electron", path="data", info=True, yguards=True)
Ti = collect("T_ion", path="data", info=True, yguards=True)
n = collect("n_ion", path="data", info=True, yguards=True)

ylength = len(Te[0,2,:,0])
tlength = len(Te[:,2,0,0])

try: dy = collect("dy", path="data", info=True, yguards=True)
except TypeError:
  print("Warning, could not find dy, setting to 1.0")
  dy=[[1.]*ylength]*5

try: g_22 = collect("g_22", path="data", info=True, yguards=True)
except TypeError:
  print("Warning, could not find g_22, setting to (80./256)^2")
  g_22=[[pow(old_div(80.,256),2)]*ylength]*5

try:
    qin = collect("heat_flux", path="data", info=True, yguards=True)
    q = [[0.]*ylength for i in range(0,old_div(tlength,t_interval)+1)]
    for i in range(0,tlength,t_interval):
      for j in range(2,ylength-2):
        q[old_div(i,t_interval)][j] = qin[i,2,j,0]
except TypeError:
  print("Calculating Braginskii heat flux")
  q = [[0.]*ylength for i in range(0,old_div(tlength,t_interval)+1)]
  tau_ei = 0
  gradT = 0
  for i in range(0,tlength,t_interval):
    for j in range(2,ylength-2):
      tau_ei = 3 * pow(pi,1.5) * pow(epsilon_0,2) * sqrt(electron_mass) * pow(2.,1.5) * pow(Te[i,2,j,0],1.5) / n[i,2,j,0] / pow(electron_charge,2) / pow(ion_charge,2) / logLambda
      gradT = (Te[i,2,j-2,0] - 8.*Te[i,2,j-1,0] + 8.*Te[i,2,j+1,0] - Te[i,2,j+2,0])/12./dy[2][j]/sqrt(g_22[2][j])
      q[old_div(i,t_interval)][j] = -3.16*n[i,2,j,0]*Te[i,2,j,0]*tau_ei/electron_mass*gradT
  stagger = False

Temax = [0.]*(old_div(tlength,t_interval)+1)
for i in range(0,tlength,t_interval):
  Temax[old_div(i,t_interval)] = max(Te[i,2,:,0])
Timax = [0.]*(old_div(tlength,t_interval)+1)
for i in range(0,tlength,t_interval):
  Timax[old_div(i,t_interval)] = max(Ti[i,2,:,0])
nmax = [0.]*(old_div(tlength,t_interval)+1)
for i in range(0,tlength,t_interval):
  nmax[old_div(i,t_interval)] = max(n[i,2,:,0])
freestreammax = [0.]*(old_div(tlength,t_interval)+1)
for i in range(old_div(tlength,t_interval)+1):
  freestreammax[i]=3./2.*nmax[i]*Temax[i]*sqrt(old_div((Temax[i]+Timax[i]),ion_mass))
  #freestreammax[i]=3./2.*nmax[i]*Temax[i]*sqrt(2*Temax[i]/electron_mass)

#freestream = [[i]*ylength for i in range(tlength/t_interval+1)]
#for i in range(0,tlength,t_interval):
  #for j in range(2,ylength-2):
    #freestream[i/t_interval][j] = sign(j-ylength/2)*3./2.*n[i,2,j,0]*Te[i,2,j,0]*sqrt((Te[i,2,j,0]+Ti[i,2,j,0])/ion_mass)

#freestream = [[i]*ylength for i in range(tlength/t_interval+1)]
#for i in range(0,tlength,t_interval):
  #for j in range(2,ylength-2):
    #freestream[i/t_interval][j] = sign(j-ylength/2)*3./2.*n[i,2,j,0]*Te[i,2,j,0]*sqrt(2*Te[i,2,j,0]/electron_mass)

#for i in range(0,tlength,t_interval):
  #pyplot.figure(i)
  #pyplot.plot(range(2,ylength-2),q[i/t_interval][2:-2],'k', range(2,ylength-2),freestream[i/t_interval][2:-2], 'r')

for i in range(0,tlength,t_interval):
  pyplot.figure(i)
  pyplot.plot(list(range(2,ylength-2)),q[old_div(i,t_interval)][2:-2],'k', list(range(2,ylength-2)),[freestreammax[old_div(i,t_interval)]]*(ylength-4), 'r',[-freestreammax[old_div(i,t_interval)]]*(ylength-4), 'r')

pyplot.show()

#pyplot.figure(5)
#pyplot.plot(logtime,q_electron_left,'r',logtime,q_ion_left,'b',logtime,q_total_left,'k')
#pyplot.title("Electron (red), Ion (blue) and Total (black) Heat Flux vs log(t)")
#pyplot.figure(6)
#pyplot.plot(q_electron_left,'r',q_ion_left,'b',q_total_left,'k')
#pyplot.title("Electron (red), Ion (blue) and Total (black) Heat Flux")
#pyplot.xlabel(u"t/μs")
#pyplot.ylabel("Q/W.m$^{-2}$")
#pyplot.show()
