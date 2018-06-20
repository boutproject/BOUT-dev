#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Runs the conduction example, produces some output
#

from __future__ import print_function
from builtins import str
from builtins import range

from boutdata.collect import collect
from sys import argv
from math import sqrt, log10, log, pi
from matplotlib import pyplot, ticker, rc

nproc = 1  # Number of processors to use

rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern']})
#rc('font',**{'family':'serif','serif':['Computer Modern'],'size':16})


gamma = float(argv[1])
if len(argv)==2:
  end_index = -1
  data_path = "data"
elif len(argv)==3:
  try:
    end_index = int(argv[2])
    data_path = "data"
  except ValueError:
    end_index = -1
    data_path = str(argv[2])
elif len(argv)==4:
  end_index = int(argv[2])
  data_path = str(argv[3])
else:
  print("Arguments: 'gamma [end_index] [data_path]' or 'gamma [data_path]'")
  Exit(1)

electron_mass = 9.10938291e-31
ion_mass =  3.34358348e-27

# Collect the data
Te = collect("T_electron", path=data_path, info=True, yguards=True)
Ti = collect("T_ion", path=data_path, info=True, yguards=True)
n = collect("n_ion", path=data_path, info=True, yguards=True)

if end_index<0:
  end_index = len(n[:,0,0,0])

logtime = []
q_electron_left = []
q_electron_right = []
q_ion_left = []
q_ion_right = []
q_total_left = []
q_total_right = []
q_target_electron_left = []
q_target_electron_right = []
q_target_ion_left = []
q_target_ion_right = []
q_target_total_left = []
q_target_total_right = []
right_index = len(Te[0,2,:,0])-4

for i in range(len(Te[:end_index,2,0,0])):
  if i>0: logtime.append(log10(i))
  else: logtime.append(0)
  Te_here = 0.5*(Te[i,2,2,0]+Te[i,2,3,0])
  n_here = 0.5*(n[i,2,2,0]+n[i,2,3,0])
  Ti_here = 0.5*(Ti[i,2,2,0]+Ti[i,2,3,0])
  sheath_potential = 0.5*Te_here*log(2*pi*electron_mass/ion_mass*(1+gamma*Ti_here/Te_here))
  q_electron_left.append((2.0*Te_here-sheath_potential)*n_here*1.602176565e-19*sqrt((Te_here+gamma*Ti_here)/3.34358348e-27*1.602176565e-19)) # in W/m^2
  q_ion_left.append(n_here*((2.5+0.5*gamma)*Ti_here+0.5*Te_here)*1.602176565e-19*sqrt((Te_here+gamma*Ti_here)/3.34358348e-27*1.602176565e-19)) # in W/m^2
  q_total_left.append(q_electron_left[i]+q_ion_left[i]) # in W/m^2

  q_target_electron_left.append((2.0*Te_here)*n_here*1.602176565e-19*sqrt((Te_here+gamma*Ti_here)/3.34358348e-27*1.602176565e-19)) # in W/m^2
  q_target_ion_left.append(n_here*((2.5+0.5*gamma)*Ti_here+0.5*Te_here-sheath_potential)*1.602176565e-19*sqrt((Te_here+gamma*Ti_here)/3.34358348e-27*1.602176565e-19)) # in W/m^2
  q_target_total_left.append(q_target_electron_left[i]+q_target_ion_left[i]) # in W/m^2

  Te_here = (Te[i,2,right_index,0]+Te[i,2,right_index+1,0])
  n_here = (n[i,2,right_index,0]+n[i,2,right_index+1,0])
  Ti_here = (Ti[i,2,right_index,0]+Ti[i,2,right_index+1,0])
  sheath_potential = 0.5*Te_here*log(2*pi*electron_mass/ion_mass*(1+gamma*Ti_here/Te_here))
  q_electron_right.append((2.0*Te_here-sheath_potential)*n_here*1.602176565e-19*sqrt((Te_here+gamma*Ti_here)/3.34358348e-27*1.602176565e-19)) # in W/m^2
  q_ion_right.append(n_here*((2.5+0.5*gamma)*Ti_here+0.5*Te_here)*1.602176565e-19*sqrt((Te_here+gamma*Ti_here)/3.34358348e-27*1.602176565e-19)) # in W/m^2
  q_total_right.append(q_electron_right[i]+q_ion_right[i]) # in W/m^2

  q_target_electron_right.append((2.0*Te_here)*n_here*1.602176565e-19*sqrt((Te_here+gamma*Ti_here)/3.34358348e-27*1.602176565e-19)) # in W/m^2
  q_target_ion_right.append(n_here*((2.5+0.5*gamma)*Ti_here+0.5*Te_here-sheath_potential)*1.602176565e-19*sqrt((Te_here+gamma*Ti_here)/3.34358348e-27*1.602176565e-19)) # in W/m^2
  q_target_total_right.append(q_target_electron_right[i]+q_target_ion_right[i]) # in W/m^2

#pyplot.figure(1)
#pyplot.plot(logtime,q_electron_left)
#pyplot.axes().set_aspect('auto')
#pyplot.title("Electron Heat Flux (lower boundary)")
#pyplot.figure(2)
#pyplot.plot(q_electron_right)
#pyplot.axes().set_aspect('auto')
#pyplot.title("Electron Heat Flux (upper boundary)")
#pyplot.figure(3)
#pyplot.plot(q_total_left)
#pyplot.axes().set_aspect('auto')
#pyplot.title("Total Heat Flux (lower boundary)")
#pyplot.figure(4)
#pyplot.plot(q_total_right)
#pyplot.axes().set_aspect('auto')
#pyplot.title("Total Heat Flux (upper boundary)")
pyplot.figure(5,dpi=80, facecolor='w')
#pyplot.plot(logtime,q_electron_left,'r',logtime,q_ion_left,'b',logtime,q_total_left,'k')
#pyplot.title("Electron (red), Ion (blue) and Total (black) Sheath Edge Heat Flux vs log(t)")
pyplot.semilogx(q_electron_left,'r',q_ion_left,'b',q_total_left,'k')
pyplot.title("Electron (red), Ion (blue) and Total (black) Sheath Edge Heat Flux\\ \\")
pyplot.xlabel("$t/\mu\mathrm{s}$")
pyplot.ylabel(r"$Q\mathrm{/W.m}^{-2}$")
pyplot.axes().xaxis.set_major_formatter(ticker.FormatStrFormatter("$%g$"))
pyplot.axes().grid(color='grey', which='both')
#pyplot.tight_layout(pad=20)
pyplot.figure(6)
pyplot.plot(q_electron_right,'r',q_ion_right,'b',q_total_right,'k')
pyplot.title("Electron (red), Ion (blue) and Total (black) Sheath Edge Heat Flux\\ \\")
pyplot.xlabel("$t/\mu\mathrm{s}$")
pyplot.ylabel(r"$Q\mathrm{/W.m}^{-2}$")
#pyplot.figure(7,dpi=80, facecolor='w')
pyplot.figure(7,dpi=800, facecolor='w')
#pyplot.plot(logtime,q_target_electron_left,'r',logtime,q_target_ion_left,'b',logtime,q_target_total_left,'k')
#pyplot.title("Electron (red), Ion (blue) and Total (black) Target Heat Flux vs log(t)")
#pyplot.semilogx(q_target_electron_left,'r',q_target_ion_left,'b',q_target_total_left,'k')
#pyplot.semilogx(q_target_electron_left,'k',q_target_ion_left,'r',q_target_total_left,'b')
pyplot.semilogx(q_target_electron_left,'r--',q_target_ion_left,'b:',q_target_total_left,'k')
#pyplot.title("Electron (red), Ion (blue) and Total (black) Target Heat Flux\\ \\")
pyplot.xlabel("$t/\mu\mathrm{s}$")
pyplot.ylabel(r"$Q\mathrm{/W.m}^{-2}$")
pyplot.ylim(0,5e9)
pyplot.axes().xaxis.set_major_formatter(ticker.FormatStrFormatter("$%g$"))
pyplot.axes().grid(color='grey', which='both')
pyplot.tight_layout(pad=20)
pyplot.figure(8)
pyplot.plot(q_target_electron_right,'r',q_target_ion_right,'b',q_target_total_right,'k')
pyplot.rc('text', usetex=True)
pyplot.title("Electron (red), Ion (blue) and Total (black) Target Heat Flux")
pyplot.xlabel("$t/\mu\mathrm{s}$")
pyplot.ylabel(r"$Q\mathrm{/W.m}^{-2}$")
pyplot.show()
