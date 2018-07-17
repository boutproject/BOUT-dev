#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Runs the conduction example, produces some output
#

from __future__ import division
from builtins import str, range
from past.utils import old_div

from boutdata.collect import collect
from sys import argv
from math import log, pi
from matplotlib import pyplot

nproc = 1  # Number of processors to use

gamma = 3.
if len(argv)>1:
  data_path = str(argv[1])
else:
  data_path = "data"

electron_mass = 9.10938291e-31
ion_mass =  3.34358348e-27

# Collect the data
Te = collect("T_electron", path=data_path, info=True, yguards=True)
Ti = collect("T_ion", path=data_path, info=True, yguards=True)
n = collect("n_ion", path=data_path, info=True, yguards=True)
V = collect("Vpar_ion", path=data_path, info=True, yguards=True)
q = collect("heat_flux", path=data_path, info=True, yguards=True)

q_electron_left = []
q_electron_right = []
right_index = len(Te[0,2,:,0])-4

for i in range(len(Te[:,2,0,0])):
  Te_left = old_div((Te[i,2,2,0]+Te[i,2,1,0]),2.)
  Ti_left = old_div((Ti[i,2,2,0]+Ti[i,2,1,0]),2.)
  n_left = old_div((n[i,2,2,0]+n[i,2,1,0]),2.)
  Te_right = old_div((Te[i,2,right_index,0]+Te[i,2,right_index+1,0]),2)
  Ti_right = old_div((Ti[i,2,right_index,0]+Ti[i,2,right_index+1,0]),2)
  n_right = old_div((n[i,2,right_index,0]+n[i,2,right_index+1,0]),2)
  sheath_potential = 0.5*Te_left*log(2*pi*electron_mass/ion_mass*(1+gamma*Ti_left/Te_left))
  q_electron_left.append((2.0*Te_left-sheath_potential)*n_left*V[i,2,2,0]) # in W/m^2

  sheath_potential = 0.5*Te_right*log(2*pi*electron_mass/ion_mass*(1+gamma*Ti_right/Te_right))
  q_electron_right.append((2.0*Te_right-sheath_potential)*n_right*V[i,2,right_index+1,0]) # in W/m^2

pyplot.figure(1)
pyplot.plot(q_electron_left,'r',q[:,2,2,0],'b',q_electron_right,'r',q[:,2,right_index+1,0],'b')
pyplot.title("Electron heat flux at the boundaries (blue) and calculated boundary value (red)\n\n")
pyplot.xlabel(u"t/μs")
pyplot.ylabel("Q/eV.m$^{-2}$")
pyplot.figure(2)
pyplot.plot(q[:,2,2,0]-q_electron_left,'b',q[:,2,right_index+1,0]-q_electron_right,'r')
pyplot.title("Difference between heat flux and its calculated boundary value at the left (blue) and right (red) boundaries\n\n")
pyplot.xlabel(u"t/μs")
pyplot.ylabel("dQ/eV.m$^{-2}$")
pyplot.show()
