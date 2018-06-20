#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division
from builtins import str, range
from past.utils import old_div

from boutdata.collect import collect

from sys import argv
from matplotlib import pyplot, ticker, rc

rc('text', usetex=True)
#rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('font',**{'family':'serif','serif':['Computer Modern'],'size':20})


if len(argv)==1:
  end_index = -1
  data_path = "data"
elif len(argv)==2:
  try:
    end_index = int(argv[1])
    data_path = "data"
  except ValueError:
    end_index = -1
    data_path = str(argv[1])
elif len(argv)==3:
  end_index = int(argv[1])
  data_path = str(argv[2])
else:
  print("Arguments: '[end_index] [data_path]' or '[data_path]'")
  Exit(1)

electron_mass = 9.10938291e-31
ion_mass =  3.34358348e-27

# Collect the data
f = collect("effective_alpha_e", path=data_path, info=True)
#Te = collect("T_electron", path=data_path, xind=2, info=True, yguards=True)
#Ti = collect("T_ion", path=data_path, xind=2, info=True, yguards=True)
#n = collect("n_ion", path="data", xind=2, info=True, yguards=True)

if end_index<0:
  end_index = len(f[:])

alpha0 = 1
alphaoveralpha0 = []
for i in range(end_index):
        alphaoveralpha0.append(old_div(f[i],alpha0))
        print(i,alphaoveralpha0[i])

# Make plot
pyplot.figure(1, facecolor='w')
#pyplot.figure(1,dpi=800, facecolor='w')

if len(argv)>2:
#  pyplot.semilogy(f[:end_index],'k')
  pyplot.semilogy(alphaoveralpha0[:end_index],'k')
else:
#  pyplot.semilogy(f[:],'k')
  pyplot.semilogy(alphaoveralpha0[:],'k')

#pyplot.title("Effective electron heat-flux limiter")
pyplot.axes().xaxis.set_major_formatter(ticker.FormatStrFormatter("$%g$"))
pyplot.axes().grid(color='grey', which='both')
#pyplot.xlabel(u"t/μs")
pyplot.xlabel("$t/\mu\mathrm{s}$")
#pyplot.ylabel(u"α_e")
pyplot.ylabel("$\langle \\alpha_e \\rangle$")
#pyplot.ylabel("$\langle \\alpha_{\\mathrm{e}} \\rangle$")
pyplot.tight_layout(pad=.1)

pyplot.savefig('fig.pdf')
pyplot.savefig('fig.eps',dpi=1200)

pyplot.show()
