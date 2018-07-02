#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
from boutdata.collect import collect
from sys import argv
from matplotlib import pyplot, ticker, rc

rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern']})



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
  print("Arguments: '[end_index] [data_path]' or 'gamma [data_path]'")
  Exit(1)

# Collect the data
Te = collect("T_electron", path=data_path, xind=2, info=True, yguards=True)
Ti = collect("T_ion", path=data_path, xind=2, info=True, yguards=True)

if end_index<0:
  end_index = len(Te[:,0,0,0])

Te_left = []
Ti_left = []
for i in range(end_index):
        Te_left.append(old_div((Te[i,0,2,0]+Te[i,0,3,0]),2))
        Ti_left.append(old_div((Ti[i,0,2,0]+Ti[i,0,3,0]),2))

# Make plot
if len(argv)>2:
  pyplot.semilogx(Te_left[:end_index],'r',Ti_left[:end_index],'b')
  pyplot.title("Te (red) and Ti (blue) at the (left) boundary")
  pyplot.axes().xaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
  pyplot.axes().grid(color='grey', which='both')
else:
  pyplot.semilogx(Te_left[:],'r',Ti_left[:],'b')
  pyplot.title("Te (red) and Ti (blue) at the (left) boundary")
  pyplot.axes().xaxis.set_major_formatter(ticker.FormatStrFormatter(r"$%g$"))
  pyplot.axes().grid(color='grey', which='both')

pyplot.show()
