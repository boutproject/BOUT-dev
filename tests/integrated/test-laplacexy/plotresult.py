#!/usr/bin/env python3

from boutdata import collect
from boutdata.data import BoutOptionsFile
from boututils.datafile import DataFile
from matplotlib import pyplot
import numpy
from sys import exit

# Set default colormap
pyplot.rcParams['image.cmap'] = 'viridis'
#pyplot.rcParams['image.cmap'] = 'plasma'

x=collect('x',path='data',xguards=False,yguards=False)
rhs=collect('rhs',path='data',xguards=False,yguards=False)
check=collect('check',path='data',xguards=False,yguards=False)
check_laplaceperp=collect('check_laplaceperp',path='data',xguards=False,yguards=False)
x2=collect('x2',path='data',xguards=False,yguards=False)
rhs2=collect('rhs2',path='data',xguards=False,yguards=False)
check2=collect('check2',path='data',xguards=False,yguards=False)
check2_laplaceperp=collect('check2_laplaceperp',path='data',xguards=False,yguards=False)

opts = BoutOptionsFile("data/BOUT.inp")
gridfilename = opts["grid"].replace('"','')
grid = DataFile(gridfilename)
R = grid['Rxy'][2:-2,:]
Z = grid['Zxy'][2:-2,:]

pyplot.subplot(241)
pyplot.pcolor(R,Z,x)
pyplot.colorbar()
pyplot.title('x')

pyplot.subplot(242)
pyplot.pcolor(R,Z,rhs)
pyplot.colorbar()
pyplot.title('rhs')

pyplot.subplot(243)
pyplot.pcolor(R,Z,check)
pyplot.colorbar()
pyplot.title('check')

pyplot.subplot(244)
pyplot.pcolor(R,Z,check_laplaceperp)
pyplot.colorbar()
pyplot.title('check_laplaceperp')

pyplot.subplot(245)
pyplot.pcolor(R,Z,x2)
pyplot.colorbar()
pyplot.title('x2')

pyplot.subplot(246)
pyplot.pcolor(R,Z,rhs2)
pyplot.colorbar()
pyplot.title('rhs2')

pyplot.subplot(247)
pyplot.pcolor(R,Z,check2)
pyplot.colorbar()
pyplot.title('check2')

pyplot.subplot(248)
pyplot.pcolor(R,Z,check2_laplaceperp)
pyplot.colorbar()
pyplot.title('check2_laplaceperp')

pyplot.show()

exit(0)
