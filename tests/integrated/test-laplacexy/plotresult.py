#!/usr/bin/env python3

from boutdata import collect
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
x2=collect('x2',path='data',xguards=False,yguards=False)
rhs2=collect('rhs2',path='data',xguards=False,yguards=False)
check2=collect('check2',path='data',xguards=False,yguards=False)

grid = DataFile("data/d3d_119919.nc")
R = grid['Rxy'][2:-2,:]
Z = grid['Zxy'][2:-2,:]

pyplot.subplot(231)
pyplot.pcolor(R,Z,x)
pyplot.colorbar()
pyplot.title('x')

pyplot.subplot(232)
pyplot.pcolor(R,Z,rhs)
pyplot.colorbar()
pyplot.title('rhs')

pyplot.subplot(233)
pyplot.pcolor(R,Z,check)
pyplot.colorbar()
pyplot.title('check')

pyplot.subplot(234)
pyplot.pcolor(R,Z,x2)
pyplot.colorbar()
pyplot.title('x2')

pyplot.subplot(235)
pyplot.pcolor(R,Z,rhs2)
pyplot.colorbar()
pyplot.title('rhs2')

pyplot.subplot(236)
pyplot.pcolor(R,Z,check2)
pyplot.colorbar()
pyplot.title('check2')

pyplot.show()

exit(0)
