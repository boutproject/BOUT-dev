#!/usr/bin/env python3

from boutdata import collect
from matplotlib import pyplot
from sys import exit

f = collect('f', path='data', yguards=True, info=False)[1:-1,1:-1]
sol = collect('sol', path='data', yguards=True, info=False)[1:-1,1:-1]
error = collect('error', path='data', yguards=True, info=False)[1:-1,1:-1]
absolute_error = collect('absolute_error', path='data', yguards=True, info=False)[1:-1,1:-1]

# Note, cells closest to x-boundary in rhs and rhs_check may be slightly different because
# of different way x-boundary cells for D2DXDY are set: in LaplaceXY corner guard cells
# are set so that a 9-point stencil can be used; in the D2DXDY function (used in
# Laplace_perp) first dfdy=DDY(f) is calculated, communicated and has free_o3 x-boundary
# conditions applied, then DDX(dfdy) is returned.
# Therefore here exclude cells closest to the x-boundary so that the difference plotted
# should be small (as controlled by rtol, atol).
rhs = collect('rhs', path='data', yguards=True, info=False)[3:-3,2:-2]
rhs_check = collect('rhs_check', path='data', yguards=True, info=False)[3:-3,2:-2]

pyplot.figure()

pyplot.subplot(231)
pyplot.pcolormesh(f)
pyplot.title('f')
pyplot.colorbar()

pyplot.subplot(232)
pyplot.pcolormesh(sol)
pyplot.title('sol')
pyplot.colorbar()

pyplot.subplot(233)
pyplot.pcolormesh(error)
pyplot.title('error')
pyplot.colorbar()

pyplot.subplot(234)
pyplot.pcolormesh(absolute_error)
pyplot.title('absolute_error')
pyplot.colorbar()

pyplot.subplot(235)
pyplot.pcolormesh(rhs)
pyplot.title('rhs')
pyplot.colorbar()

pyplot.subplot(236)
pyplot.pcolormesh(rhs - rhs_check)
pyplot.title('rhs diff')
pyplot.colorbar()

pyplot.show()

exit(0)
