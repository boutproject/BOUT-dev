from __future__ import division
from past.utils import old_div
import numpy
from gen_surface import gen_surface
import copy
from boututils.fft_deriv import fft_deriv

# Take derivative in y, taking into account branch-cuts
def DDY( var, mesh):
  f = copy.deepcopy(var)

  dtheta = 2.*numpy.pi / numpy.float(numpy.sum(mesh.npol))

  status = gen_surface(mesh=mesh) # Start generator
  while True:
    period, yi, xi, last = gen_surface(last=None, xi=None, period=None)
    if period :
        f[xi,yi] = numpy.real(fft_deriv(var[xi,yi]))
    else:
        f[xi,yi] = numpy.gradient(var[xi,yi])
    if last : break
  return old_div(f, dtheta)
