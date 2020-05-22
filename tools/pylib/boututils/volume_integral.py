"""Integrate over a volume

"""

from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
import numpy as np
from boututils.calculus import deriv

def volume_integral(var, grid, xr=False):
    """Integrate a variable over a volume

    Parameters
    ----------
    var : array_like
        Variable to integrate
    grid : dict
        A dictionary of various grid quantities
    xr : (int, int), optional
        Range of x indices (default: all of x)

    Returns
    -------
    float
        Volumne integral of variable

    """

    s = np.ndim(var)

    if s == 4 :
        # 4D [t,x,y,z] - integrate for each t
        nx = np.shape(var)[1]
        ny = np.shape(var)[2]
        nt = np.shape(var)[0]

        result = np.zeros(nt)
        for t in range(nt) :
            result[t] = volume_integral(var[t,:,:,:],g,xr=xr)
        return result

    elif s == 3 :
        # 3D [x,y,z] - average in Z
        nx = np.shape(var)[0]
        ny = np.shape(var)[1]
 #       nz = np.shape(var)[2]

        zi = np.zeros((nx, ny))
        for x in range(nx):
            for y in range(ny):
                zi[x,y] = np.mean(var[x,y,:])

        return volume_integral(zi, g, xr=xr)


    elif s != 2 :
        print("ERROR: volume_integral var must be 2, 3 or 4D")


    # 2D [x,y]
    nx = np.shape(var)[0]
    ny = np.shape(var)[1]

    if xr == False : xr=[0,nx-1]

    result = 0.0

    #status = gen_surface(mesh=grid) ; Start generator
    xi = -1
    yi = np.arange(0,ny,dtype=int)
    last = 0
  #  iy = np.zeros(nx)
    while True:
        #yi = gen_surface(last=last, xi=xi, period=periodic)
        xi = xi + 1
        if xi == nx-1 : last = 1

        if (xi >= np.min(xr)) & (xi <= np.max(xr)) :
            dtheta = 2.*np.pi / np.float(ny)
            r = grid['Rxy'][xi,yi]
            z = grid['Zxy'][xi,yi]
            n = np.size(r)
            dl = old_div(np.sqrt( deriv(r)**2 + deriv(z)**2 ), dtheta)

      # Area of flux-surface
            dA = (grid['Bxy'][xi,yi]/grid['Bpxy'][xi,yi]*dl) * (r*2.*np.pi)
      # Volume
            if xi == nx-1 :
                dpsi = (grid['psixy'][xi,yi] - grid['psixy'][xi-1,yi])
            else:
                dpsi = (grid['psixy'][xi+1,yi] - grid['psixy'][xi,yi])

            dV = dA * dpsi / (r*(grid['Bpxy'][xi,yi])) # May need factor of 2pi
            dV = np.abs(dV)

            result = result + np.sum(var[xi,yi] * dV)

        if last==1 : break

    return result
