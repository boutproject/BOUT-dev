"""Average over a surface

"""

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import range
from past.utils import old_div
import numpy as np
from boututils.calculus import deriv
from boututils.int_func import int_func
from .idl_tabulate import idl_tabulate


def surface_average(var, grid, area=None):
    """Average a variable over a surface

    Parameters
    ----------
    var : array_like
        3-D or 4D variable to integrate (either [x,y,z] or [t,x,y,z])
    grid : dict
        A dictionary of various grid quantities
    area : bool
        Average by flux-surface area = (B/Bp)*dl * R*dz

    Returns
    -------
    float
        Surface average of variable

    """

    s = np.ndim(var)

    if s == 4 :
        nx = np.shape(var)[1]
        ny = np.shape(var)[2]
        nt = np.shape(var)[0]

        result = np.zeros((nx,nt))
        for t in range (nt):
            result[:,t] = surface_average(var[t,:,:,:], grid, area=area)

        return result
    elif s != 3 :
        raise RuntimeError("ERROR: surface_average var must be 3D or 4D")

    # 3D [x,y,z]
    nx = np.shape(var)[0]
    ny = np.shape(var)[1]


    # Calculate poloidal angle from grid
    theta = np.zeros((nx,ny))

    #status = gen_surface(mesh=grid) ; Start generator
    xi = -1
    yi = np.arange(0,ny,dtype=int)
    last = 0
    while True:
    #yi = gen_surface(last=last, xi=xi, period=periodic)
        xi = xi + 1
        if xi == nx-1 :
            last = 1

        dtheta = 2.*np.pi / np.float(ny)
        r = grid['Rxy'][xi,yi]
        z = grid['Zxy'][xi,yi]
        n = np.size(r)

        dl = old_div(np.sqrt( deriv(r)**2 + deriv(z)**2 ), dtheta)
        if area:
            dA = (old_div(grid['Bxy'][xi,yi],grid['Bpxy'][xi,yi]))*r*dl
            A = int_func(np.arange(n),dA)
            theta[xi,yi] = 2.*np.pi*A/A[n-1]
        else:
            nu = dl * (grid['Btxy'][xi,yi]) / ((grid['Bpxy'][xi,yi]) * r )
            theta[xi,yi] = int_func(np.arange(n)*dtheta,nu)
            theta[xi,yi] = 2.*np.pi*theta[xi,yi] / theta[xi,yi[n-1]]

        if last==1 : break

    vy = np.zeros(ny)
    result = np.zeros(nx)
    for x in range(nx) :
        for y in range(ny) :
            vy[y] = np.mean(var[x,y,:])

        result[x] = old_div(idl_tabulate(theta[x,:], vy), (2.*np.pi))

    return result
