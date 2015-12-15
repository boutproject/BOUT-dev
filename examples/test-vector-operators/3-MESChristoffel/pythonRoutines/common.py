#!/usr/bin/env python

"""Common files used in the python routines of this folder"""

from boutdata.collect import collect
import numpy as np

#{{{getChristoffel
def getChristoffel(alpha, zero):
    """
    Returns the analytic Christoffel symbols as a dictionary
    zero is added to the coefficients not varying with alpha in order
    to make them 2D
    """

    C = {\
         'G1_11':{'analytical':0 +zero},\
         'G1_12':{'analytical':0 +zero},\
         'G1_13':{'analytical':0 +zero},\
         'G1_22':{'analytical':0 +zero},\
         'G1_23':{'analytical':0 +zero},\
         'G1_33':{'analytical':0 +zero},\
         'G2_11':{'analytical':-2+zero},\
         'G2_12':{'analytical':0 +zero},\
         'G2_13':{'analytical':0 +zero},\
         'G2_22':{'analytical':0 +zero},\
         'G2_23':{'analytical':0 +zero},\
         'G2_33':{'analytical':0 +zero},\
         'G3_11':{'analytical':-2+zero},\
         'G3_12':{'analytical':0 +zero},\
         'G3_13':{'analytical':0 +zero},\
         'G3_22':{'analytical':0 +zero},\
         'G3_23':{'analytical':0 +zero},\
         'G3_33':{'analytical':0 +zero},\
        }
    return C
#}}}

#{{{obtainGeometry
def obtainGeometry(path, xguards, yguards):
    #{{{ Obtain the geometry
    dx  = collect('dx', path=path, yguards=yguards, xguards=xguards, info=False)
    dy  = collect('dy', path=path, yguards=yguards, xguards=xguards, info=False)
    MXG = collect('MXG', path=path, info=False)
    MYG = collect('MYG', path=path, info=False)
    dims = dx.shape
    dx  = dx[0,0]
    dy  = dy[0,0]

    #{{{ Create the y array
    if yguards:
        inner_y_points = dims[1] - 2*MYG
    else:
        inner_y_points = dims[1]
    y      = dy * np.array(np.arange(0.5, inner_y_points))
    # Insert the first and the last ghost point
    if yguards:
        y = np.insert(y,  0, - 0.5*dy)
        y = np.append(y, y[-1] + dy)
    #}}}

    #{{{Create the x array
    if xguards:
        inner_x_points = dims[0] - 2*MXG
    else:
        inner_x_points = dims[0]
    # Boundaries half between grid points
    x    = (dx * np.array(np.arange(0.5, inner_x_points)))
    # Insert the first and last ghost point
    if xguards:
        x = np.insert(x,  0, - 0.5*dx)
        x = np.append(x, x[-1] + dx)
    #}}}

    # Generate the grids
    # NOTE: If the grids are not made in these ways, the plots
    #       can look rather wierd due to mismatch of expected counting
    #       of grid points and actual counting of grid points)
    # The x, y plane
    y_YX, x_YX = np.meshgrid(y, x)
    #}}}

    #{{{Creating alpha and zero
    zero  = np.empty(x_YX.shape)
    alpha = x_YX
    #}}}

    return y_YX, x_YX, zero, alpha
#}}}
