#!/usr/bin/env python

"""Vector plot in cartesian geometry"""

from boutdata.collect import collect
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pylab import plot
import numpy as np

#{{{ Set the plot style
title_size = 30
plt.rc("font", size = 30)
plt.rc("axes", labelsize = 25, titlesize = title_size)
plt.rc("xtick", labelsize = 25)
plt.rc("ytick", labelsize = 25)
plt.rc("legend", fontsize = 30)
plt.rc("lines", linewidth = 2)
plt_size = (18, 12)
#}}}

# All post processing functions called by bout_runners must accept the
# first argument from bout_runners (called 'folder' in
# __call_post_processing_function)
#{{{plotVecField
def plotVecField(path,\
                 xguards = False,\
                 yguards = False,\
                 scale = None,\
                 nCount = 100):
    """Function which plots the vector field."""

    #{{{ Collect the data
    print("Showing data from " + path)
    # NOTE: We have specified that myVec_ has contravariant components,
    # so the x component of myvec will be myVec_x, in the covariant
    # formalism there would've been a underscore between myVec_ and x
    myVec_x = collect('myVec_x', xguards=False, yguards=False, path=path, info=False)
    myVec_y = collect('myVec_y', xguards=False, yguards=False, path=path, info=False)
    myVec_z = collect('myVec_z', xguards=False, yguards=False, path=path, info=False)
    #}}}

    #{{{Slice the data in time
    vecs = [myVec_x, myVec_y, myVec_z]

    # Slice in time
    vecs = [vec[0,:,:,:] for vec in vecs]
    #}}}

    #{{{ Obtain the geometry
    dx  = collect('dx', path=path, yguards=yguards, xguards=xguards, info=False)
    dy  = collect('dy', path=path, yguards=yguards, xguards=xguards, info=False)
    dz  = collect('dz',  path=path, info=False)
    MXG = collect('MXG', path=path, info=False)
    MYG = collect('MYG', path=path, info=False)
    dx  = dx[0,0]
    dy  = dy[0,0]

    # Obtaining the dimension of the data
    dims=vecs[0].shape

    # Create the z array
    z  = dz * np.array(range(dims[2]))

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
    # Insert the last ghost point (the first ghost point is already
    # present diametrically oposite of the first inner point)
    if xguards:
        x = np.append(x, x[-1] + dx)
    #}}}

    # Generate the grids
    # NOTE: If the grids are not made in these ways, the plots
    #       (especially the quiver plot) can look rather wierd due to
    #       mismatch of expected counting of grid points and actual
    #       counting of grid points)
    # The x, z plane
    z_ZX, x_ZX = np.meshgrid(z, x)
    # The x, y plane
    y_YX, x_YX = np.meshgrid(y, x)
    # The y, z plane
    z_ZY, y_ZY = np.meshgrid(z, y)
    #}}}

    #{{{ Create the axes
    # Figure for the vector plots
    f1 = plt.figure(figsize = plt_size)
    ncols = 3
    gs  = gridspec.GridSpec(nrows = 1, ncols = ncols)
    f1ax1 = plt.subplot(gs[0])
    f1ax2 = plt.subplot(gs[1])
    f1ax3 = plt.subplot(gs[2])
    all_ax = [f1ax1, f1ax2, f1ax3]
    #}}}

    #{{{Slice the data in space
    # Make the cuts
    vecsXCut = [vec[int(vec.shape[0]/2.0), :, :] for vec in vecs]
    vecsYCut = [vec[:,int(vec.shape[1]/2.0),:] for vec in vecs]
    vecsZCut = [vec[:,:,int(vec.shape[2]/2.0)] for vec in vecs]
    #}}}

    #{{{ Plot
    f1ax1.contourf(x_YX, y_YX, np.sqrt(vecsZCut[0]**2 + vecsZCut[1]**2),\
                 nCount, cmap='RdYlBu_r')
    f1ax1.quiver(x_YX, y_YX, vecsZCut[0], vecsZCut[1], scale=scale)

    f1ax2.contourf(x_ZX, z_ZX, np.sqrt(vecsYCut[0]**2 + vecsYCut[2]**2),\
                 nCount, cmap='RdYlBu_r')
    f1ax2.quiver(x_ZX, z_ZX, vecsYCut[0], vecsYCut[2], scale=scale)

    f1ax3.contourf(y_ZY, z_ZY, np.sqrt(vecsXCut[1]**2 + vecsXCut[2]**2),\
                 nCount, cmap='RdYlBu_r')
    f1ax3.quiver(y_ZY, z_ZY, vecsXCut[1], vecsXCut[2], scale=scale)

    #{{{ Decorations
    ax1xl = 'x'
    ax1yl = 'y'
    ax1t =  'z-slice'

    ax2xl = 'x'
    ax2yl = 'z'
    ax2t =  'y-slice'

    ax3xl = 'y'
    ax3yl = 'z'
    ax3t =  'x-slice'

    f1ax1.set_xlabel(ax1xl)
    f1ax1.set_ylabel(ax1yl)
    f1ax1.set_title (ax1t)

    f1ax2.set_xlabel(ax2xl)
    f1ax2.set_ylabel(ax2yl)
    f1ax2.set_title (ax2t)

    f1ax3.set_xlabel(ax3xl)
    f1ax3.set_ylabel(ax3yl)
    f1ax3.set_title (ax3t)

    for ax in all_ax:
        ax.grid()
        ax.axis('equal')
    #}}}

    plt.tight_layout()
    plt.show()
    #}}}
#}}}
