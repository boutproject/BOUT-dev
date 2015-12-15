#!/usr/bin/env python

"""Vector plot in cylindrical geometry"""

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
    """Function which plots a vector field."""

    #{{{ Collect the data
    print("Showing data from " + path)
    # NOTE: We have specified that myVec has covariant components,
    # so the rho component of myvec will be myVec_x, in the contravariant
    # formalism there would not havve been a underscore between myVec and rho
    myVec_x = collect('myVec_x', xguards=False, yguards=False, path=path, info=False)
    myVec_y = collect('myVec_y', xguards=False, yguards=False, path=path, info=False)
    myVec_z = collect('myVec_z', xguards=False, yguards=False, path=path, info=False)
    #}}}

    #{{{Slice the data in time and add theta slice
    vecs = [myVec_x, myVec_y, myVec_z]

    # Slice in time
    vecs = [vec[0,:,:,:] for vec in vecs]

    # Add the last theta slice
    vecs = addSlice(vecs)
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

    # Create the theta array
    theta  = dz * np.array(range(dims[2]))

    #{{{ Create the z array
    if yguards:
        inner_y_points = dims[1] - 2*MYG
    else:
        inner_y_points = dims[1]
    z = dy * np.array(np.arange(0.5, inner_y_points))
    # Insert the first and the last ghost point
    if yguards:
        z = np.insert(z,  0, - 0.5*dy)
        z = np.append(z, z[-1] + dy)
    #}}}

    #{{{Create the rho array
    if xguards:
        inner_x_points = dims[0] - 2*MXG
    else:
        inner_x_points = dims[0]
    # Boundaries half between grid points
    rho = (dx * np.array(np.arange(0.5, inner_x_points)))
    # Insert the last ghost point (the first ghost point is already
    # present diametrically oposite of the first inner point)
    if xguards:
        rho = np.append(rho, rho[-1] + dx)
    #}}}

    # Generate the grids
    # NOTE: If the grids are not made in these ways, the plots
    #       (especially the quiver plot) can look rather wierd due to
    #       mismatch of expected counting of grid points and actual
    #       counting of grid points)
    # The rho, theta plane
    THETA_TR, RHO_TR = np.meshgrid(theta, rho)
    # The rho, z plane
    Z_ZR, RHO_ZR = np.meshgrid(z, rho)
    # The z, theta plane
    THETA_TZ, Z_TZ = np.meshgrid(theta, z)
    #}}}

    #{{{ Geometry mapping
    X_TR = RHO_TR*np.cos(THETA_TR)
    Y_TR = RHO_TR*np.sin(THETA_TR)
    #}}}

    #{{{ Create the axes
    # Figure for the vector plots
    f1 = plt.figure(figsize = plt_size)
    ncols = 2
    gs  = gridspec.GridSpec(nrows = 1, ncols = ncols)
    f1ax1 = plt.subplot(gs[0])
    f1ax2 = plt.subplot(gs[1])
    all_ax = [f1ax1, f1ax2]
    #}}}

    #{{{Slice the data in space
    vecRho, vecZ, vecTheta = vecs
    # Make the cuts for the RhoTheta plots (slice in Z)
    zInd = int(vecRho.shape[1]/2.0)
    vecRhoRT   = vecRho  [:, zInd, :]
    vecThetaRT = vecTheta[:, zInd, :]
    # Make the cuts for the RhoZ plots (slice in Theta)
    thetaInd   = 0
    vecRhoRZ   = vecRho  [:, :, thetaInd]
    vecZRZ     = vecZ    [:, :, thetaInd]
    vecThetaRZ = vecTheta[:, :, thetaInd]
    # Negative side of rho-theta plot
    # Find the opposite side
    piInd = round(vecRho.shape[2]/2.0)
    if thetaInd > piInd:
        vecRhoRZNeg   = vecRho  [:, :, thetaInd - piInd]
        vecZRZNeg     = vecZ    [:, :, thetaInd - piInd]
        vecThetaRZNeg = vecTheta[:, :, thetaInd - piInd]
    else:
        vecRhoRZNeg   = vecRho  [:, :, thetaInd + piInd]
        vecZRZNeg     = vecZ    [:, :, thetaInd + piInd]
        vecThetaRZNeg = vecTheta[:, :, thetaInd + piInd]
    #}}}

    #{{{ Vector mappings
    # In the rho-z plane, rho will map to x naturally
    # However, in the rho-theta plot, rho should not map to x, as the
    # direction of rho is changing from point to point. Thus, we map
    # vecRhoRT and vecThetaRT to vecX and vecY
    vecX = vecRhoRT*(X_TR/np.sqrt(X_TR**2 + Y_TR**2)) - vecThetaRT*(Y_TR/(X_TR**2 + Y_TR**2))
    vecY = vecRhoRT*(Y_TR/np.sqrt(X_TR**2 + Y_TR**2)) + vecThetaRT*(X_TR/(X_TR**2 + Y_TR**2))
    # vecRhoRZNeg and vecZRZNeg are rotated, the rotation matrix should
    # be the same for both co and contravariant vectors
    # If we multiply the rotation matrix
    # https://en.wikipedia.org/wiki/Rotation_matrix
    # to the vector components, we get
    vecRhoRZNeg = -vecRhoRZNeg
    #}}}

    #{{{ Get the limits to set in the rho-z plot
    #{{{ For the vector plot
    # The original rho-z fields
    # Calculate the length
    lenRzVecs = np.sqrt(vecRhoRZ**2 + vecZRZ**2)
    # Find the max and min
    posMax = np.max(lenRzVecs)
    posMin = np.min(lenRzVecs)

    # rho-z fields on the negative side
    # Calculate the length
    lenRzVecsNeg = np.sqrt(vecRhoRZNeg**2 + vecZRZNeg**2)
    # Find the max and min
    negMax = np.max(lenRzVecsNeg)
    negMin = np.min(lenRzVecsNeg)

    # Find vmax and vmin
    vmaxRzVec = posMax if posMax > negMax else negMax
    vminRzVec = posMin if posMin < negMin else negMin
    #}}}

    #{{{ For the function plot
    # The original rho-z fields
    # Find the max and min
    posMax = np.max(vecThetaRZ)
    posMin = np.min(vecThetaRZ)

    # rho-z fields on the negative side
    # Find the max and min
    negMax = np.max(vecThetaRZNeg)
    negMin = np.min(vecThetaRZNeg)

    # Find vmax and vmin
    vmaxRzFunc = posMax if posMax > negMax else negMax
    vminRzFunc = posMin if posMin < negMin else negMin
    #}}}
    #}}}

    #{{{ Plot
    #{{{ Vector plots
    # The rho, rho plot
    f1ax1.contourf(X_TR, Y_TR, np.sqrt(vecX**2 + vecY**2),\
                 nCount, cmap='RdYlBu_r')
    f1ax1.quiver(X_TR, Y_TR, vecX, vecY, scale=scale)

    # The rho-z plot
    f1ax2.contourf(RHO_ZR, Z_ZR, lenRzVecs,\
                 nCount, vmax = vmaxRzVec, vmin=vminRzVec, cmap='RdYlBu_r')
    f1ax2.quiver(RHO_ZR, Z_ZR, vecRhoRZ, vecZRZ, scale=scale)

    f1ax2.contourf(-RHO_ZR, Z_ZR, lenRzVecsNeg,\
                 nCount, vmax = vmaxRzVec, vmin=vminRzVec, cmap='RdYlBu_r')
    f1ax2.quiver(-RHO_ZR, Z_ZR, vecRhoRZNeg, vecZRZNeg, scale=scale)
    #}}}

    #{{{ Decorations
    ax1xl = r'$\rho$'
    ax1yl = r'$\rho$'
    ax1t =  r'$z$'+'-slice'

    ax2xl = r'$\rho$'
    ax2yl = r'$z$'
    ax2t =  r'$\theta$'+'-slice'

    f1ax1.set_xlabel(ax1xl)
    f1ax1.set_ylabel(ax1yl)
    f1ax1.set_title (ax1t)

    f1ax2.set_xlabel(ax2xl)
    f1ax2.set_ylabel(ax2yl)
    f1ax2.set_title (ax2t)

    for ax in all_ax:
        ax.grid()
        ax.axis('equal')
    #}}}

    plt.tight_layout()
    plt.show()
    #}}}
#}}}

#{{{addSlice
def addSlice(vecs):
    """Adds the values in theta=0 in the new point theta=2*pi"""

    # Determine sizes
    oldSize     = np.array(vecs[0].shape)
    newSize     = oldSize.copy()
    newSize[-1] = oldSize[-1]+1

    # Appendable list
    newVecs = []
    for vec in vecs:
        # Create the new field
        newField = np.empty(newSize)

        # Fill the new field with the old data
        # (NOTE: End-point index does not start counting on 0)
        newField[ :, :, :oldSize[-1]] = vec
        # And add the values from point theta = 0 to theta = 2*pi
        # (-1 since start count on 0)
        newField[:, :, newSize[-1]-1] = vec[:,:,0]
        newVecs.append(newField)

    return newVecs
#}}}
