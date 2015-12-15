#!/usr/bin/env python

"""Shows the Christoffel symbols"""

import matplotlib.pyplot as plt
from matplotlib import gridspec, get_backend
from pythonRoutines.common import getChristoffel, obtainGeometry
from boutdata.collect import collect
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
plt.rc("axes.formatter", useoffset = False)
plt_size = (18, 12)
#}}}

# All post processing functions called by bout_runners must accept the
# first argument from bout_runners (called 'folder' in
# __call_post_processing_function)
#{{{plotChristoffel
def plotChristoffel(path,\
                    showRelErr = True, showBOUT = False, showAnalytical = False,\
                    xguards = False, yguards = False, nCount = 100,
                    showOnlyThese = False):
    """
    Check the difference between the analytical calculated and BOUT++
    Christoffel symbols

    showRelErr     - show the relative error
    showBOUT       - show the Christoffel symbols calculated in BOUT++
    showAnalytical - show analytically calculated Christoffel symbols
    showOnlyThese  - array of Christoffel symbols going to be shown
    """
    # Set numpy to raise all errors
    np.seterr(all='raise')

    print("Showing data from " + path)
    # Get the geometry
    y_YX, x_YX, zero, alpha = obtainGeometry(path, xguards, yguards)
    # Get the Christoffel symbols
    C = getChristoffel(alpha, zero)

    if showOnlyThese != False:
        Cnew = {}
        for coeff in showOnlyThese:
            Cnew[coeff] = C.pop(coeff)
        C = Cnew


    #{{{ Loop over the coefficients
    for coeff in C.keys():
        #{{{ Collect the data
        C[coeff]['bout'] =\
                collect(coeff, xguards=xguards, yguards=yguards, path=path, info=False)
        # Pick the zeroth z-point (the Field2D quantities are isotropic in
        # BOUT++)
        C[coeff]['bout'] = C[coeff]['bout'][:,:,0]
        try:
            C[coeff]['relErr'] =\
                    np.abs(C[coeff]['analytical'] - C[coeff]['bout'])\
                           /np.abs(C[coeff]['analytical'])
        except FloatingPointError:
            C[coeff]['relErr'] = zero
        #}}}

        #{{{ Make the plots
        if showRelErr:
            C[coeff]['relErrFig'] = {'fig':plt.figure(figsize = plt_size)}
            C[coeff]['relErrFig']['ax'] = plt.subplot()
            C[coeff]['relErrFig']['plot'] =\
                    C[coeff]['relErrFig']['ax'].contourf(x_YX, y_YX,\
                                                      C[coeff]['relErr'], nCount,\
                                                      cmap='RdYlBu_r')
            C[coeff]['relErrFig']['ax'].set_xlabel('x')
            C[coeff]['relErrFig']['ax'].set_ylabel('y')
            C[coeff]['relErrFig']['cbar'] = C[coeff]['relErrFig']['fig'].colorbar(C[coeff]['relErrFig']['plot'])
            C[coeff]['relErrFig']['cbar'].set_label(\
                    r'$\frac{|\mathrm{analytical} - \mathrm{BOUT++}|}'+\
                    r'{|\mathrm{analytical}|}$', size = 50)
            C[coeff]['relErrFig']['fig'].canvas.set_window_title(coeff)
        if showBOUT:
            C[coeff]['boutFig'] = {'fig':plt.figure(figsize = plt_size)}
            C[coeff]['boutFig']['ax'] = plt.subplot()
            C[coeff]['boutFig']['plot'] =\
                    C[coeff]['boutFig']['ax'].contourf(x_YX, y_YX,\
                                                       C[coeff]['bout'], nCount,\
                                                       cmap='RdYlBu_r')
            C[coeff]['boutFig']['ax'].set_xlabel('x')
            C[coeff]['boutFig']['ax'].set_ylabel('y')
            C[coeff]['boutFig']['cbar'] = C[coeff]['boutFig']['fig'].colorbar(C[coeff]['boutFig']['plot'])
            C[coeff]['boutFig']['cbar'].set_label('BOUT++')
            C[coeff]['boutFig']['fig'].canvas.set_window_title(coeff)
        if showAnalytical:
            C[coeff]['analyticalFig'] = {'fig':plt.figure(figsize = plt_size)}
            C[coeff]['analyticalFig']['ax'] = plt.subplot()
            C[coeff]['analyticalFig']['plot'] =\
                    C[coeff]['analyticalFig']['ax'].contourf(x_YX, y_YX,\
                                                      C[coeff]['analytical'], nCount,\
                                                      cmap='RdYlBu_r')
            C[coeff]['analyticalFig']['ax'].set_xlabel('x')
            C[coeff]['analyticalFig']['ax'].set_ylabel('y')
            C[coeff]['analyticalFig']['cbar'] = C[coeff]['analyticalFig']['fig'].colorbar(C[coeff]['analyticalFig']['plot'])
            C[coeff]['analyticalFig']['cbar'].set_label('Analytical')
            C[coeff]['analyticalFig']['fig'].canvas.set_window_title(coeff)
        #}}}
    #}}}

    plt.show()
#}}}
