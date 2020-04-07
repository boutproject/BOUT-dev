from __future__ import print_function
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# Find the gradient at a given point by fitting
# 
# This code is used a LOT so quite time-critical
# 
# interp_data  - Structure containing data for interpolation
# ri, zi - R and Z indices
#

import numpy  
import sys
from boututils.bunch import Bunch
from scipy.interpolate import RectBivariateSpline

def local_gradient( interp_data, ri, zi, status=0, f=None, dfdr=None, dfdz=None):
                    
    
    nr = interp_data.nx
    nz = interp_data.ny
    
    if ri < 0  or ri > nr-1 or zi < 0 or zi > nz-1 :
        status = 1
        return Bunch(status=status)
    
    if hasattr(interp_data, 'method2'):

        d = interp_data.method2
 
    else:               
                                      
        #Calculate derivatives
        
        print("Calculating derivatives for local gradient (method 2)")
        
        dd  = numpy.gradient(interp_data.f)
        ddr = dd[0]
        ddz = dd[1]
        
        xi = numpy.arange(0.,nr).astype(float)
        yi = numpy.arange(0.,nz).astype(float)
        
        fspline = RectBivariateSpline(xi, yi,interp_data.f)
        fddr = RectBivariateSpline(xi, yi, ddr)
        fddz = RectBivariateSpline(xi, yi, ddz)
        
        d = Bunch(ddr=ddr, ddz=ddz,
                  fspline = fspline,
                  fddr = fddr,
                  fddz = fddz)
        
        #print d.ddr
        #print d.ddz
        
        setattr(interp_data, 'method2', d)
        
    comp=Bunch(status=status)
    
    if f != None :
        comp.f=d.fspline(ri,zi)
    if dfdr != None :
        comp.dfdr=d.fddr(ri,zi)
    if dfdz != None :
        comp.dfdz=d.fddz(ri,zi)
        
    return comp
    
