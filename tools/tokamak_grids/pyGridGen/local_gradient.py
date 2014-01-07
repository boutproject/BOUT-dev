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
from bunch import Bunch
from scipy.interpolate import RectBivariateSpline

def local_gradient( interp_data, ri, zi, status=0, f=None, dfdr=None, dfdz=None):
                    
    
    nr = interp_data.nx
    nz = interp_data.ny
    xi=numpy.arange(0.,nr).astype(float)
    yi=numpy.arange(0.,nz).astype(float)
    
      
  
    if ri < 0  or ri > nr-1 or zi < 0 or zi > nz-1 :
        status = 1
        return Bunch(status=status)
             
    if hasattr(interp_data, 'method2'):

        d = interp_data.method2
 
    else:               
                                      
    #Calculate derivatives
        
        print "Calculating derivatives for local gradient (method 2)"

        dd = numpy.gradient(interp_data.f)
        ddr=dd[0]
        ddz=dd[1]
       
        d = Bunch(ddr=ddr, ddz=ddz)
        
        #print d.ddr
        #print d.ddz
        
        setattr(interp_data, 'method2', d)
        
          
    
    comp=Bunch(status=status)
    
#    print ri, zi
    if f != None :
        func = RectBivariateSpline(xi, yi,interp_data.f)
        f=func(ri,zi)
        comp.f=f
      #  print 'f',f
    if dfdr != None :
         
        func = RectBivariateSpline(xi, yi, d.ddr)
        dfdr = func(ri,zi)
        comp.dfdr=dfdr
        #print 'dfdr',dfdr
    if dfdz != None :
       
        func = RectBivariateSpline(xi, yi,  d.ddz)
        dfdz = func(ri,zi)
        comp.dfdz=dfdz
        #print 'dfdz',dfdz
    
    

    return comp
    
