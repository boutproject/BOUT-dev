"""Equilibrium analysis routine

Takes a RZ psi grid, and finds x-points and o-points
"""

from __future__ import print_function
from __future__ import division

from builtins import zip
from builtins import str
from builtins import range
from past.utils import old_div

import numpy
from . import local_min_max
from scipy.interpolate import RectBivariateSpline
from matplotlib.pyplot import contour, gradient,  annotate, plot, draw
from crosslines import find_inter


def analyse_equil(F, R, Z):
    """Takes an RZ psi grid, and finds x-points and o-points

    Parameters
    ----------
    F : array_like
        2-D array of psi values
    R : array_like
        1-D array of major radii, its length should be the same as the
        first dimension of F
    Z : array_like
        1-D array of heights, its length should be the same as the
        second dimension of F

    Returns
    -------
    object
        An object of critical points containing:

          n_opoint, n_xpoint   - Number of O- and X-points
          primary_opt          - Index of plasma centre O-point
          inner_sep            - X-point index of inner separatrix
          opt_ri, opt_zi       - R and Z indices for each O-point
          opt_f                - Psi value at each O-point
          xpt_ri, xpt_zi       - R and Z indices for each X-point
          xpt_f                - Psi value of each X-point

    """
    s = numpy.shape(F)
    nx = s[0]
    ny = s[1]
  
  #;;;;;;;;;;;;;;; Find critical points ;;;;;;;;;;;;;
  #
  # Need to find starting locations for O-points (minima/maxima)
  # and X-points (saddle points)
  #
    Rr=numpy.tile(R,nx).reshape(nx,ny).T
    Zz=numpy.tile(Z,ny).reshape(nx,ny)
    
    contour1=contour(Rr,Zz,gradient(F)[0], levels=[0.0], colors='r')
    contour2=contour(Rr,Zz,gradient(F)[1], levels=[0.0], colors='r')    

    draw()


### 1st method - line crossings ---------------------------
    res=find_inter( contour1, contour2)
            
    #rex1=numpy.interp(res[0],  R, numpy.arange(R.size)).astype(int)
    #zex1=numpy.interp(res[1],  Z, numpy.arange(Z.size)).astype(int)

    rex1=res[0]
    zex1=res[1]    
            
    w=numpy.where((rex1 > R[2]) & (rex1 < R[nx-3]) & (zex1 > Z[2]) & (zex1 < Z[nx-3]))
    nextrema = numpy.size(w)
    rex1=rex1[w].flatten()
    zex1=zex1[w].flatten()
    

### 2nd method - local maxima_minima -----------------------
    res1=local_min_max.detect_local_minima(F)
    res2=local_min_max.detect_local_maxima(F)
    res=numpy.append(res1,res2,1)

    rex2=res[0,:].flatten()
    zex2=res[1,:].flatten()
     
      
    w=numpy.where((rex2 > 2) & (rex2 < nx-3) & (zex2 >2) & (zex2 < nx-3))
    nextrema = numpy.size(w)
    rex2=rex2[w].flatten()
    zex2=zex2[w].flatten()
    

    n_opoint=nextrema
    n_xpoint=numpy.size(rex1)-n_opoint
    
 # Needed for interp below   
    
    Rx=numpy.arange(numpy.size(R))
    Zx=numpy.arange(numpy.size(Z))
                                                               

                                                                                          
    print("Number of O-points: "+numpy.str(n_opoint))
    print("Number of X-points: "+numpy.str(n_xpoint))
    
 # Deduce the O & X points   
    
    x=R[rex2]
    y=Z[zex2]
    
    dr=old_div((R[numpy.size(R)-1]-R[0]),numpy.size(R))
    dz=old_div((Z[numpy.size(Z)-1]-Z[0]),numpy.size(Z))


    repeated=set()
    for i in range(numpy.size(rex1)):
        for j in range(numpy.size(x)):
            if numpy.abs(rex1[i]-x[j]) < 2*dr and numpy.abs(zex1[i]-y[j]) < 2*dz : repeated.add(i)
        
 # o-points
 
    o_ri=numpy.take(rex1,numpy.array(list(repeated)))
    opt_ri=numpy.interp(o_ri,R,Rx)
    o_zi=numpy.take(zex1,numpy.array(list(repeated)))  
    opt_zi=numpy.interp(o_zi,Z,Zx)      
    opt_f=numpy.zeros(numpy.size(opt_ri))
    func = RectBivariateSpline(Rx, Zx, F)
    for i in range(numpy.size(opt_ri)): opt_f[i]=func(opt_ri[i], opt_zi[i])

    n_opoint=numpy.size(opt_ri)
    
 # x-points
                       
    x_ri=numpy.delete(rex1, numpy.array(list(repeated)))
    xpt_ri=numpy.interp(x_ri,R,Rx)
    x_zi=numpy.delete(zex1, numpy.array(list(repeated)))
    xpt_zi=numpy.interp(x_zi,Z,Zx)
    xpt_f=numpy.zeros(numpy.size(xpt_ri))
    func = RectBivariateSpline(Rx, Zx, F)
    for i in range(numpy.size(xpt_ri)): xpt_f[i]=func(xpt_ri[i], xpt_zi[i])
    
    n_xpoint=numpy.size(xpt_ri)
    
 #  plot o-points 

    plot(o_ri,o_zi,'o', markersize=10)
     
    labels = ['{0}'.format(i) for i in range(o_ri.size)]
    for label, xp, yp in zip(labels, o_ri, o_zi):
        annotate(label,  xy = (xp, yp), xytext = (10, 10), textcoords = 'offset points',size='large', color='b')

    draw()
  
 #  plot x-points     
    
    plot(x_ri,x_zi,'x', markersize=10)
     
    labels = ['{0}'.format(i) for i in range(x_ri.size)]
    for label, xp, yp in zip(labels, x_ri, x_zi):
        annotate(label,  xy = (xp, yp), xytext = (10, 10), textcoords = 'offset points',size='large', color='r')

    draw()

    print("Number of O-points: "+str(n_opoint))

    if n_opoint == 0 :
        raise RuntimeError("No O-points! Giving up on this equilibrium")


  #;;;;;;;;;;;;;; Find plasma centre ;;;;;;;;;;;;;;;;;;;
  # Find the O-point closest to the middle of the grid
  
    mind = (opt_ri[0] - (old_div(numpy.float(nx),2.)))**2 + (opt_zi[0] - (old_div(numpy.float(ny),2.)))**2
    ind = 0
    for i in range (1, n_opoint) :
        d = (opt_ri[i] - (old_div(numpy.float(nx),2.)))**2 + (opt_zi[i] - (old_div(numpy.float(ny),2.)))**2
        if d < mind :
            ind = i
            mind = d
    
    primary_opt = ind
    print("Primary O-point is at "+ numpy.str(numpy.interp(opt_ri[ind],numpy.arange(numpy.size(R)),R)) + ", " + numpy.str(numpy.interp(opt_zi[ind],numpy.arange(numpy.size(Z)),Z)))
    print("")
  
    if n_xpoint > 0 :

    # Find the primary separatrix

    # First remove non-monotonic separatrices
        nkeep = 0
        for i in range (n_xpoint) :
      # Draw a line between the O-point and X-point

            n = 100 # Number of points
            farr = numpy.zeros(n)
            dr = old_div((xpt_ri[i] - opt_ri[ind]), numpy.float(n))
            dz = old_div((xpt_zi[i] - opt_zi[ind]), numpy.float(n))
            for j in range (n) :
                # interpolate f at this location
                func = RectBivariateSpline(Rx, Zx, F)

                farr[j] = func(opt_ri[ind] + dr*numpy.float(j), opt_zi[ind] + dz*numpy.float(j))
       

            # farr should be monotonic, and shouldn't cross any other separatrices

            maxind = numpy.argmax(farr)
            minind = numpy.argmin(farr)
            if (maxind < minind) : maxind, minind = minind, maxind

        # Allow a little leeway to account for errors
        # NOTE: This needs a bit of refining
            if (maxind > (n-3)) and (minind < 3) :
            # Monotonic, so add this to a list of x-points to keep
                if nkeep == 0 :
                    keep = [i]
                else:
                    keep = numpy.append(keep, i)
                
                
                nkeep = nkeep + 1
       

        if nkeep > 0 :
            print("Keeping x-points ", keep)
            xpt_ri = xpt_ri[keep]
            xpt_zi = xpt_zi[keep]
            xpt_f = xpt_f[keep]
        else:
            "No x-points kept"
        
        n_xpoint = nkeep

     
        # Now find x-point closest to primary O-point
        s = numpy.argsort(numpy.abs(opt_f[ind] - xpt_f))
        xpt_ri = xpt_ri[s]
        xpt_zi = xpt_zi[s]
        xpt_f = xpt_f[s]
        inner_sep = 0
       
    else:

    # No x-points. Pick mid-point in f
   
        xpt_f = 0.5*(numpy.max(F) + numpy.min(F))
    
        print("WARNING: No X-points. Setting separatrix to F = "+str(xpt_f))

        xpt_ri = 0
        xpt_zi = 0
        inner_sep = 0
    


    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    # Put results into a structure
  
    result = Bunch(n_opoint=n_opoint, n_xpoint=n_xpoint, # Number of O- and X-points
            primary_opt=primary_opt, # Which O-point is the plasma centre
            inner_sep=inner_sep, #Innermost X-point separatrix
            opt_ri=opt_ri, opt_zi=opt_zi, opt_f=opt_f, # O-point location (indices) and psi values
            xpt_ri=xpt_ri, xpt_zi=xpt_zi, xpt_f=xpt_f) # X-point locations and psi values
  
    return result

