from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# Equilibrium analysis routine
# 
# Takes a RZ psi grid, and finds x-points and o-points
# 
# F - F(nr, nz) 2D array of psi values
# R - R(nr) 1D array of major radii
# Z - Z(nz) 1D array of heights
#
# Returns a structure of critical points containing:
#
#   n_opoint, n_xpoint   - Number of O- and X-points
#   primary_opt          - Index of plasma centre O-point
#   inner_sep            - X-point index of inner separatrix
#   opt_ri, opt_zi       - R and Z indices for each O-point
#   opt_f                - Psi value at each O-point
#   xpt_ri, xpt_zi       - R and Z indices for each X-point
#   xpt_f                - Psi value of each X-point
# 

import numpy 
from boututils.bunch import Bunch
from scipy.interpolate import RectBivariateSpline
from pylab import contour, gradient, scatter, annotate, draw
from crosslines import find_inter

def analyse_equil ( F, R, Z):
    s = numpy.shape(F)
    nx = s[0]
    ny = s[1]
  
  #;;;;;;;;;;;;;;; Find critical points ;;;;;;;;;;;;;
  #
  # Need to find starting locations for O-points (minima/maxima)
  # and X-points (saddle points)
  #
    Rr=numpy.tile(R,nx).reshape(nx,ny).T  # needed for contour
    Zz=numpy.tile(Z,ny).reshape(nx,ny)
    
 
    contour1=contour(Rr,Zz,gradient(F)[0], levels=[0.0], colors='r')
    contour2=contour(Rr,Zz,gradient(F)[1], levels=[0.0], colors='r')    

    draw()
    

### --- line crossings ---------------------------
    res=find_inter( contour1, contour2)
            
    rex=numpy.interp(res[0],  R, numpy.arange(R.size)).astype(int)
    zex=numpy.interp(res[1],  Z, numpy.arange(Z.size)).astype(int)
   
      
    w=numpy.where((rex > 2) & (rex < nx-3) & (zex >2) & (zex < nx-3))
    nextrema = numpy.size(w)
    rex=rex[w].flatten()
    zex=zex[w].flatten()
  
  #;;;;;;;;;;;;;; Characterise extrema ;;;;;;;;;;;;;;;;;
  # Fit a surface through local points using 6x6 matrix
  # This is to determine the type of extrema, and to
  # refine the location
  # 
  
    
    n_opoint = 0
    n_xpoint = 0

  # Calculate inverse matrix
    rio = numpy.array([-1, 0, 0, 0, 1, 1]) # R index offsets
    zio = numpy.array([ 0,-1, 0, 1, 0, 1]) # Z index offsets
  
  # Fitting a + br + cz + drz + er^2 + fz^2
    A = numpy.transpose([[numpy.zeros(6,numpy.int)+1],
                 [rio], 
                 [zio], 
                 [rio*zio], 
                 [rio**2], 
                 [zio**2]])
                 
    A=A[:,0,:]    
  
    for e in range (nextrema) :
      
      # Fit in index space so result is index number
        print("Critical point "+str(e))
    
     
        valid = 1

        localf = numpy.zeros(6)
        for i in range (6) :
      # Get the f value in a stencil around this point
            xi = (rex[e]+rio[i]) #> 0) < (nx-1) # Zero-gradient at edges
            yi = (zex[e]+zio[i]) #> 0) < (ny-1)
            localf[i] = F[xi, yi]
   
        res, _, _, _ = numpy.linalg.lstsq(A,localf)
  
          
    # Res now contains [a,b,c,d,e,f]
    #                  [0,1,2,3,4,5]

    # This determines whether saddle or extremum
        det = 4.*res[4]*res[5] - res[3]**2
    
        if det < 0.0 :
            print("   X-point")
        else:
            print("   O-point")
    
    
    # Get location (2x2 matrix of coefficients)
    
        rinew = old_div((res[3]*res[2] - 2.*res[1]*res[5]), det)
        zinew = old_div((res[3]*res[1] - 2.*res[4]*res[2]), det)

        if (numpy.abs(rinew) > 1.) or (numpy.abs(zinew) > 1.0) :
      #; Method has gone slightly wrong. Try a different method.
      #; Get a contour line starting at this point. Should
      #; produce a circle around the real o-point. 
            print("   Fitted location deviates too much")
            if det < 0.0 :
                print("   => X-point probably not valid")
                print("      deviation = "+numpy.str(rinew)+","+numpy.str(zinew))
                #valid = 0
#            else:
#    
#                contour_lines, F, findgen(nx), findgen(ny), levels=[F[rex[e], zex[e]]], \
#                    path_info=info, path_xy=xy
#    
#                if np.size(info) > 1 :
#                # More than one contour. Select the one closest
#                    ind = closest_line(info, xy, rex[e], zex[e])
#                    info = info[ind]
#                else: info = info[0]
#    
#                rinew = 0.5*(MAX(xy[0, info.offset:(info.offset + info.n - 1)]) + \
#                     MIN(xy[0, info.offset:(info.offset + info.n - 1)])) - rex[e]
#                zinew = 0.5*(MAX(xy[1, info.offset:(info.offset + info.n - 1)]) + \
#                     MIN(xy[1, info.offset:(info.offset + info.n - 1)])) - zex[e]
#    
#                if (np.abs(rinew) > 2.) or (np.abs(zinew) > 2.0) :
#                    print "   Backup method also failed. Keeping initial guess"
#                    rinew = 0.
#                    zinew = 0.
#                 
#

        if valid :
            fnew = (res[0] + res[1]*rinew + res[2]*zinew + 
                res[3]*rinew*zinew + res[4]*rinew**2 + res[5]*zinew**2)
      
            rinew = rinew + rex[e]
            zinew = zinew + zex[e]
      
            print("   Starting index: " + str(rex[e])+", "+str(zex[e]))
            print("   Refined  index: " + str(rinew)+", "+str(zinew))
      
            x=numpy.arange(numpy.size(R))
            y=numpy.arange(numpy.size(Z))
            rnew = numpy.interp(rinew,x,R)
            znew = numpy.interp(zinew, y, Z)
      
            print("   Position: " + str(rnew)+", "+str(znew))
            print("   F = "+str(fnew))
      
            if det < 0.0 :

                if n_xpoint == 0 :
                    xpt_ri = [rinew]
                    xpt_zi = [zinew]
                    xpt_f = [fnew]
                    n_xpoint = n_xpoint + 1
                else:
            # Check if this duplicates an existing point
                    
                    if rinew in xpt_ri and zinew in xpt_zi :
                        print("   Duplicates existing X-point.")
                    else:
                        xpt_ri = numpy.append(xpt_ri, rinew)
                        xpt_zi = numpy.append(xpt_zi, zinew)
                        xpt_f = numpy.append(xpt_f, fnew)
                        n_xpoint = n_xpoint + 1
                
                                                                                
                scatter(rnew,znew,s=100, marker='x', color='r')
                
                annotate(numpy.str(n_xpoint-1),  xy = (rnew, znew), xytext = (10, 10), textcoords = 'offset points',size='large', color='r')

                draw()
            else:

                if n_opoint == 0 :
                    opt_ri = [rinew]
                    opt_zi = [zinew]
                    opt_f = [fnew]
                    n_opoint = n_opoint + 1
                else:
            # Check if this duplicates an existing point
        
                    if rinew in opt_ri and zinew in opt_zi :
                        print("   Duplicates existing O-point")
                    else:
                        opt_ri = numpy.append(opt_ri, rinew)
                        opt_zi = numpy.append(opt_zi, zinew)
                        opt_f = numpy.append(opt_f, fnew)
                        n_opoint = n_opoint + 1
 
                scatter(rnew,znew,s=100, marker='o',color='r')
               
                annotate(numpy.str(n_opoint-1),  xy = (rnew, znew), xytext = (10, 10), textcoords = 'offset points', size='large', color='b')
                draw()
     
                      
    print("Number of O-points: "+numpy.str(n_opoint))
    print("Number of X-points: "+numpy.str(n_xpoint))

    
    if n_opoint == 0 :
            opt_ri = [rinew]
            opt_zi = [zinew]
            opt_f = [fnew]
            n_opoint = n_opoint + 1

    print("Number of O-points: "+str(n_opoint))

    if n_opoint == 0 :
        print("No O-points! Giving up on this equilibrium")
        return Bunch(n_opoint=0, n_xpoint=0, primary_opt=-1)
  

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
    print("Primary O-point is at "+str(numpy.interp(opt_ri[ind],x,R)) + ", " + str(numpy.interp(opt_zi[ind],y,Z)))
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
                func = RectBivariateSpline(x, y, F)

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

