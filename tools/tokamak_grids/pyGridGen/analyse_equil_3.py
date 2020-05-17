# >>>>> for python3 >>>>>
from __future__ import print_function 
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
# <<<<< for python3 <<<<<
import numpy as np 
from boututils.bunch import Bunch
from scipy.interpolate import RectBivariateSpline
from pylab import contour, scatter, annotate, draw
from crosslines import find_inter

def analyse_equil ( F, R, Z):
    """
    Equilibrium analysis routine
    Takes a RZ psi grid, and finds x-points and o-points
    
    F - F(nr, nz) 2D array of psi values
    R - R(nr) 1D array of major radii
    Z - Z(nz) 1D array of heights
    
    Returns a structure of critical points containing:

    n_opoint, n_xpoint   - Number of O- and X-points
    primary_opt          - Index of plasma centre O-point
    inner_sep            - X-point index of inner separatrix
    opt_ri, opt_zi       - R and Z indices for each O-point
    opt_f                - Psi value at each O-point
    xpt_ri, xpt_zi       - R and Z indices for each X-point
    xpt_f                - Psi value of each X-point

    change log:
    2017-04-27 H.SETO (QST)
    * numpy -> np, nx -> nr, ny -> nz, Rr -> Rrz, Zz -> Zrz 
    * fix indent of comment 
    * add some comments 
    """

    # get array size
    [nr,nz] = F.shape
    
    #;;;;;;;;;;;;;;; Find critical points ;;;;;;;;;;;;;
    #
    # Need to find starting locations for O-points (minima/maxima)
    # and X-points (saddle points)
    
    Rrz = np.tile(R,nz).reshape(nz,nr).T 
    Zrz = np.tile(Z,nr).reshape(nr,nz)
        
    #; Use contour to get crossing-points where dfdr = dfdz = 0
    contour1=contour(Rrz,Zrz,np.gradient(F)[0], levels=[0.0], colors='r') 
    contour2=contour(Rrz,Zrz,np.gradient(F)[1], levels=[0.0], colors='r')   

    draw()
    
    #; find where these two cross 

    lines=find_inter(contour1, contour2)
    rex=np.interp(lines[0],  R, np.arange(R.size)).astype(int) 
    zex=np.interp(lines[1],  Z, np.arange(Z.size)).astype(int)
    
    #; check for points too close to the edge
    w=np.where((rex > 2) & (rex < nr-3) & (zex >2) & (zex < nr-3))
    nextrema = np.size(w)
    rex=rex[w].flatten()
    zex=zex[w].flatten()
    rp=lines[0][w]
    zp=lines[1][w]
    
    #;;;;;;;;;;;;;; Characterise extrema ;;;;;;;;;;;;;;;;;
    # Fit a surface through local points using 6x6 matrix
    # This is to determine the type of extrema, and to
    # refine the location
    # 
    
    n_opoint = 0
    n_xpoint = 0

    # Calculate inverse matrix
    rio = np.array([-1, 0, 0, 0, 1, 1]) # R index offsets
    zio = np.array([ 0,-1, 0, 1, 0, 1]) # Z index offsets
  
    # Fitting flux funcion by quadradic formulation in RZ-plane 
    #
    #  F = C0 + C1*ri + C2*zi + C3*ri*zi + C4*ri*ri + C5*zi*zi
    #                       
    #             x---3---5  Px=(ri,zi)
    #             |   |   |  P0=(-1, 0), P1=( 0,-1)
    #             0---2---4  P2=( 0, 0), P3=( 0, 1) 
    #             |   |   |  P4=( 1, 0), P5=( 1, 1) 
    #             x---1---x
    #
    A = np.transpose([[np.ones(6,dtype=np.float)],
                      [rio], 
                      [zio], 
                      [rio*zio], 
                      [rio**2], 
                      [zio**2]])[:,0,:]
    
    for e in range (nextrema) :
      
        # Fit in index space so result is index number
        print("Critical point "+str(e))
    
        localf = np.zeros(6)
        for i in range (6) :
            # Get the f value in a stencil around this point
            
            ri = rex[e]+rio[i] # Zero-gradient at edges
            zi = zex[e]+zio[i] 
            localf[i] = F[ri, zi]
        
        # solve Ax=b to obtain coefficients Cx
        # res:    x=[C0,C1,C2,C3,C4,C5].T
        # localf: b=[F0,F1,F2,F3,F4,F5].T
        
        res, _, _, _ = np.linalg.lstsq(A,localf)
  
        # This determines whether saddle or extremum
        #  
        # det = d2Fdri2*dF2dzi2 - d2Fdridzi*d2Fdzidri
        #

        det = 4.*res[4]*res[5] - res[3]**2

        rinew = (res[3]*res[2]-2.*res[1]*res[5])/det
        zinew = (res[3]*res[1]-2.*res[4]*res[2])/det

        if det < 0.0 :
            print("   X-point:",end="")
        else:
            print("   O-point:",end="")
                
        rnew = rp[e]
        znew = zp[e]
        
        func = RectBivariateSpline(R, Z, F)
        
        fnew =func(rnew, znew)
        
        x=np.arange(np.size(R))
        y=np.arange(np.size(Z))
        rinew = np.interp(rnew,R,x)
        zinew = np.interp(znew, Z, y)
        
        print("   R = %15.6e, Z= %15.6e, F = %15.6e" % (rnew,znew,fnew[0][0]))
        
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
                        xpt_ri = np.append(xpt_ri, rinew)
                        xpt_zi = np.append(xpt_zi, zinew)
                        xpt_f  = np.append(xpt_f,  fnew)
                        n_xpoint = n_xpoint + 1
                
                                                                                
                scatter(rnew,znew,s=100, marker='x', color='r')
                
                annotate(np.str(n_xpoint-1),  xy = (rnew, znew), xytext = (10, 10),
                         textcoords = 'offset points',size='large', color='r')
                
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
                        opt_ri = np.append(opt_ri, rinew)
                        opt_zi = np.append(opt_zi, zinew)
                        opt_f  = np.append(opt_f, fnew)
                        n_opoint = n_opoint + 1
 
                scatter(rnew,znew,s=100, marker='o',color='r')
               
                annotate(np.str(n_opoint-1),  xy = (rnew, znew), xytext = (10, 10), 
                         textcoords = 'offset points', size='large', color='b')
                draw()
     
                      
    print("Number of O-points: "+np.str(n_opoint))
    print("Number of X-points: "+np.str(n_xpoint))

    
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
  
    mind = (opt_ri[0] - (old_div(np.float(nr),2.)))**2 + (opt_zi[0] - (old_div(np.float(nz),2.)))**2
    ind = 0
    for i in range (1, n_opoint) :
        d = (opt_ri[i] - (old_div(np.float(nr),2.)))**2 + (opt_zi[i] - (old_div(np.float(nz),2.)))**2
        if d < mind :
            ind = i
            mind = d
            
    primary_opt = ind
    print("Primary O-point is at "+str(np.interp(opt_ri[ind],x,R)) + ", " + str(np.interp(opt_zi[ind],y,Z)))
    print("")
  
    if n_xpoint > 0 :

    # Find the primary separatrix

    # First remove non-monotonic separatrices
        nkeep = 0
        for i in range (n_xpoint) :
      # Draw a line between the O-point and X-point

            n = 100 # Number of points
            farr = np.zeros(n)
            dr = old_div((xpt_ri[i] - opt_ri[ind]), np.float(n))
            dz = old_div((xpt_zi[i] - opt_zi[ind]), np.float(n))
            for j in range (n) :
                # interpolate f at this location
                func = RectBivariateSpline(x, y, F)

                farr[j] = func(opt_ri[ind] + dr*np.float(j), opt_zi[ind] + dz*np.float(j))
       

            # farr should be monotonic, and shouldn't cross any other separatrices

            maxind = np.argmax(farr)
            minind = np.argmin(farr)
            if (maxind < minind) : maxind, minind = minind, maxind

        # Allow a little leeway to account for errors
        # NOTE: This needs a bit of refining
            if (maxind > (n-3)) and (minind < 3) :
            # Monotonic, so add this to a list of x-points to keep
                if nkeep == 0 :
                    keep = [i]
                else:
                    keep = np.append(keep, i)
                
                
                nkeep = nkeep + 1
       

        if nkeep > 0 :
            print("Keeping x-points ", keep)
            xpt_ri = xpt_ri[keep]
            xpt_zi = xpt_zi[keep]
            xpt_f  = xpt_f[keep]
        else:
            "No x-points kept"
        
        n_xpoint = nkeep

        # Now find x-point closest to primary O-point
        s = np.argsort(np.abs(opt_f[ind] - xpt_f))
        xpt_ri = xpt_ri[s]
        xpt_zi = xpt_zi[s]
        xpt_f  = xpt_f[s]
        inner_sep = 0
    else:

    # No x-points. Pick mid-point in f
   
        xpt_f = 0.5*(np.max(F) + np.min(F))
    
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

