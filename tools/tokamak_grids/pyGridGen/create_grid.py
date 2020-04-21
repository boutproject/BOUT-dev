from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
# Tokamak grid generator
# ======================
# 
# Generates a flux-surface aligned grid from
# an R-Z mesh of psi values. 
#
# Features:
# --------
#  o An arbitrary number of X-points
#  o Automatic default settings when not
#    supplied. 
#
# Author: Ben Dudson, University of York, Nov 2009
# 
#
# NOTE: Throughout, "F" means un-normalised psi,
#       and "psi" means normalised psi
#
# Useage:
# ------
# 
# FUNCTION create_grid, F, R, Z, settings, critical=critical,
#                                boundary=boundary, fpsi=fpsi
#
# F - psi(nr, nz) 2D array
# R - R(nr)  1D array
# Z - Z(nz)  1D array
# 
# Settings is a structure containing:
#   
#   psi_inner, psi_outer  - Range of normalised psi
#   nrad                  - Number of radial grid points
#                           Scalar -> Total number. Distributed
#                                     automatically
#                           Array -> Specified for each section
#   rad_peaking           - Radial separatrix peaking factor
#                           Not supplied -> 0
#                           scalar -> same for all regions
#                           array -> Different for each region
#   npol                  - Number of poloidal points.
#                           Scalar -> Total number. Distributed
#                                     automatically
#                           Array -> Specified for each section
#   pol_peaking           - Poloidal peaking factor
#                           Not supplied -> 0
#                           Scalar -> Same for all poloidal sections
#                           Array -> Specified for each section
# 
# A minimal settings structure must contain 4 scalars:
#    psi_inner, psi_outer, nrad and npol.
#
# Critical is a structure containing critical points:
#  
#   n_opoint, n_xpoint   - Number of O- and X-points
#   primary_opt          - Index of plasma centre O-point
#   inner_sep            - X-point index of inner separatrix
#   opt_ri, opt_zi       - R and Z indices for each O-point
#   opt_f                - Psi value at each O-point
#   xpt_ri, xpt_zi       - R and Z indices for each X-point
#   xpt_f                - Psi value of each X-point
# 
# If critical is not supplied, it is calculated
# 
# Boundary is a 2D array of indices
#   boundary[0, *]       - R values
#   boundary[1, *]       - Z values
#
# fpsi is an (optional) current function as a 2D array
#   fpsi[0,*]            - Psi values
#   fpsi[1,*]            - f values
#
# Return structure contains:
# 
#   error                       - Non-zero if an error occurred
#   psi_inner, psi_outer        - Normalised ranges of psi used
#   nrad[d], npol[d]            - Number of grid points used in each domain
#   yup_xsplit[d]               - Upper Y edge of domain. x < xsplit -> inner
#   yup_xin[d], yup_xout[d]     - Inner and outer domain connections. -1 = none
#   ydown_xsplit[d]             - Lower Y edge of domain. x < xsplit -> inner
#   ydown_xin[d], ydown_xout[d] - Inner and outer domain connections
# 
#   Rixy[x,y], Zixy[x,y]   - 2D arrays of indices into R and Z array
#   Rxy[x,y], Zxy[x,y]     - 2D arrays of (rad,pol) grid-point locations
#   psixy[x,y]             - 2D array of normalised psi at each point
#   faxis, fnorm           - Psi normalisation factors
#   settings               - Structure containing final settings used
#   critical               - Structure describing O- and X-points
#                            (from analyse_equil.pro)
#

import sys
#import inspect
from pylab import figure, show, draw, plot, contour, setp, clabel
import numpy
from boututils.bunch import Bunch
from scipy import integrate
from local_gradient import local_gradient
from follow_gradient import follow_gradient
from ask import query_yes_no
from boututils.calculus import deriv
from smooth import SMOOTH



#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# 
# Poloidal grid
#
# Divide up a poloidal arc
#

def poloidal_grid(interp_data, R, Z, ri, zi, n, fpsi=None, parweight=None,
                    ydown_dist=None, yup_dist=None, 
                    ydown_space=None, yup_space=None):

    if parweight  == None:
        parweight = 0.0  # Default is poloidal distance

    np = numpy.size(ri)

    #if np == 0 :
    #    # Calculate poloidal distance along starting line
    #    drdi = numpy.gradient(numpy.interp(ri, numpy.arange(R.size).astype(float), R))
    #    dzdi = numpy.gradient(numpy.interp(zi, numpy.arange(Z.size).astype(float), Z))
    # 
    #    dldi = numpy.sqrt(drdi**2 + dzdi**2)
    #    poldist = int_func(numpy.arange(0.,np), dldi) # Poloidal distance along line
    #else:
    rpos = numpy.interp(ri,numpy.arange(R.size).astype(float),R)
    zpos = numpy.interp(zi,numpy.arange(Z.size).astype(float),Z)
        
    dd = numpy.sqrt((zpos[1::] - zpos[0:(np-1)])**2 + (rpos[1::] - rpos[0:(np-1)])**2)
    dd = numpy.append(dd, numpy.sqrt((zpos[0] - zpos[np-1])**2 + (rpos[0] - rpos[np-1])**2))
    poldist = numpy.zeros(np)
    for i in range (1,np) :
        poldist[i] = poldist[i-1] + dd[i-1]
    
    if numpy.ndim(fpsi) == 2:
        # Parallel distance along line
        # Need poloidal and toroidal field
        ni = numpy.size(ri)
        bp = numpy.zeros(ni)
        bt = numpy.zeros(ni)
        m = interp_data.method
        interp_data.method = 2
        
        for i in range (ni) :
           
            out=local_gradient( interp_data, ri[i], zi[i], status=0, f=0., dfdr=0., dfdz=0.)
            f=out.f[0][0]
            dfdr=out.dfdr[0][0]
            dfdz=out.dfdz[0][0]
            status=out.status
            

            # dfd* are derivatives wrt the indices. Need to multiply by dr/di etc
            dfdr /= numpy.interp(ri[i],numpy.arange(R.size).astype(float),numpy.gradient(R))
            dfdz /= numpy.interp(zi[i],numpy.arange(Z.size).astype(float),numpy.gradient(Z))
      
            if i == 0 : # F=Bt*R at primary O-point
                btr = numpy.interp(f, fpsi[0,:], fpsi[1,:])
    
      
            bp[i] = old_div(numpy.sqrt(dfdr**2 + dfdz**2), numpy.interp(ri[i], numpy.arange(R.size).astype(float), R))
            bt[i] = numpy.abs( old_div(btr, numpy.interp( ri[i], numpy.arange(R.size).astype(float),R )))
     
        interp_data.method = m
        b = numpy.sqrt(bt**2 + bp**2)
        ddpar = dd * b / bp
        pardist = numpy.zeros(np)
        for i in range (1,np) :
            ip = (i + 1) % np
            pardist[i] = pardist[i-1] + ddpar[i] #0.5*(ddpar[i-1] + ddpar[ip])
     
    else:
        pardist = poldist # Just use the same poloidal distance

    print("PARWEIGHT: ", parweight)

    dist = parweight*pardist + (1. - parweight)*poldist
    
    # Divide up distance. No points at the end (could be x-point)
    if n >= 2 :
        if ydown_dist == None :
            ydown_dist = dist[np-1]* 0.5 / numpy.float(n)
        if yup_dist == None :
            yup_dist = dist[np-1] * 0.5 / numpy.float(n)

        if ydown_space == None :
            ydown_space = ydown_dist
        if yup_space == None :
            yup_space = ydown_dist
    #dloc = (dist[np-1] - ydown_dist - yup_dist) * FINDGEN(n)/FLOAT(n-1) + ydown_dist  # Distance locations
    
        fn = numpy.float(n-1)
        d = (dist[np-1] - ydown_dist - yup_dist) # Distance between first and last
        i = numpy.arange(0.,n)
        
        #yd = numpy.min( ydown_space, 0.5*d/fn ) # np.min requires int-type inputs H.SETO
        #yu = numpy.min( yup_space , 0.5*d/fn )  # np.min requires int-type inputs H.SETO
        yd = numpy.minimum(ydown_space, 0.5*d/fn )
        yu = numpy.minimum(yup_space , 0.5*d/fn )

        # Fit to ai + bi^2 + c[i-sin(2pi*i/(n-1))*(n-1)/(2pi)]
        a = yd*2.
        b = old_div((2.*yu - a), fn)
        c = old_div(d,fn) - a - 0.5*b*fn
                
        dloc =  ydown_dist + a*i + 0.5*b*i**2 + c*(i - numpy.sin(2.*numpy.pi*i / fn)*fn/(2.*numpy.pi))
        
        ddloc = a  + b*i + c*(1. - numpy.cos(2.*numpy.pi*i / fn))
       
        #; Fit to dist = a*i^3 + b*i^2 + c*i
        #;c = ydown_dist*2.
        #;b = 3.*(d/fn^2 - c/fn) - 2.*yup_dist/fn + c/fn
        #;a = d/fn^3 - c/fn^2 - b/fn
        #;dloc = ydown_dist + c*i + b*i^2 + a*i^3
    
    else :
        print("SORRY; Need 2 points in each region")
        sys.exit("Error message")


    # Get indices in ri, zi
    ind = numpy.interp(dloc, dist, numpy.arange(0.,np))

    return ind


#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# Create a grid around a given line
#
# Arguments:
#   interp_data  Data used for interpolation of psi(x,y)
#   R, Z         
#   ri, zi       1D indices into F for starting line
#   f0           Starting f for the line
#   fin, fout    Range of f to grid 
#   npar         Number of perpendicular lines to generate
#   nin, nout    Number of points inside and outside line
def grid_region ( interp_data, R, Z, 
                  ri, zi,        # Starting line to grid.
                  fvals,         # Location of the surfaces
                  sind,          # Index in fvals of the starting line
                  npar,          # Number of points along the line
                  slast=None,    # Index in fvals of last successful point
                  sfirst=None, 
                  oplot=None, 
                  boundary=None, 
                  ffirst=None, flast=None, 
                  fpsi=None,  # f(psi) = R*Bt optional current function
                  parweight=None, # Space equally in parallel (1) or poloidal (0) distance
                  ydown_dist=None, yup_dist=None, 
                  ydown_space=None, yup_space=None ):
  
    nsurf = numpy.size(fvals)
  
    if sind >= 0 :
        # starting position is on one of the output surfaces
        f0 = fvals[sind]
        nin = sind
    else:
        # Starting position between surfaces
        n = old_div(numpy.size(ri),2)
        out=local_gradient (interp_data, ri[n], zi[n], status=0, f=f0)
        status=out.status
        f0=out.f[0][0]
    
        if fvals[0] < fvals[nsurf-1] :
            w = numpy.where(fvals < f0, nin)
        else:
            w = numpy.where(fvals >= f0, nin)
     
   
    nout = nsurf - nin - 1

    sfirst = 0      # Innermost successful index
    slast = nsurf-1 # Last successful index 

    ffirst = fvals[sfirst]
    flast = fvals[slast]

    print("    => Gridding range: ", numpy.min(fvals), numpy.max(fvals))
  
    nr = interp_data.nx
    nz = interp_data.ny


    ind = poloidal_grid(interp_data, R, Z, ri, zi, npar, fpsi=fpsi, 
                      ydown_dist=ydown_dist, yup_dist=yup_dist, 
                      ydown_space=ydown_space, yup_space=yup_space, 
                      parweight=parweight)
      
    rii = numpy.interp(ind, numpy.arange(ri.size).astype(float), ri)
    zii = numpy.interp(ind, numpy.arange(zi.size).astype(float), zi)
    

    #rii = int_func(SMOOTH(deriv(rii), 3)) + rii[0]
    #zii = int_func(SMOOTH(deriv(zii), 3)) + zii[0]
    #STOP
  
    # Refine the location of the starting point
    for i in range (npar) :
        ri1=0.
        zi1=0.
        out=follow_gradient( interp_data, R, Z, rii[i], zii[i], f0, ri1, zi1 )
        ri1=out.rinext
        zi1=out.zinext
        
        rii[i] = ri1
        zii[i] = zi1
        

    # From each starting point, follow gradient in both directions
    
    rixy = numpy.zeros((nsurf, npar))
    zixy = numpy.zeros((nsurf, npar))
    status=0
    
    print('Starting')

    for i in range (npar) :
    
      #  print 'i=', i

        if sind >= 0 :
            rixy[nin, i] = rii[i]
            zixy[nin, i] = zii[i]
        else:
            # fvals[nin] should be just outside the starting position
            ftarg = fvals[nin]
            
            rinext=0.
            zinext=0.
            
            out=follow_gradient( interp_data, R, Z, rii[i], zii[i],
                            ftarg, rinext, zinext, status=0 )
            status=out.status
            rinext=out.rinext
            zinext=out.zinext
            rixy[nin, i] = rinext
            zixy[nin, i] = zinext
        
            
            
        for j in range (nout) :
            #print nin+j+1
            ftarg = fvals[nin+j+1]
      
            rinext=0.
            zinext=0.

            out=follow_gradient ( interp_data, R, Z, rixy[nin+j, i], zixy[nin+j, i],  
                                ftarg, rinext, zinext, status=status,  
                                boundary=boundary, fbndry=0.)
               
            status=out.status
            rinext=out.rinext
            zinext=out.zinext
            fbndry=out.fbndry
            
            
            #print "create_grid"
            #print ftarg, rinext, zinext

             

            if status == 1 :
                rixy[nin+j+1, i] = -1.0
                if nin+j < slast :
                    slast = nin+j # last good surface index
                fbndry = fvals[slast]
                if (fvals[1] - fvals[0])*(flast - fbndry) > 0 :
                    flast = 0.95*fbndry + 0.05*f0
                    break
                
            elif status == 2 :
        # Hit a boundary 
                rixy[nin+j+1, i] = rinext
                zixy[nin+j+1, i] = zinext
                if nin+j < slast :
                    slast = nin+j # Set the last point
                if (fvals[1] - fvals[0])*(flast - fbndry) > 0 :
                    flast = 0.95*fbndry + 0.05*f0
                break
                
            else :
                rixy[nin+j+1, i] = rinext
                zixy[nin+j+1, i] = zinext
            
            
            
            
      
        for j in range (nin) :
           # print nin-j-1
            ftarg = fvals[nin-j-1]
                 
                  
            rinext=0.
            zinext=0.

 
            out=follow_gradient ( interp_data, R, Z, rixy[nin-j, i], zixy[nin-j, i],  
                            ftarg, rinext, zinext, status=status,  
                                boundary=boundary, fbndry=fbndry )
            status=out.status
            rinext=out.rinext
            zinext=out.zinext
            fbndry=out.fbndry
            
            #print ftarg, rinext, zinext
       
       
            if status == 1 :
                rixy[nin-j-1, i] = -1.0
                if nin-j > sfirst :
                    sfirst = nin-j
                fbndry = fvals[sfirst]
                if (fvals[1] - fvals[0])*(ffirst - fbndry) < 0 :
                    ffirst = 0.95*fbndry + 0.05*f0
                break
      

            rixy[nin-j-1, i] = rinext
            zixy[nin-j-1, i] = zinext


            if status == 2 :
                if nin-j > sfirst :
                    sfirst = nin-j
                if (fvals[1] - fvals[0])*(ffirst - fbndry) < 0 :
                    ffirst = 0.95*fbndry + 0.05*f0
                break
           
            

        #print rixy[:,i], zixy[:,i]
        
        
        if oplot != None:
            numpy.interp(rixy[:, i], numpy.arange(R.size).astype(float), R)
            numpy.interp(zixy[:, i], numpy.arange(Z.size).astype(float), Z)
            figure (0)
            oplot_contour(rixy[:, i], zixy[:, i] , R, Z)

             
    # Set [:,0] to branch cut locations (for core only case) by H.Seto
    y_in = numpy.argmin(rixy[0,:])
    rixy = numpy.concatenate((rixy[:,y_in:],rixy[:,:y_in]),axis=1)
    zixy = numpy.concatenate((zixy[:,y_in:],zixy[:,:y_in]),axis=1)

    return Bunch(rixy=rixy, zixy=zixy, rxy=numpy.interp(rixy, numpy.arange(R.size), R), zxy=numpy.interp(zixy, numpy.arange(Z.size), Z))


    
   
    
    
# Return the contour lines of a given set of levels
# Just a small wrapper around the built-in CONTOUR
def contour_lines( z, x, y, levels=20):
                      
    if numpy.size(y) != numpy.size(z) :
        contour( z, levels=levels)
    else:
        contour(x, y, z, levels=levels)
        
    return
    
    
    
    
    
    
    
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# Main gridding function
#
def create_grid( F, R, Z, in_settings, critical,
                      boundary=None,  
                      iter=None, 
                      fpsi=None, fast=None):#, # f(psi) = R*Bt current function
                      #nrad_flexible,
                      #single_rad_grid, 
                      #xpt_mindist, xpt_mul, strictbndry,debug):

    # if size(nrad_flexible) == 0 :
    #    nrad_flexible = 0

    if iter==None:
        iter = 0
        
    if iter > 3:
        print("ERROR: Too many iterations")
        return #, {error:1}

    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    # Check the settings
    # If a setting is missing, set a default value
  
    # inspect.getargspec(create_grid)
    
    # if N_PARAMS() LT 3 THEN BEGIN
    #  PRINT, "ERROR: Need at least a 2D array of psi values, R and Z arrays"
    #  RETURN, {error:1}
    #ENDIF ELSE IF N_PARAMS() LT 4 THEN BEGIN
    #  ; Settings omitted. Set defaults
    #  print "Settings not given -> using default values"
    # settings = {psi_inner:0.9, 
    #            psi_outer:1.1, 
    #           nrad:36, 
    #          npol:64, 
    #         rad_peaking:0.0, 
    #        pol_peaking:0.0, 
    #       parweight:0.0}
    # ENDIF ELSE BEGIN
    # print "Checking settings"
    settings = in_settings # So the input isn't changed
    #  str_check_present, settings, 'psi_inner', 0.9
    #  str_check_present, settings, 'psi_outer', 1.1
    #  str_check_present, settings, 'nrad', 36
    #  str_check_present, settings, 'npol', 64
    #  str_check_present, settings, 'rad_peaking', 0.0
    #  str_check_present, settings, 'pol_peaking', 0.0
    #  str_check_present, settings, 'parweight', 0.0
    #ENDELSE

    s  = numpy.ndim(F)
    s1 = numpy.shape(F)
    if s != 2:
        print("ERROR: First argument must be 2D array of psi values")
        return #, {error:1}
    nx = s1[0]
    ny = s1[1]
  
    s  = numpy.ndim(R)
    s1 = numpy.size(R)
    if s != 1  or s1 != nx :
        print("ERROR: Second argument must be 1D array of major radii")
        return # {error:1}
  
    s  = numpy.ndim(Z)
    s1 = numpy.size(Z)
    if s != 1  or s1 != ny:
        print("ERROR: Second argument must be 1D array of heights")
        return # {error:1}


    # Get an even number of points for efficient FFTs
    if nx % 2 == 1:
        # odd number of points in R. Cut out last point
        #R = R[0:(nx-1)]    # H.SETO (QST)
        #F = F[0:(nx-1), :] # H.SETO (QST)
        R = R[:-1]
        F = F[:-1, :]
        nx = nx - 1

    if ny % 2 == 1:
        # odd number of points in Z. Cut out last point
        #Z = Z[0:(ny-1)]   # H.SETO (QST)
        #F = F[:,0:(ny-1)] # H.SETO (QST)
        Z = Z[:-1]
        F = F[:,:-1]
        ny = ny - 1

    #if boundary != None:# for python3 H.SETO (QST)
    if not boundary is None: 
        s = numpy.ndim(boundary)
        s1= numpy.shape(boundary)
        if s != 2  or s1[0] != 2:
            print("WARNING: boundary must be a 2D array: [2, n]. Ignoring")
            boundary = 0
        else:       
            # Calculate indices 
            #bndryi = numpy.zeros((2,1188)) 
            #bndryi[0,:] = numpy.interp(bndryi[0,:], R, numpy.arange(0.,nx))
            #bndryi[1,:] = numpy.interp(bndryi[1,:], Z, numpy.arange(0.,ny))
            bndryi = numpy.zeros(boundary.shape,dtype=float)  # H.SETO (QST)
            bndryi[0,:] = numpy.interp(boundary[0,:], R, numpy.arange(0.,nx))
            bndryi[1,:] = numpy.interp(boundary[1,:], Z, numpy.arange(0.,ny))
            #print(bndryi[0,:]*(R[-1]-R[0])/(nx-1.)+R[0]-boundary[0,:])
            #print(bndryi[1,:]*(Z[-1]-Z[0])/(ny-1.)+Z[0]-boundary[1,:])

        #if bndryi == None : # for python3 H.SETO (QST)
        if bndryi is None :
            bndryi = numpy.zeros((2,4))
            #bndryi[0,:] = [1, nx-1, nx-1, 1] # H.SETO (QST)
            #bndryi[1,:] = [1, 1, ny-1, ny-1] # H.SETO (QST)
            bndryi[0,:] = [1, nx-2, nx-2, 1]
            bndryi[1,:] = [1, 1, ny-2, ny-2]

    #;;;;;;;;;;;;;; Psi interpolation data ;;;;;;;;;;;;;;
    
    interp_data = Bunch(nx=nx, ny=ny, 
                        method=0, 
                        f= F)       # Always include function
               
    if fast == 'fast':
        print("Using Fast settings")
        interp_data.method = 2
        
    #;;;;;;;;;;;;;;; First plot ;;;;;;;;;;;;;;;;
    
    #nlev = 100
    #minf = numpy.min(F)
    #maxf = numpy.max(F)
    #levels = numpy.arange(numpy.float(nlev))*(maxf-minf)/numpy.float(nlev-1) + minf

    Rr=numpy.tile(R,ny).reshape(ny,nx).T
    Zz=numpy.tile(Z,nx).reshape(nx,ny)
        
    #contour( Rr, Zz, F, levels=levels)
    
    #  arrange the plot on the screen      
    #mngr = get_current_fig_manager()
    #geom = mngr.window.geometry()
    #x,y,dx,dy = geom.getRect()
    #mngr.window.setGeometry(100, 100, dx, dy)


  
    #if boundary != None :
    #if not boundary is None :
    #    plot(boundary[0,:],boundary[1,:],'r--')
  
    #show(block=False)  
    
    
    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    n_opoint = critical.n_opoint
    n_xpoint = critical.n_xpoint
    primary_opt = critical.primary_opt
    inner_sep   = critical.inner_sep
    opt_ri = critical.opt_ri
    opt_zi = critical.opt_zi
    opt_f  = critical.opt_f
    xpt_ri = numpy.array(critical.xpt_ri).flatten()
    xpt_zi = numpy.array(critical.xpt_zi).flatten()
    xpt_f  = numpy.array(critical.xpt_f).flatten()
    
    # Refining x-point location is omitted H.SETO (QST)
    

    # Overplot the separatrices, O-points
    # oplot_critical, F, R, Z, critical

    # Psi normalisation factors

    faxis = opt_f[primary_opt]
    
    fnorm = xpt_f[inner_sep] - opt_f[primary_opt]
    


    # From normalised psi, get range of f
    f_inner = faxis + numpy.min(settings.psi_inner)*fnorm
    f_outer = faxis + numpy.max(settings.psi_outer)*fnorm
    
    # Check the number of x-points
    #if critical.n_xpoint == 0 :
    if n_xpoint == 0 :
        #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        # Grid entirely in the core
    	print("Generating grid entirely in the core")

    print("Now pyGridGen can hadle only the core") #H.SETO (QST)
    nrad = numpy.sum(settings.nrad) # Add up all points
    npol = numpy.sum(settings.npol)
    rad_peaking = settings.rad_peaking[0] # Just the first region
    pol_peaking = settings.pol_peaking[0]

    # work out where to put the surfaces in the core
    fvals = radial_grid(nrad, f_inner, f_outer, 1, 1, [xpt_f[inner_sep]], rad_peaking)

    fvals = fvals.flatten()
    # Create a starting surface (midpoint in radial direction)
    sind = numpy.int(old_div(nrad, 2))
    start_f = fvals[sind]
    

    #  contour_lines( F, numpy.arange(nx).astype(float), numpy.arange(ny).astype(float), levels=[start_f])
    cs=contour( Rr, Zz, F,  levels=[start_f])

    p = cs.collections[0].get_paths()

    #
    #  You might get more than one contours for the same start_f. We need to keep the closed one  
    vn=numpy.zeros(numpy.size(p))
    
      
    v = p[0].vertices
    vn[0]=numpy.shape(v)[0]
    xx=[v[:,0]]
    yy=[v[:,1]]

    if numpy.shape(vn)[0] > 1:
        for i in range(1,numpy.shape(vn)[0]):
            v = p[i].vertices
            vn[i]=numpy.shape(v)[0]
            xx.append(v[:,0])
            yy.append(v[:,1])
            #xx = [xx,v[:,0]]
            #yy = [yy,v[:,1]]
    
    print("PRIMARY: ", primary_opt, opt_ri[primary_opt], opt_zi[primary_opt])

    if numpy.shape(vn)[0] > 1 :
        # Find the surface closest to the o-point
        opt_r = numpy.interp(opt_ri[primary_opt], numpy.arange(len(R)), R)
        opt_z = numpy.interp(opt_zi[primary_opt], numpy.arange(len(Z)), Z)
        
        ind = closest_line(xx, yy, opt_r, opt_z)
        
        x=xx[ind]
        y=yy[ind]
        print("Contour: ", ind,opt_r,opt_z)
    else:
        ind = 0
        x=xx[0]
        y=yy[0]
          

    # plot the start_f line     
    zc = cs.collections[0]
    setp(zc, linewidth=1)
    clabel(cs, [start_f],  # label the level
           inline=1,
           fmt='%9.6f',
           fontsize=14)
     
    draw()             
   
    show(block=False)  
   
#

    ans=query_yes_no('Press enter to create grid')  
    
    if ans != 1 : 
        #show() # remove unnecessary show() to prevent hang-up H.SETO
	sys.exit()
        
    start_ri, start_zi=transform_xy(x,y,R,Z)

    ## Make sure that the line goes clockwise
    #
    m = numpy.argmax(numpy.interp(start_zi,numpy.arange(Z.size).astype(float), Z))
    if (numpy.gradient(numpy.interp(start_ri, numpy.arange(R.size).astype(float), R)))[m] < 0.0:
        # R should be increasing at the top. Need to reverse
        start_ri = start_ri[::-1]
        start_zi = start_zi[::-1]
        print('points reversed')
        

    ## Last point should be the same as the first
    #
    # Smooth and refine the starting location
    np = numpy.size(start_ri)
    s = 3
    ar=numpy.append(numpy.append(start_ri[(np-s-1):(np-1)], start_ri), start_ri[1:s+1])
    start_ri = SMOOTH(ar, window_len=s)[s+1:(np+s+1)]
    az=numpy.append(numpy.append(start_zi[(np-s-1):(np-1)], start_zi), start_zi[1:s+1])
    start_zi = SMOOTH(az, window_len=s)[s+1:(np+s+1)]
    
    #r_smooth = numpy.interp(start_ri, numpy.arange(len(R)), R)
    #z_smooth = numpy.interp(start_zi, numpy.arange(len(Z)), Z)
    #plot(r_smooth,z_smooth,'r--')
    #draw()
    #raw_input()

    for i in range (np) :
        ri1=0. # not used H.SETO 
        zi1=0. # not used H.SETO 
        out=follow_gradient( interp_data, R, Z, start_ri[i], start_zi[i], start_f, ri1, zi1 )
        status=out.status
        ri1=out.rinext # rinext = ri + integral of dri/df from start_f to target_f
        zi1=out.zinext # zinext = zi + integral of dzi/df from start_f to target_f

        start_ri[i] = ri1
        start_zi[i] = zi1

    # now start_ri and start_zi on target_f-line
    a = grid_region(interp_data, R, Z,
                    start_ri, start_zi, 
                    fvals, 
                    sind, 
                    npol, 
                    boundary=boundary,
                    fpsi=fpsi, 
                    parweight=settings.parweight, 
                    oplot='oplot')
   
    # H.SETO 

    plot( numpy.append(a.rxy[0,:], a.rxy[0,0]), numpy.append(a.zxy[0,:], a.zxy[0,0]), 'r')
    
          
    for i in range (1, nrad) :
        plot( numpy.append(a.rxy[i,:], a.rxy[i,0]), numpy.append(a.zxy[i,:], a.zxy[i,0]), 'r')
    
    
    for i in range (0, npol-1) :
        plot( a.rxy[:,i], a.zxy[:,i], 'r')
    
    draw()
    

    # Get other useful variables
    psixy = numpy.zeros((nrad, npol))
    for i in range (0, npol) :
        psixy[:,i] = old_div((fvals - faxis),fnorm) # to get normalised psi
    
    # Calculate magnetic field components
    dpsidR = numpy.zeros((nrad, npol))
    dpsidZ = numpy.zeros((nrad, npol))
    
    interp_data.method = 2
    
    

    for i in range (nrad) :
        for j in range (npol) :
            out = local_gradient( interp_data, a.rixy[i,j], a.zixy[i,j], status=0, dfdr=0., dfdz=0.)
            status=out.status
            dfdr=out.dfdr[0][0]
            dfdz=out.dfdz[0][0]
    # dfd* are derivatives wrt the indices. Need to multiply by dr/di etc
            dpsidR[i,j] = old_div(dfdr,numpy.interp(a.rixy[i,j], numpy.arange(R.size).astype(float),numpy.gradient(R))) 
            dpsidZ[i,j] = old_div(dfdz,numpy.interp(a.zixy[i,j], numpy.arange(Z.size).astype(float),numpy.gradient(Z))) 
                 

    # Set topology to connect in the core
    yup_xsplit = [nrad]
    ydown_xsplit = [nrad]
    yup_xin = [0]
    yup_xout = [-1]
    ydown_xin = [0]
    ydown_xout = [-1]

    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    # Create result structure

    result = Bunch(error=0, # Signal success
              psi_inner=settings.psi_inner, psi_outer=settings.psi_outer, # Unchanged psi range
              #nrad=nrad, npol=npol, #Number of points in radial and poloidal direction
              nrad=[nrad], npol=[npol], # Number of points in radial and poloidal direction (must be array)
              Rixy=a.rixy, Zixy=a.zixy, # Indices into R and Z of each point
              Rxy=a.rxy, Zxy=a.zxy, # Location of each grid point
              psixy=psixy, # Normalised psi for each point
              dpsidR=dpsidR, dpsidZ=dpsidZ, # Psi derivatives (for Bpol)
              faxis=faxis, fnorm=fnorm, # Psi normalisation factors
              settings=settings, # Settings used to create grid
              critical=critical, # Critical points
              yup_xsplit=yup_xsplit, # X index where domain splits (number of points in xin)
              ydown_xsplit=ydown_xsplit, 
              yup_xin=yup_xin, yup_xout=yup_xout, # Domain index to connect to
              ydown_xin=ydown_xin, ydown_xout=ydown_xout)

    return result

#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# 
# radial grid
#
# n             - number of grid points
# pin, pout     - range of psi
# seps          - locations of separatrices
# sep_factor    - separatrix peaking
# in_dp=in_dp   - Fix the dx on the lower side
# out_dp=out_dp - Fix the dx on the upper side

def radial_grid( n, pin, pout, include_in, include_out, seps, sep_factor, 
         in_dp=None, out_dp=None):
                      
    if n == 1 :
        return [0.5*(pin+pout)]


    x = numpy.arange(0.,n)
    m = numpy.float(n-1)
    if include_in==None :
        x = x + 0.5
        m = m + 0.5
  
    
    if include_out==None:
        m = m + 0.5
 
    x = old_div(x, m) 
  

    if in_dp==None and out_dp==None :
    # Neither inner or outer gradients set. Just return equal spacing
        return pin + (pout - pin)*x
  
    
    norm = (x[1] - x[0])*(pout - pin)

    if in_dp != None and out_dp != None :
    # Fit to dist = a*i^3 + b*i^2 + c*i
        c = old_div(in_dp,norm)
        b = 3.*(1. - c) - old_div(out_dp,norm) + c
        a = 1. - c - b
    elif in_dp != None :
    # Only inner set
        c = old_div(in_dp,norm)
        a = 0.5*(c-1.)
        b = 1. - c - a
    
        #a = 0      
        #c = in_dp/norm
        #b = 1. - c 
    else:
    # Only outer set. Used in PF region
    # Fit to (1-b)*x^a + bx for fixed b
        df = old_div(out_dp, norm)
        b = 0.25 < df  # Make sure a > 0
        a = old_div((df - b), (1. - b))
        vals = pin + (pout - pin)*( (1.-b)*x^a + b*x )
        return vals
  
    
    vals = pin + (pout - pin)*(c*x + b*x^2 + a*x^3)
    #STOP
    return vals

# integrate a function, always using the maximum
# number of grid-points possible for highest accuracy
#
# Changelog
# ---------
#
# 2010-05-24 Ben Dudson <bd512@york.ac.uk>
#
#    * Modified to allow calls with only one argument
#

def int_func(**kwargs):
    if iter(list(kwargs.items())) == 1 :
        f = kwargs.xin
        x = numpy.arrange(numpy.size(f))*1.
    else:
        f = kwargs.fin
        x = kwargs.xin
    
    n = numpy.size(f)
  
    g = numpy.zeros(n)
    
    if kwargs.simple != None :
     # Just use trapezium rule
     
        g[0] = 0.0
        for i in range (1, n):
            g[i] = g[i-1] + 0.5*(x[i] - x[i-1])*(f[i] + f[i-1])
         
    else:
     
        n2 = numpy.int(old_div(n,2))
     
        g[0] = 0.0
        for i in range (n2, n) :
            g[i] = integrate.simps( f[0:i], x[0:i] )
      
      
        for i in range (1, n2) :
            g[i] = g[n-1] - integrate.simps( f[i:], x[i:] )
      
   
    
    return g 




# Find the closest contour line to a given point
def closest_line(x, y, ri, zi, mind=None):
    n = len(x)
    if len(y) != n:
        raise ValueError("Length of x and y lists must be equal")
    
    mind = numpy.min( (x[0] - ri)**2 + (y[0] - zi)**2 )
    ind = 0
    for i in range (1, n) :
        d = numpy.min( (x[i] - ri)**2 + (y[i] - zi)**2 )
        if d < mind :
            mind = d
            ind = i
    return ind 


def transform_xy(ri, zi, R, Z):
    
    nx=numpy.size(R)
    ny=numpy.size(Z)
    x=numpy.arange(nx).astype(float)
    y=numpy.arange(ny).astype(float)    

    
    xp=numpy.interp(ri, R, x)
    yp=numpy.interp(zi, Z, y)
      
    return xp,yp
    
    
def oplot_contour( ri, zi, R, Z):
                      
    x1=numpy.interp(ri, numpy.arange(numpy.size(R)).astype(float), R)
    y1=numpy.interp(zi, numpy.arange(numpy.size(Z)).astype(float), Z)
  
    plot(x1, y1, 'ro')
    draw()

