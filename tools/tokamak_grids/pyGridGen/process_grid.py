from __future__ import print_function
from __future__ import division
from builtins import input
from builtins import range
from past.utils import old_div

# Takes the original R-Z data (from G-EQDSK), and the flux mesh
# from create_grid.pro
#
# o derives additional quantities on the mesh
# o Enforces force-balance
# o Calculates integrated quantities for field-aligned codes
#
# Inputs
# ------
#
# rz_grid - a structure containing
#    npsigrid - 1D normalised psi grid
#    fpol     - Poloidal current function
#    pres     - Plasma pressure in nt/m^2
#    qpsi     - q values
#
# mesh - Structure produced by create_grid.pro
#
#
# Keywords
# --------
#
# poorquality (output) - set to 1 if grid is poor quality
# gui         (input switch) - If set, uses dialogs to question users

import numpy
import time
from bunch import Bunch
import sys
from netcdf_io import file_open, file_write, file_close
from scipy.optimize import curve_fit
from gen_surface import gen_surface
from scipy import interpolate
from pylab import plot, figure, show, title, subplots_adjust, subplot, ylim, get_current_fig_manager

from ddx import DDX
from ddy import DDY
from int_y import int_y
from boututils.calculus import deriv
from boututils.int_func import int_func

from ask import query_yes_no
from surface import SURFACE
from scipy.optimize import root
from smooth import SMOOTH
from curvature import curvature
from adjust_jpar import adjust_jpar
from spline_mono import spline_mono
from rz_curvature import rz_curvature
import copy


#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# Average over flux-surfaces
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

def surface_average( var, mesh ):
    f = copy.deepcopy(var)
    status = gen_surface(mesh=mesh) # Start generator
    while True:
        period, yi, xi, last = gen_surface(last=None, xi=None)
        f[xi,yi] = numpy.mean(var[xi,yi]) # Average over this surface
        if last == 1: break

    return f


#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# Calculate f = R * Bt
# Using LSODE method
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

## Function for use by LSODE to integrate Bt
#def bt_differential( x, Bt):
#    global bt_com, psi, ax, bx
#    a = numpy.interp(ax, psi, x)
#    b = numpy.interp(bx, psi, x)
#
#    return a*Bt + b/Bt


def solve_f ( Rxy, psixy, pxy, Bpxy, hthe):

    MU = 4.e-7*numpy.pi

    s = numpy.shape(Rxy)
    nx = s[0]
    ny = s[1]

    a = old_div(-DDX(psixy, Rxy), Rxy)
    b = -MU*DDX(psixy, pxy) - Bpxy*DDX(Bpxy*hthe)/hthe

#    CATCH, theError
#    IF theError EQ 0 THEN BEGIN
#    # Call LSODE to follow gradient
#
#    ENDIF ELSE BEGIN
#
#    ENDELSE


def force_balance ( psixy, Rxy, Bpxy, Btxy, hthe, pxy):
    MU =4.e-7*numpy.pi

    a = old_div(DDX(psixy, Rxy), Rxy)
    b = MU*DDX(psixy, pxy) - Bpxy*DDX(psixy, Bpxy*hthe)/hthe

    return DDX(psixy, Btxy) + a*Btxy + old_div(b,Btxy)


#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# Calculate toroidal field
# Using NEWTON method
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

def Bt_func ( Bt , psi, a, b):
    #global  psi, a, b

    return deriv( psi, Bt ) + a*Bt + old_div(b, Bt)


def newton_Bt ( psixy, Rxy, Btxy, Bpxy, pxy, hthe, mesh):
    #global  psi, a, b
    MU = 4.e-7*numpy.pi

    s = numpy.shape(Rxy)
    nx = s[0]
    ny = s[1]

    axy = old_div(DDX(psixy, Rxy), Rxy)
    bxy = MU*DDX(psixy, pxy) - Bpxy*DDX(psixy, Bpxy*hthe)/hthe

    Btxy2 = numpy.zeros((nx, ny))
    for i in range (ny) :
        psi = psixy[:,i]
        a = axy[:,i]
        b = bxy[:,i]
        print("Solving f for y=", i)
        sol=root(Bt_func, Btxy[:,i], args=(psi, a, b) )
        Btxy2[:,i] = sol.x



    # Average f over flux surfaces
    fxy = surface_average(Btxy2*Rxy, mesh)

    return old_div(fxy, Rxy)


#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# Solve for pressure and f=R*Bt using force balance
# Using CURVEFIT routine
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

# Calculate jxb force at every point, given f on flux-surfaces
# X = X  - dummy, not used
# profiles = [mu0*p(surface), f(surface)]
# force = force(nx, ny), reformed into single vector
#def jxb_funct( X, profiles, force, pder) :
#    global jxb_com, nx, ny, indxy, psixy, Rxy, hthe, axy
#
#    nsurf = numpy.size(profiles) / 2
#
#    dpdx = profiles[0:(nsurf-1)]
#    f = profiles[nsurf::]
#
#    # Put profiles into 2D arrays
#    fxy = numpy.zeros((nx, ny)).astype(float)
#    dpxy = numpy.zeros((nx, ny)).astype(float)
#
#    for x in range (nx) :
#        for y in range (ny) :
#            i = indxy[x,y]
#            fxy[x,y] = f[i]
#            dpxy[x,y] = dpdx[i]
#
#    # Components of B
#    Btxy = fxy / Rxy
#
#    force = ( Btxy*hthe*DDX(psixy, Btxy)
#             + fxy^2*axy
#             + hthe*dpxy )
#
#    pder = numpy.zeros((nx, ny, 2*nsurf))
#
#    # Diagonal dependencies (dF / dfi)
#    dFdfi = ( (hthe/Rxy)*DDX(psixy, Btxy)
#            + 2.*fxy*axy )
#
#    # Set the elements of pder
#    for x in range (nx) :
#        for y in range (ny) :
#            # Get indices into profiles
#            xp = x+1 < nx-1
#            xm = x-1 > 0
#            i = indxy[x,y]
#            ip = indxy[xp,y]
#            im = indxy[xm,y]
#
#            # f components
#            pder[x,y, nsurf+i]  = dFdfi[x,y]
#            dx = psixy[xp,y] - psixy[xm,y]
#            pder[x,y, nsurf+ip] = pder[x,y, nsurf+ip] + hthe[x,y]*fxy[x,y]/ (Rxy[x,y]*Rxy[xp,y]*dx)
#            pder[x,y, nsurf+im] = pder[x,y, nsurf+im] - hthe[x,y]*fxy[x,y]/ (Rxy[x,y]*Rxy[xm,y]*dx)
#
#            # p component
#            pder[x,y, i] = hthe[x,y]
#
#    force = numpy.reshape(force, nx*ny)
#    pder = numpy.reshape(pder, nx*ny, 2*nsurf)
#
#
#def fit_profiles( mesh, psixy, Rxy, hthe, Bpxy, Btxy, dpdx ):
#    global jxb_com, nx, ny, indxy, psi, R, h, axy
#
#    MU = 4.e-7*numpy.pi
#
#    psi = psixy
#    r = Rxy
#    h = hthe
#
#    s = numpy.shape(Rxy)
#    nx = s[0]
#    ny = s[1]
#
#    # Map between location in xy and surface number
#    indxy = numpy.zeros((nx, ny)).astype(int)
#
#    status = gen_surface(mesh=mesh) # Start generator
#    i = 0
#    while True:
#        period, yi, xi, last = gen_surface(last=None, xi=None, period=None)
#        indxy[xi,yi] = i
#
#        if i == 0 :
#            farr = [numpy.mean(Btxy[xi,yi]*Rxy[xi,yi])]
#            parr = [numpy.mean(dpdx[xi,yi])]
#        else:
#            farr = [farr, numpy.mean(Btxy[xi,yi]*Rxy[xi,yi])]
#            parr = [parr, numpy.mean(dpdx[xi,yi])]
#
#
#        i = i + 1
#        if last == 1: break
#
#
#    nsurf = numpy.size(farr)
#
#    profiles = [MU*parr, farr]
#
#    # Calculate useful quantities
#    axy = hthe*DDX(psixy, Rxy)/(Rxy^3)
#
#    fit = curve_fit(jxb_funct,numpy.arrange(nx*ny).astype(float),
#                 numpy.reshape(-Bpxy*DDX(psixy, Bpxy*hthe), nx*ny))
#
#                # weights,
#                 #profiles,
#                 #function_name='jxb_funct', noder)
#
#    Btxy2 = numpy.zeros((nx, ny))
#    dpdx2 = numpy.zeros((nx, ny))
#
#    status = gen_surface(mesh=mesh) # Start generator
#    i = 0
#    while True:
#        period, yi, xi, last = gen_surface(last=None, xi=None, period=None)
#        Btxy2[xi, yi] = profiles[nsurf+i] / Rxy[xi,yi]
#        dpdx2[xi, yi] = profiles[i]
#        i = i + 1
#        if last == 1 : break


#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#; Correct hthe using force balance

def new_hfunc ( h, psi, a, b, h0, fixpos ):
   # global psi, fixpos, h0, a, b


    if fixpos == 0 :
        h2 = numpy.append(h0, h)
    elif fixpos == numpy.size(psi)-1 :
        h2 = numpy.append(h, h0)
    else:
        h2 = numpy.append(numpy.append(h[0:(fixpos)], h0), h[fixpos::])

    f = a*h2 + b*deriv( psi, h2)

    if fixpos == 0 :
        f = f[1::]
    elif fixpos == numpy.size(psi)-1 :
        f = f[0:(numpy.size(f)-1)]
    else:
        f = numpy.append(f[0:(fixpos)], f[(fixpos+1)::])


    return f


def correct_hthe ( Rxy, psixy, Btxy, Bpxy, hthe, pressure, fixhthe=None):
   # global xarr, fixpos, h0, a, b

    s = numpy.shape(Rxy)
    nx = s[0]
    ny = s[1]

    MU = 4.e-7*numpy.pi

    if fixhthe==None : fixhthe = 0
    if fixhthe < 0 : fixhthe = 0
    if fixhthe > nx-1 : fixhthe = nx-1

    fixpos = fixhthe
    print("FIX = ", fixhthe)

    axy =( Btxy*DDX(psixy, Btxy) + Bpxy*DDX(psixy, Bpxy)
        + Btxy**2*DDX(psixy, Rxy)/Rxy + MU*DDX(psixy, pressure))
    bxy = Bpxy**2

    nh = numpy.zeros((nx, ny))
    nh[fixhthe,:] = hthe[fixhthe,:]
    for i in range (ny) :
        print("Correcting y index ", i)
        xarr = psixy[:,i]
        a = axy[:,i]
        b = bxy[:,i]
        h0 = hthe[fixhthe,i]

        if fixhthe == 0 :
            htmp = hthe[1::,i]
        elif fixhthe >= nx-1 :
            # fix last point
            htmp = hthe[0:(nx-1),i]
            fixhthe = nx-1
        else:
            # fix somewhere in the middle
            htmp = numpy.append(hthe[0:(fixhthe),i], hthe[(fixhthe+1)::,i])


        sol = root(new_hfunc, htmp , args=(xarr,a,b,h0,fixpos))
        htmp=sol.x

        if fixhthe == 0 :
            nh[1::] = htmp
        elif fixhthe >= nx-1 :
            nh[0:(nx-1), i] = htmp
        else:
            nh[0:(fixhthe), i] = htmp[0:(fixhthe)]
            nh[(fixhthe+1)::, i]  = htmp[fixhthe::]

        w = numpy.size(numpy.where(nh[:,i] < 0.0))
        if w > 0 :
            print("Error in hthe solver: Negative solution at y = ", i)
            #sys.exit()



    return nh


#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#; Refine an equilibrium
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#def grid_newt ( data ):
#
#    global  nx, ny, psixy, gs_f, dpdpsi, R0, Z0, min_f, xfix
#
#    n = nx*ny
#    dxin = numpy.reshape(data, nx-1, ny)
#    dx = numpy.zeros((nx, ny))
#    if xfix <= 0 :
#        dx[1::,:] = dxin
#    elif xfix >= (nx-1) :
#        dx[0:(nx-2),:] = dxin
#    else:
#        dx[0:(xfix-1),:] = dxin[0:(nfix-1),:]
#        dx[(xfix+1)::,:] = dxin[nfix::,:]
#
#
#    xpos = dx
#    for i in range (nx) : xpos[i,:] = xpos[i,:] + i
#
#    Rxy = numpy.zeros(nx, ny)
#    Zxy = Rxy
#
#    for y in range (ny) :
#        Rxy[:,y] = numpy.interp(R0[:,y], numpy.arrange(nx).astype(float), xpos[:,y], spline)
#        Zxy[:,y] = numpy.interp(Z0[:,y], numpy.arrange(nx).astype(float), xpos[:,y], spline)
#
#
#    # calculate Bpxy, Btxy and hthe
#
#    Btxy = numpy.zeros((nx, ny))
#    for x in range (nx) : Btxy[x,:] = gs_f[x] / Rxy[x,:]
#    hthe = calc_hthe(Rxy, Zxy)
#    Bpxy = calc_bp(psixy, Rxy, Zxy)
#
#    F = -1.0*calc_force(psixy, Bpxy, Btxy, hthe, Rxy, dpdpsi)
#
#    fm = numpy.max(numpy.abs(F))
#
#    if (fm < min_f) or (min_f < 0.0) :
#        min_f = fm
#        print numpy.max(numpy.abs(Rxy - R0)), numpy.max(numpy.abs(Zxy - Z0)), numpy.max(numpy.abs(F))
#
#
#    #!P.multi=[0,0,2,0,0]
#    #surface, Bpxy, chars=3
#    #surface, F, chars=3
#
#    if yind < 0 :
#        val = numpy.reshape(F, n)
#    else: val = F[:,yind]
#
#    return val
#

#FUNCTION refine_equlibrium, mesh
#
#END

#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#; Main grid processing routine
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

def process_grid( rz_grid, mesh, output=None, poorquality=None,
                  gui=None, parent=None, reverse_bt=None,
                  curv=None, smoothpressure=None,
                  smoothhthe=None, smoothcurv=None,
                  settings=None):

    if settings==None :
        # Create an empty structure
        settings = Bunch(dummy=0)

    # Check settings
    #    settings.calcp= -1
    #    settings.calcbt= -1
    #    settings.calchthe= -1
    #    settings.calcjpar= -1

    settings.calcp   = -1# H.Seto (QST) indent is fixed 
    settings.calcbt  = -1# H.Seto (QST) indent is fixed 
    settings.calchthe= -1# H.Seto (QST) indent is fixed 
    settings.calcjpar= -1# H.Seto (QST) indent is fixed 

    # ;CATCH, err
    # ;IF err NE 0 THEN BEGIN
    # ;  PRINT, "PROCESS_GRID failed"
    # ;  PRINT, "   Error message: "+!ERROR_STATE.MSG
    # ;  CATCH, /cancel
    # ;  RETURN
    # ;ENDIF
    
    MU = 4.e-7*numpy.pi
    
    poorquality = 0
    
    if output==None : output="bout.grd.nc"
    
    # Size of the mesh
    nx = numpy.int(numpy.sum(mesh.nrad))
    ny = numpy.int(numpy.sum(mesh.npol))

    # Find the midplane
    ymid = 0
    status = gen_surface(mesh=mesh) # Start generator

    while True:
        period, yi, xi, last = gen_surface(period=None, last=None, xi=None)

        if period :
            rm = numpy.max(mesh.Rxy[xi,yi])
            ymidindx = numpy.argmax(mesh.Rxy[xi,yi])
            ymid = yi[ymidindx]
            break

        if last==1: break


    Rxy = numpy.asarray(mesh.Rxy)
    Zxy = numpy.asarray(mesh.Zxy)

    psixy = mesh.psixy*mesh.fnorm + mesh.faxis # Non-normalised psi

    pressure = numpy.zeros((nx, ny))


    
    # Use splines to interpolate pressure profile
    status = gen_surface(mesh=mesh) # Start generator
    while True:
        # Get the next domain
        period, yi, xi, last = gen_surface(period=period, last=last, xi=xi)

        if period :
            # Pressure only given on core surfaces
            # pressure[xi,yi] = SPLINE(rz_grid.npsigrid, rz_grid.pres, mesh.psixy[xi,yi[0]], /double)
            sol=interpolate.UnivariateSpline(rz_grid.npsigrid, rz_grid.pres,s=1)
            pressure[xi,yi] =sol(mesh.psixy[xi,yi[0]])
            
        else:
            # pressure is set to that at LCFS in SOL region H.SETO (QST)
            pressure[xi,yi] = rz_grid.pres[numpy.size(rz_grid.pres)-1]

        if last==1 : break


    # Add a minimum amount
    if numpy.min(pressure) < 1.0e-2*numpy.max(pressure) :
        print("****Minimum pressure is very small:", numpy.min(pressure))
        print("****Setting minimum pressure to 1% of maximum")
        pressure = pressure + 1e-2*numpy.max(pressure)


    if smoothpressure != None :
        p0 = pressure[:,ymid] # Keep initial pressure for comparison
        while True :
            #!P.multi=[0,0,2,0,0]
            fig=figure()
            plot( p0, xtitle="X index", ytitle="pressure at y="+numpy.strip(numpy.str(ymid),2)+" dashed=original", color=1, lines=1)
            plot( pressure[:,ymid], color=1)
            plot( deriv(p0), xtitle="X index", ytitle="DERIV(pressure)", color=1, lines=1)
            plot( deriv(pressure[:,ymid]), color=1 )
            sm = query_yes_no("Smooth pressure profile?")#, gui=gui, dialog_parent=parent)
            if sm :
                # Smooth the pressure profile

                p2 = pressure
                for i in range (6) :
                    status = gen_surface(mesh=mesh) # Start generator
                    while True :
                        # Get the next domain
                        period, yi, xi, last = gen_surface(period=period, last=last, xi=xi)

                        if (xi > 0) and (xi < (nx-1)) :
                            for j in range (numpy.size(yi)) :
                                p2[xi,yi[j]] = ( 0.5*pressure[xi,yi[j]] +
                                                0.25*(pressure[xi-1,yi[j]] + pressure[xi+1,yi[j]])
                                                )



                        # Make sure it's still constant on flux surfaces
                        p2[xi,yi] = numpy.mean(p2[xi,yi])
                        if last != None : break
                    pressure = p2


            if sm == 0 : break


    if numpy.min(pressure) < 0.0 :
        print("")
        print("============= WARNING ==============")
        print("Poor quality equilibrium: Pressure is negative")
        print("")
        poorquality = 1


    dpdpsi = DDX(psixy, pressure)


    #;IF MAX(dpdpsi)*mesh.fnorm GT 0.0 THEN BEGIN
    #;  PRINT, ""
    #;  PRINT, "============= WARNING =============="
    #;  PRINT, "Poor quality equilibrium: Pressure is increasing radially"
    #;  PRINT, ""
    #;  poorquality = 1
    #;ENDIF

    # Grid spacing
    dx = numpy.zeros((nx, ny))
    for y in range (ny) :
        dx[0:(nx-1),y] = psixy[1::,y] - psixy[0:(nx-1),y]
        dx[nx-1,y] = dx[nx-2,y]


    # Sign
    bpsign = 1.
    xcoord = psixy
    if numpy.min(dx) < 0. :
        bpsign = -1.
        dx = -dx # dx always positive
        xcoord = -xcoord


    dtheta = 2.*numpy.pi / numpy.float(ny)
    dy = numpy.zeros((nx, ny)) + dtheta


    # B field components
    # Following signs mean that psi increasing outwards from
    # core to edge results in Bp clockwise in the poloidal plane
    # i.e. in the positive Grad Theta direction.

    Brxy = old_div(mesh.dpsidZ, Rxy)
    Bzxy = old_div(-mesh.dpsidR, Rxy)
    Bpxy = numpy.sqrt(Brxy**2 + Bzxy**2)


    # Determine direction (dot B with grad y vector)

    dot = ( Brxy[0,ymid]*(Rxy[0,ymid+1] - Rxy[0,ymid-1]) +
            Bzxy[0,ymid]*(Zxy[0,ymid+1] - Zxy[0,ymid-1])
            )

    if dot < 0. :
        print("**** Poloidal field is in opposite direction to Grad Theta -> Bp negative")
        Bpxy = -Bpxy
        if bpsign > 0 : sys.exit() # Should be negative
        bpsign = -1.0
    else:
        if bpsign < 0 : sys.exit() # Should be positive
        bpsign = 1.


    # Get toroidal field from poloidal current function fpol
    Btxy = numpy.zeros((nx, ny))
    fprime = numpy.zeros((nx, ny))
    fp = deriv(rz_grid.npsigrid*(rz_grid.sibdry - rz_grid.simagx), rz_grid.fpol) # x= (psi-psia)/(psib-psia)


    status = gen_surface(mesh=mesh) # Start generator
    while True:
        # Get the next domain
        period, yi, xi, last = gen_surface(period=period, last=period, xi=xi)

        if period :
            # In the core
            #fpol = numpy.interp(rz_grid.fpol, rz_grid.npsigrid, mesh.psixy[xi,yi], /spline)

            sol=interpolate.UnivariateSpline(rz_grid.npsigrid, rz_grid.fpol,s=1)
         #  fpol = SPLINE(rz_grid.npsigrid, rz_grid.fpol, mesh.psixy[xi,yi[0]], 'double')
            fpol = sol(mesh.psixy[xi,yi[0]])

            sol=interpolate.UnivariateSpline(rz_grid.npsigrid, fp ,s=1)
           # fprime[xi,yi] = SPLINE(rz_grid.npsigrid, fp, mesh.psixy[xi,yi[0]], 'double')
            fprime[xi,yi] = sol(mesh.psixy[xi,yi[0]])

        else:
            # Outside core. Could be PF or SOL
            fpol = rz_grid.fpol[numpy.size(rz_grid.fpol)-1]
            fprime[xi,yi] = 0.

        Btxy[xi,yi] = old_div(fpol, Rxy[xi,yi])

        if last ==1 : break

    # Total B field
    Bxy = numpy.sqrt(Btxy**2 + Bpxy**2)


    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    # Go through the domains to get a starting estimate
    # of hthe
    hthe = numpy.zeros((nx, ny))

    #   Pick a midplane index
    status = gen_surface(mesh=mesh) # Start generator
    while True:
        # Get the next domain
        period, yi, xi, last = gen_surface(period=period, last=last, xi=xi)

        if period :
            # In the core
            rmax = numpy.argmax(Rxy[xi,yi])
            ymidplane = yi[rmax]
            break

        if last == 1: break

    status = gen_surface(mesh=mesh) # Start generator
    while True:
    # Get the next domain
        period, yi, xi, last = gen_surface(period=period, last=last, xi=xi)

        n = numpy.size(yi)

        # Get distance along this line

        if period :

            # Periodic, so can use FFT
            #drdi = REAL_PART(fft_deriv(Rxy[xi, yi]))
            #dzdi = REAL_PART(fft_deriv(Zxy[xi, yi]))
            line=numpy.append(Rxy[xi,yi[n-1::]], Rxy[xi,yi])
            line=numpy.append(line,Rxy[xi,yi[0:1]])
            # evaluate drdi by FDM rather than FFT H.SETO (QST)
            drdi = deriv(line)[1:n+1]

            line=numpy.append(Zxy[xi,yi[n-1::]], Zxy[xi,yi])
            line=numpy.append(line,Zxy[xi,yi[0:1]])

            dzdi = deriv(line)[1:n+1]
        else:
            # Non-periodic
            drdi = numpy.gradient(Rxy[xi, yi])
            dzdi = numpy.gradient(Zxy[xi, yi])


        dldi = numpy.sqrt(drdi**2 + dzdi**2)


        if 0 : # always false 

        # Need to smooth to get sensible results
            if period :
                n = numpy.size(dldi)
                line=numpy.append(dldi[(n-2)::], dldi) # once
                line=numpy.append(line,dldi[0:2])
                dldi = SMOOTH(line, 5)[4:(n+4)]

                line=numpy.append(dldi[(n-2)::], dldi) #twice
                line=numpy.append(line,dldi[0:2])
                dldi = SMOOTH(line, 5)[4:(n+4)]

                line=numpy.append(dldi[(n-2)::], dldi) # three
                line=numpy.append(line,dldi[0:2])
                dldi = SMOOTH(line, 5)[4:(n+4)]

            else:
                line = dldi
                dldi = SMOOTH(line, 5)[2:n+2]
                line = dldi
                dldi = SMOOTH(line, 5)[2:n+2]
                line = dldi
                dldi = SMOOTH(dldi, 5)[2:n+2]

        
        hthe[xi, yi] = old_div(dldi, dtheta) # First estimate of hthe

        # Get outboard midplane
        if period and xi == 0 :
            m = numpy.argmax(Rxy[0,yi])
            ymidplane = yi[m]

        if last == 1 : break

    print("Midplane index ", ymidplane)

    fb0 = force_balance(psixy, Rxy, Bpxy, Btxy, hthe, pressure)
    print("Force imbalance: ", numpy.mean(numpy.abs(fb0)), numpy.max(numpy.abs(fb0)))




    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    # Correct pressure using hthe
    
    print("Calculating pressure profile from force balance")
    
    try:
        
        # Calculate force balance
        dpdx = old_div(( -Bpxy*DDX(xcoord, Bpxy * hthe) - Btxy*hthe*DDX(xcoord, Btxy) - (Btxy*Btxy*hthe/Rxy)*DDX(xcoord, Rxy) ), (MU*hthe))

        # Surface average
        dpdx2 = surface_average(dpdx, mesh)

        pres = numpy.zeros((nx, ny))
        # Integrate to get pressure
        for i in range (ny) :
            pres[:,i] = int_func(psixy[:,i], dpdx2[:,i])
            pres[:,i] = pres[:,i] - pres[nx-1,i]



        status = gen_surface(mesh=mesh) # Start generator
        while True:
            # Get the next domain
            period, yi, xi, last = gen_surface(period=None, last=None, xi=None)

            ma = numpy.max(pres[xi,yi])

            for i in range (numpy.size(yi)) :
                pres[:,yi[i]] = pres[:,yi[i]] - pres[xi,yi[i]] + ma

            if last == 1 : break


        pres = pres - numpy.min(pres)

        # Some sort of smoothing here?


        fb0 = force_balance(psixy, Rxy, Bpxy, Btxy, hthe, pres)
        print("Force imbalance: ", numpy.mean(numpy.abs(fb0)), numpy.max(numpy.abs(fb0)))


        #!P.MULTI=[0,0,2,0,0]
        fig=figure(figsize=(7, 11))
        subplots_adjust(left=.07, bottom=.07, right=0.95, top=0.95,
                wspace=.3, hspace=.25)

        SURFACE( pressure, fig, xtitle="X", ytitle="Y", var='Pa', sub=[2,1,1])
        title("Input pressure")
        SURFACE( pres, fig, xtitle="X", ytitle="Y", var='Pa', sub=[2,1,2])
        title("New pressure")
        #  arrange the plot on the screen
        #      mngr = get_current_fig_manager()
        #      geom = mngr.window.geometry()
        #      x,y,dx,dy = geom.getRect()
        #      mngr.window.setGeometry(0, 0, dx, dy)
        #
        show(block=False)


        calcp = settings.calcp

        if calcp == -1 :
            calcp = query_yes_no("Keep new pressure?")#, gui=gui, dialog_parent=parent)
        else: time.sleep( 2 )
        if calcp == 1 :
            pressure = pres
            dpdpsi = dpdx2


    except Exception:
        print("WARNING: Pressure profile calculation failed: ")#, !ERROR_STATE.MSG
        pass

    #CATCH, /cancel

    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    # Correct f = RBt using force balance

    calcbt = settings.calcbt
    if calcbt == -1 : calcbt = query_yes_no("Correct f=RBt using force balance?")#, gui=gui, dialog_parent=parent)
    if calcbt == 1 :

        new_Btxy = newton_Bt(psixy, Rxy, Btxy, Bpxy, pres, hthe, mesh)

        fb0 = force_balance(psixy, Rxy, Bpxy, new_Btxy, hthe, pressure)
        print("force imbalance: ", numpy.mean(numpy.abs(fb0)), numpy.max(numpy.abs(fb0)))


        fig=figure(figsize=(7, 11))
        subplots_adjust(left=.07, bottom=.07, right=0.95, top=0.95,
                wspace=.3, hspace=.25)

        subplot(211)
        SURFACE( Btxy, fig, xtitle="X", ytitle="Y", var='T', sub=[2,1,1])
        title("Input Bt")
        subplot(212)
        SURFACE( new_Btxy, fig, xtitle="X", ytitle="Y", var='T', sub=[2,1,2])
        title("New Bt")
        #  arrange the plot on the screen
        #mngr = get_current_fig_manager()
        #geom = mngr.window.geometry()
        #x,y,dx,dy = geom.getRect()
        #mngr.window.setGeometry(600, 0, dx, dy)


        show(block=False)

        calcbt = settings.calcbt
        if calcbt == -1 : calcbt = query_yes_no("Keep new Bt?")#, gui=gui, dialog_parent=parent)
        if calcbt == 1 :
            Btxy = new_Btxy
            Bxy = numpy.sqrt(Btxy**2 + Bpxy**2)

  #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  # CALCULATE HTHE
  # Modify hthe to fit force balance using initial guess
  # Does not depend on signs
  #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    calchthe = settings.calchthe
    if calchthe == -1 : calchthe = query_yes_no("Adjust hthe using force balance?")#, gui=gui, dialog_parent=parent)
    if calchthe == 1 :
        # This doesn't behave well close to the x-points
        fixhthe = numpy.int(old_div(nx, 2))
        nh = correct_hthe(Rxy, psixy, Btxy, Bpxy, hthe, pressure, fixhthe=fixhthe)

        fb0 = force_balance(psixy, Rxy, Bpxy, Btxy, nh, pressure)
        print("Force imbalance: ", numpy.mean(numpy.abs(fb0)), numpy.max(numpy.abs(fb0)))

        print("numpy.maximum difference in hthe: ", numpy.max(numpy.abs(hthe - nh)))
        print("numpy.maximum percentage difference: ", 100.*numpy.max(numpy.abs(old_div((hthe - nh),hthe))))

       #!P.multi=[0,0,1,0,0]
        fig=figure(figsize=(7, 4))
        title("Poloidal arc length at midplane. line is initial estimate")
        plot( hthe[:,0], '-' )
        plot( nh[:,0], 'r-+' )
                  #  arrange the plot on the screen
        #mngr = get_current_fig_manager()
        #geom = mngr.window.geometry()
        #x,y,dx,dy = geom.getRect()
        #mngr.window.setGeometry(0, 1150, dx, dy)

        show(block=False)

        if query_yes_no("Keep new hthe?") == 1:#, gui=gui, dialog_parent=parent) :
            hthe = nh



    if smoothhthe != None :
        # Smooth hthe to prevent large jumps in X or Y. This
        # should be done by creating a better mesh in the first place

        # Need to smooth in Y and X otherwise smoothing in X
        # produces discontinuities in Y
        hold = hthe

        if 1 : # always done but not implemented so far H.SETO (QST)
            # Nonlinear smoothing. Tries to smooth only regions with large
            # changes in gradient 

            hthe =0.# smooth_nl(hthe, mesh)

        else: # never done wrong indents have been fixed H.SETO (QST)
            # Just use smooth in both directions

            for i in range (ny) :
                hthe[:,i] = SMOOTH(SMOOTH(hthe[:,i],10),10)


            status = gen_surface(mesh=mesh) # Start generator
            while True:
                # Get the next domain
                period, yi, xi, last = gen_surface(period=None, last=None, xi=None)

                n = numpy.size(yi)

                if period :
                    hthe[xi,yi] = (SMOOTH([hthe[xi,yi[(n-4):(n-1)]], hthe[xi,yi], hthe[xi,yi[0:3]]], 4))[4:(n+3)]
                else:
                    hthe[xi,yi] = SMOOTH(hthe[xi,yi], 4)

                if last == 1: break
    


    # Calculate field-line pitch
    pitch = hthe * Btxy / (Bpxy * Rxy)


    # derivative with psi
    dqdpsi = DDX(psixy, pitch)


    qinty, qloop = int_y(pitch, mesh, loop=0, nosmooth='nosmooth', simple='simple')
    qinty = qinty * dtheta
    qloop = qloop * dtheta


    sinty = int_y(dqdpsi, mesh, nosmooth='nosmooth', simple='simple') * dtheta



    # NOTE: This is only valid in the core
    pol_angle = numpy.zeros((nx,ny))
    for i in range (nx) :  pol_angle[i, :] = 2.0*numpy.pi * qinty[i,:] / qloop[i]


  #;;;;;;;;;;;;;;;;;;; THETA_ZERO ;;;;;;;;;;;;;;;;;;;;;;
  # re-set zshift to be zero at the outboard midplane

    print("MIDPLANE INDEX = ", ymidplane)

    status = gen_surface(mesh=mesh) # Start generator
    while True:
    # Get the next domain
        period, yi, xi, last = gen_surface(period=None, last=None, xi=None)

        w = numpy.size(numpy.where(yi == ymidplane))
        if w > 0 :
      # Crosses the midplane
            qinty[xi, yi] = qinty[xi, yi] - qinty[xi, ymidplane]
            sinty[xi, yi] = sinty[xi, yi] - sinty[xi, ymidplane]
        else:
      # Doesn't include a point at the midplane
            qinty[xi, yi] = qinty[xi, yi] - qinty[xi,yi[0]]
            sinty[xi, yi] = sinty[xi, yi] - sinty[xi,yi[0]]

        if last ==1 : break

    print("")
    print("==== Calculating curvature ====")

  #;;;;;;;;;;;;;;;;;;; CURVATURE ;;;;;;;;;;;;;;;;;;;;;;;
  # Calculating b x kappa

    if curv == None :

        print("*** Calculating curvature in toroidal coordinates")

        thetaxy = numpy.zeros((nx, ny))
        status = gen_surface(mesh=mesh) # Start generator
        while True:
            # Get the next domain
            period, yi, xi, last = gen_surface(period=None, last=None, xi=None)
            thetaxy[xi,yi] = numpy.arange(numpy.size(yi)).astype(float)*dtheta
            if last ==1 : break


        bxcv = curvature( nx, ny, Rxy,Zxy, Brxy, Bzxy, Btxy,
                    psixy, thetaxy, hthe,
                     mesh=mesh)

        bxcvx = bpsign*bxcv.psi
        bxcvy= bxcv.theta
        bxcvz = bpsign*(bxcv.phi - sinty*bxcv.psi - pitch*bxcv.theta)


        # x borders
        bxcvx[0,:] = bxcvx[1,:]
        bxcvx[nx-1,:] = bxcvx[nx-2,:]

        bxcvy[0,:] = bxcvy[1,:]
        bxcvy[nx-1,:] = bxcvy[nx-2,:]

        bxcvz[0,:] = bxcvz[1,:]
        bxcvz[nx-1,:] = bxcvz[nx-2,:]

    elif curv == 1 :
        # Calculate on R-Z mesh and then interpolate onto grid
        # ( cylindrical coordinates)

        print("*** Calculating curvature in cylindrical coordinates")
        print("liner interpolation is use insted of bilinear interpolation: exit")
        exit()
        bxcv = rz_curvature(rz_grid)

        # DCT methods cause spurious oscillations
        # Linear interpolation seems to be more robust # H.SETO (QST) question
        bxcv_psi = numpy.interp(bxcv.psi, mesh.Rixy, mesh.Zixy)
        bxcv_theta = old_div(numpy.interp(bxcv.theta, mesh.Rixy, mesh.Zixy), hthe)
        bxcv_phi = numpy.interp(bxcv.phi, mesh.Rixy, mesh.Zixy)

        # If Bp is reversed, then Grad x = - Grad psi
        bxcvx = bpsign*bxcv_psi
        bxcvy = bxcv_theta
        bxcvz = bpsign*(bxcv_phi - sinty*bxcv_psi - pitch*bxcv_theta)
    elif curv == 2 :
        # Curvature from Curl(b/B)

        bxcvx = bpsign*(Bpxy * Btxy*Rxy * DDY(old_div(1., Bxy), mesh) / hthe)
        bxcvy = -bpsign*Bxy*Bpxy * DDX(xcoord, Btxy*Rxy/Bxy^2) / (2.*hthe)
        bxcvz = Bpxy^3 * DDX(xcoord, old_div(hthe,Bpxy)) / (2.*hthe*Bxy) - Btxy*Rxy*DDX(xcoord, old_div(Btxy,Rxy)) / (2.*Bxy) - sinty*bxcvx

    else:
        # calculate in flux coordinates.

        print("*** Calculating curvature in flux coordinates")

        dpb = numpy.zeros((nx, ny))      # quantity used for y and z components

        for i in range (ny) :
            dpb[:,i] = MU*dpdpsi/Bxy[:,i]

        dpb = dpb + DDX(xcoord, Bxy)

        bxcvx = bpsign*(Bpxy * Btxy*Rxy * DDY(old_div(1., Bxy), mesh) / hthe)
        bxcvy = bpsign*(Bpxy*Btxy*Rxy*dpb / (hthe*Bxy^2))
        bxcvz = -dpb - sinty*bxcvx



    if smoothcurv:
        # Smooth curvature to prevent large jumps

        # Nonlinear smoothing. Tries to smooth only regions with large
        # changes in gradient

        bz = bxcvz + sinty * bxcvx

        print("Smoothing bxcvx...")
        bxcvx = 0.#smooth_nl(bxcvx, mesh)
        print("Smoothing bxcvy...")
        bxcvy = 0.#smooth_nl(bxcvy, mesh)
        print("Smoothing bxcvz...")
        bz = 0.#smooth_nl(bz, mesh)

        bxcvz = bz - sinty * bxcvx


  #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  # CALCULATE PARALLEL CURRENT
  #
  # Three ways to calculate Jpar0:
  # 1. From fprime and pprime
  # 2. From Curl(B) in field-aligned coords
  # 3. From the curvature
  #
  # Provides a way to check if Btor should be reversed
  #
  #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    print("")
    print("==== Calculating parallel current ====")

    jpar0 = - Bxy * fprime / MU - Rxy*Btxy * dpdpsi / Bxy


    # Set to zero in PF and SOL
    status = gen_surface(mesh=mesh) # Start generator
    while True:
    # Get the next domain
        period, yi, xi, last = gen_surface(period=None, last=None, xi=None)

        if period == None : jpar0[xi,yi] = 0.0
        if last == 1 : break

  # Curl(B) expression for Jpar0 (very noisy usually)
    j0 = ( bpsign*((Bpxy*Btxy*Rxy/(Bxy*hthe))*( DDX(xcoord, Bxy**2*hthe/Bpxy) - bpsign*Btxy*Rxy*DDX(xcoord,Btxy*hthe/(Rxy*Bpxy)) )
        - Bxy*DDX(xcoord, Btxy*Rxy)) / MU )



  # Create a temporary mesh structure to send to adjust_jpar
    tmp_mesh = Bunch(mesh,
                           bxcvx=bxcvx, bxcvy=bxcvy,  bxcvz=bxcvz,
                            Bpxy=Bpxy,  Btxy=Btxy,  Bxy=Bxy,
                            dx=dx,  dy=dy,
                            hthe=hthe,  jpar0=jpar0,  pressure=pressure)
    tmp_mesh.psixy = psixy

    jpar = adjust_jpar( tmp_mesh, noplot='noplot')


   #!P.multi=[0,2,2,0,0]

    fig=figure(figsize=(15, 11))
    subplots_adjust(left=.07, bottom=.07, right=0.95, top=0.95,
                wspace=.3, hspace=.25)

    subplot(221)
    SURFACE( jpar0, fig, xtitle="X", ytitle="Y", var='A', sub=[2,2,1])
    title("Jpar from F' and P'")

    subplot(222)
    SURFACE( jpar, fig, xtitle="X", ytitle="Y", var='A', sub=[2,2,2])
    title("Jpar from curvature")

    subplot(223)
    plot( jpar0[0,:],'-', jpar[0,:] ,'+' )
    ylim([numpy.min([jpar0[0,:],jpar[0,:]]), numpy.max([jpar0[0,:],jpar[0,:]])])
    title("jpar at x=0. Solid from f' and p'")

    subplot(224)
    plot(jpar0[:,ymidplane],'-' , jpar[:,ymidplane] , '+' )
    ylim([numpy.min([jpar0[:,ymidplane],jpar[:,ymidplane]]),numpy.max([jpar0[:,ymidplane],jpar[:,ymidplane]])])

    title("Jpar at y="+numpy.str(ymidplane)+" Solid from f' and p'")

        #  arrange the plot on the screen
    #mngr = get_current_fig_manager()
    #geom = mngr.window.geometry()
    #x,y,dx,dy = geom.getRect()
    #mngr.window.setGeometry(1350, 0, dx, dy)


    show(block=False)

 # !P.multi=0

    calcjpar = settings.calcjpar
    if calcjpar == -1 : calcjpar = query_yes_no("Use Jpar from curvature?")#, gui=gui, dialog_parent=parent)
    if calcjpar == True :
        jpar0 = jpar


    if 0 :

    # Try smoothing jpar0 in psi, preserving zero points and maxima
        jps = jpar0
        for y in range ( ny ):
            j = jpar0[:,y]
            js = j
            ma = numpy.max(numpy.abs(j))
            ip = numpy.argmax(numpy.abs(j))
            if (ma < 1.e-4) or (ip == 0) :
                jps[:,y] = j

            level = 1.
            #i0 = MAX(WHERE(ABS(j[0:ip]) LT level))
            i1 = numpy.min(numpy.where(numpy.abs(j[ip::]) < level))

            #IF i0 LE 0 THEN i0 = 1
            i0 = 1

            if i1 == -1 :
                i1 = nx-2
            else:
                i1 = i1 + ip

            if (ip <= i0) or (ip >= i1) :

      # Now preserve starting and end points, and peak value
                div = numpy.int(old_div((i1-i0),10))+1 # reduce number of points by this factor

                inds = [i0] # first point
                for i in [i0+div, ip-div, div] : inds = [inds, i]
                inds = [inds, ip] # Put in the peak point

        # Calculate spline interpolation of inner part

                js[0:ip] = spline_mono(inds, j[inds], numpy.arange(ip+1),
                             yp0=(j[i0] - j[i0-1]), ypn_1=0.0)

                inds = [ip] # peak point
                for i in [ip+div, i1-div, div] :
                    inds = [inds, i]


                inds = [inds, i1]  # Last point
                js[ip:i1] = spline_mono(inds, j[inds], ip+numpy.arange(i1-ip+1),
                              yp0=0.0, ypn_1=(j[i1+1]-j[i1]))

                jps[:,y] = js



  #;;;;;;;;;;;;;;;;;;; TOPOLOGY ;;;;;;;;;;;;;;;;;;;;;;;
  # Calculate indices for backwards-compatibility

    nr = numpy.size(mesh.nrad)
    np = numpy.size(mesh.npol)
    if (nr == 2) and (np == 3) :
        print("Single null equilibrium")

        ixseps1 = mesh.nrad[0]
        ixseps2 = nx

        jyseps1_1 = mesh.npol[0]-1
        jyseps1_2 = mesh.npol[0] + numpy.int(old_div(mesh.npol[1],2))
        ny_inner = jyseps1_2
        jyseps2_1 = jyseps1_2
        jyseps2_2 = ny - mesh.npol[2]-1

    elif (nr == 3) and (np == 6) :
        print("Double null equilibrium")

        ixseps1 = mesh.nrad[0]
        ixseps2 = ixseps1 + mesh.nrad[1]

        jyseps1_1 = mesh.npol[0]-1
        jyseps2_1 = jyseps1_1 + mesh.npol[1]

        ny_inner = jyseps2_1 + mesh.npol[2] + 1

        jyseps1_2 = ny_inner + mesh.npol[3] - 1
        jyseps2_2 = jyseps1_2 + mesh.npol[4]

    elif (nr == 1) and (np == 1) :

        print("Single domain")

        ixseps1 = nx
        ixseps2 = nx

        jyseps1_1 = -1
        jyseps1_2 = numpy.int(old_div(ny,2))
        jyseps2_1 = numpy.int(old_div(ny,2))
        ny_inner = numpy.int(old_div(ny,2))
        jyseps2_2 = ny - 1

    else:
        print("***************************************")
        print("* WARNING: Equilibrium not recognised *")
        print("*                                     *")
        print("*  Check mesh carefully!              *")
        print("*                                     *")
        print("*  Contact Ben Dudson                 *")
        print("*      benjamin.dudson@york.ac.uk     *")
        print("***************************************")
        ixseps1 = -1
        ixseps2 = -1

        jyseps1_1 = -1
        jyseps1_2 = numpy.int(old_div(ny,2))
        jyseps2_1 = numpy.int(old_div(ny,2))
        ny_inner = numpy.int(old_div(ny,2))
        jyseps2_2 = ny - 1


    print("Generating plasma profiles:")

    print("  1. Flat temperature profile")
    print("  2. Flat density profile")
    print("  3. Te proportional to density")
    while True:
        opt = input("Profile option:")
        if eval(opt) >= 1 and eval(opt) <= 3 : break


    if eval(opt) == 1 :
        # flat temperature profile

        print("Setting flat temperature profile")
        while True:
            Te_x = eval(input("Temperature (eV):"))


        # get density
            Ni = old_div(pressure, (2.* Te_x* 1.602e-19*1.0e20))

            print("numpy.maximum density (10^20 m^-3):", numpy.max(Ni))

            done = query_yes_no("Is this ok?")
            if done == 1 : break

        Te = numpy.zeros((nx, ny))+Te_x
        Ti = Te
        Ni_x = numpy.max(Ni)
        Ti_x = Te_x
    elif eval(opt) == 2 :
        print("Setting flat density profile")

        while True:
            Ni_x = eval(input("Density [10^20 m^-3]:"))

            # get temperature
            Te = old_div(pressure, (2.* Ni_x * 1.602e-19*1.0e20))

            print("numpy.maximum temperature (eV):", numpy.max(Te))
            if query_yes_no("Is this ok?") == 1 : break

        Ti = Te
        Ni = numpy.zeros((nx, ny)) + Ni_x
        Te_x = numpy.max(Te)
        Ti_x = Te_x
    else:
        print("Setting te proportional to density")

        while True:
            Te_x = eval(input("Maximum temperature [eV]:"))


            Ni_x = old_div(numpy.max(pressure), (2.*Te_x * 1.602e-19*1.0e20))

            print("Maximum density [10^20 m^-3]:", Ni_x)

            shape = numpy.sqrt(pressure / numpy.max(pressure))
            Te = Te_x * shape
            Ni = Ni_x * shape
            if query_yes_no("Is this ok?") == 1 : break
        Ti = Te
        Ti_x =  Te_x


    rmag = numpy.max(numpy.abs(Rxy))
    print("Setting rmag = ", rmag)

    bmag = numpy.max(numpy.abs(Bxy))
    print("Setting bmag = ", bmag)

    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    # save to file
    # open a new netCDF file for writing.
    handle = file_open(output)

    print("Writing grid to file "+output)

    # Size of the grid

    s = file_write(handle, "nx", nx)
    s = file_write(handle, "ny", ny)


    # Topology for original scheme
    s = file_write(handle, "ixseps1", ixseps1)
    s = file_write(handle, "ixseps2", ixseps2)
    s = file_write(handle, "jyseps1_1", jyseps1_1)
    s = file_write(handle, "jyseps1_2", jyseps1_2)
    s = file_write(handle, "jyseps2_1", jyseps2_1)
    s = file_write(handle, "jyseps2_2", jyseps2_2)
    s = file_write(handle, "ny_inner", ny_inner);

    # Grid spacing

    s = file_write(handle, "dx", dx)
    s = file_write(handle, "dy", dy)

    s = file_write(handle, "ShiftAngle", qloop)
    s = file_write(handle, "zShift", qinty)
    s = file_write(handle, "pol_angle", pol_angle)
    s = file_write(handle, "ShiftTorsion", dqdpsi)

    s = file_write(handle, "Rxy",  Rxy)
    s = file_write(handle, "Zxy",  Zxy)
    s = file_write(handle, "Bpxy", Bpxy)
    s = file_write(handle, "Btxy", Btxy)
    s = file_write(handle, "Bxy",  Bxy)
    s = file_write(handle, "hthe", hthe)
    s = file_write(handle, "sinty", sinty)
    s = file_write(handle, "psixy", psixy)

    # Topology for general configurations
    s = file_write(handle, "yup_xsplit", mesh.yup_xsplit)
    s = file_write(handle, "ydown_xsplit", mesh.ydown_xsplit)
    s = file_write(handle, "yup_xin", mesh.yup_xin)
    s = file_write(handle, "yup_xout", mesh.yup_xout)
    s = file_write(handle, "ydown_xin", mesh.ydown_xin)
    s = file_write(handle, "ydown_xout", mesh.ydown_xout)
    s = file_write(handle, "nrad", mesh.nrad)
    s = file_write(handle, "npol", mesh.npol)

    # plasma profiles

    s = file_write(handle, "pressure", pressure)
    s = file_write(handle, "Jpar0", jpar0)
    s = file_write(handle, "Ni0", Ni)
    s = file_write(handle, "Te0", Te)
    s = file_write(handle, "Ti0", Ti)


    s = file_write(handle, "Ni_x", Ni_x)
    s = file_write(handle, "Te_x", Te_x)
    s = file_write(handle, "Ti_x", Ti_x)
    s = file_write(handle, "bmag", bmag)
    s = file_write(handle, "rmag", rmag)

    # Curvature
    s = file_write(handle, "bxcvx", bxcvx)
    s = file_write(handle, "bxcvy", bxcvy)
    s = file_write(handle, "bxcvz", bxcvz)

    # Psi range
    s = file_write(handle, "psi_axis", mesh.faxis)
    psi_bndry = mesh.faxis + mesh.fnorm
    s = file_write(handle, "psi_bndry", psi_bndry)

    file_close, handle
    print("DONE")

    #!P.multi=[0,0,1,0,0]
