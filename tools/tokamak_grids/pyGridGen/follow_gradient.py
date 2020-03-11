from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
from builtins import object

import sys
import numpy
from boututils.bunch import Bunch


try:
    from ode.lsode import lsode
except ImportError:
    print("No ode.lsode available. Trying scipy.integrate.odeint")
    from scipy.integrate import odeint

    # Define a class to emulate lsode behavior
    class lsode(object):
        def __init__(self, func, f0, rz0):
            # Function for odeint needs to have reversed inputs
            self._func = lambda pos, fcur : func(fcur, pos)
            self._f0 = f0
            self._rz0 = rz0

        def integrate(self, ftarget):
            solode = odeint(self._func, self._rz0, [self._f0,ftarget], full_output=True)
            self._f0 = ftarget
            self.steps = solode[1]['nst'][0]
            return solode[0][1,:]

from local_gradient import local_gradient
from itertools import chain
from boututils.calculus import deriv
from saveobject import saveobject

global rd_com, idata, lastgoodf, lastgoodpos, Rpos, Zpos, ood, tol, Ri, Zi, dR, dZ

#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# Follow the gradient from a given point to a target f
#

# Calculate dR/df and dZ/df for use by LSODE
# Input: pos[0] = R, pos[1] = Z
# Output [0] = dR/df = -Bz/B^2 , [1] = dZ/df = Br/B^2
def radial_differential( fcur, pos):


    global rd_com, idata, lastgoodf, lastgoodpos, Rpos, Zpos, ood, tol, Ri, Zi, dR, dZ

    out=local_gradient( idata, pos[0], pos[1], status=0, dfdr=0., dfdz=0. )
    status=out.status
    dfdr=out.dfdr[0][0]
    dfdz=out.dfdz[0][0]

    ood = 0

    if status ==0 :
        # No error in localgradient.
        lastgoodf = fcur
        lastgoodpos = pos

    else:
        # If status NE 0 then an error occurred.
        # Allow dfdz to cause an error so escape LSODE
        ood = 1 # Out Of Domain

  #if numpy.size(boundary) > 1 :
  #  # Got a boundary - check for crossings
  #
  #  cpos , ncross, inds2 = line_crossings([ri0, pos[0]], [zi0, pos[1]], 0,
  #                        boundary[0,:], boundary[1,:], 1, ncross=0, inds2=0)
  #  if (ncross % 2) == 1 : # Odd number of boundary crossings
  #    # Check how far away the crossing is
  #      if numpy.sqrt( (pos[0] - cpos[0,0])**2 + (pos[1] - cpos[1,0])**2 ) > tol :
  #      # Need to trigger an error
  #      #PRINT, "HIT BOUNDARY:", SQRT( (pos[0] - cpos[0,0])^2 + (pos[1] - cpos[1,0])^2 )
  #          status = 'a.bcd' # gibberish
  #

    #dRdi = numpy.interp(pos[0],numpy.arange(Rpos.size).astype(float),deriv(Rpos))
    #dZdi = numpy.interp(pos[1],numpy.arange(Zpos.size).astype(float),deriv(Zpos))
    dRdi = numpy.interp(pos[0],Ri,dR)
    dZdi = numpy.interp(pos[1],Zi,dZ)



  # Check mismatch between fcur and f ?
    Br = old_div(dfdz,dZdi)
    Bz = old_div(-dfdr,dRdi)
    B2 = Br**2 + Bz**2

    return [-Bz/B2/dRdi, Br/B2/dZdi]


#
# interp_data (in)  Data for interpolation of 2D Psi
# R, Z      (in)    1D arrays of position vs. index
# ri0,zi0   (in)    Starting indices
# ftarget   (in)    The f to aim for
# ri,zi     (out)   Final position
# status    (out)   Non-zero if hits a boundary. 1 if out of range at
#                   the start, 2 if before reaching ftarget
# boundary  (in)    Optional 2D array [2,n] of points on the boundary
# fbndry    (out)   If hits boundary, gives final f value
# ibndry    (out)   If hits boundary, index where hit
#
def follow_gradient( interp_data, R, Z, ri0, zi0, ftarget, ri, zi, status=0,
                     boundary=None, fbndry=None, ibndry=None ):

    global rd_com, idata, lastgoodf, lastgoodpos, Rpos, Zpos, ood, tol, Ri, Zi, dR, dZ

    tol = 0.1

    Rpos = R
    Zpos = Z

    Ri=numpy.arange(Rpos.size).astype(float)
    Zi=numpy.arange(Zpos.size).astype(float)
    dR=deriv(Rpos)
    dZ=deriv(Zpos)

    ibndry = -1

    idata = interp_data

    #if boundary != None :
    if not boundary is None :
        bndry = boundary
        ri0c = ri0
        zi0c = zi0
    else:
        bndry = 0
    
    ood = 0

    if ftarget==None : print(ftarget)

    # Get starting f
    out=local_gradient( interp_data, ri0, zi0, status=status, f=0., dfdr=None, dfdz=None)
    status=out.status
    f0=out.f
    if status == 1 :
        ri = ri0
        zi = zi0
        status = 1
        return Bunch(status=status, ri=ri, zi=zi)


    fmax = ftarget # Target (with maybe boundary in the way)

    # Call LSODE to follow gradient
    rzold = [ri0, zi0]
    rcount = 0

    solver = lsode(radial_differential, f0, rzold)
    rznew=solver.integrate(ftarget)
    nsteps = solver.steps

    lstat=0
    #print 'nsteps=',nsteps
    #print rzold, rznew
    #
    #sys.exit()
#
#   #         if nsteps > 100 : lstat = -1
#            if lstat == -1 :
#                print "  -> Excessive work "+str(f0)+" to "+str(ftarget)+" Trying to continue..."
#                lstat = 2 # continue
#                rcount = rcount + 1
#                if rcount > 3 :
#                    print "   -> Too many repeats. Giving Up."
#
#                    ri = lastgoodpos[0]
#                    zi = lastgoodpos[1]
#                    fmax = lastgoodf
#
#                    return Bunch(status=status,ri=ri,zi=zi)
#
#            # Get f at this new location
#                out=local_gradient( interp_data, rznew[0], rznew[1], status=status, f=f0, dfdr=None, dfdz=None)
#                status=out.status
#                f0=out.f
#
#                if status == 1 :
#                    ri = ri0
#                    zi = zi0
#                    status = 1
#                    return Bunch(status=status, rinext=ri, zinext=zi)
#
#                rzold = rznew
#
#
#            else :
#                print "I break"
#                break
#
#        print 'am I here?'
#
#        if status==0:
#            print 'I break from try?'
#            break
#
#        if lstat < 0 :
#            print "Error in LSODE routine when following psi gradient."
#            print "LSODE status: ", lstat
#           #STOP
#
#
#    except Exception as theError:
#        print theError


    ri = rznew[0]
    zi = rznew[1]

 #   else:
 #   # An error occurred in LSODE.
 #   # lastgoodf contains the last known good f value
 #   #PRINT, "CAUGHT ERROR "+!ERROR_STATE.MSG
 #   #CATCH, /cancel
 #       ri = lastgoodpos[0]
 #       zi = lastgoodpos[1]
 #       fmax = lastgoodf
 #       if ood :
 #       # Gone Out Of Domain
 #           status = 2
 #           fbndry = lastgoodf
 #       #PRINT, "Out of domain at f = ", fbndry
 #       # Repeat to verify that this does work
 #           rzold = [ri0, zi0]
 #           try :
 #               fbndry = lastgoodf - 0.1*(ftarget - f0)
 #               if theError != 0 :
 #                   print "   Error again at ", fbndry
 #
 #
 #               solver=lsode(radial_differential, f0, rzold)
 #               rznew=solver.integrate(fbndry - f0)
 #           except Exception as theError:
 #               print theError
 #
 #           return Bunch(status=status, rinext=ri, zinext=zi)
 #
 #   # Otherwise just crossed a boundary
 #
 #   #CATCH, /cancel



    #if boundary != None:
    ## Check if the line crossed a boundary
    ##PRINT, "Checking boundary ", boundary[*,1:2], [ri0, ri], [zi0, zi]
    #    cpos, ncross, inds2 = line_crossings([ri0, ri], [zi0, zi], 0,
    #                      boundary[0,:], boundary[1,:], 1, ncross=0, inds2=0)
    #    if (ncross % 2) == 1 : # Odd number of boundary crossings
    #        if numpy.sqrt( (ri - cpos[0,0])**2 + (zi - cpos[1,0])**2 ) > 0.1 :
    #    #PRINT, "FINDING BOUNDARY", SQRT( (ri - cpos[0,0])^2 + (zi - cpos[1,0])^2 )
    #    # Use divide-and-conquer to find crossing point
    #
    #            tol = 1e-4 # Make the boundary crossing stricter
    #
    #            ibndry = inds2[0] # Index in boundary where hit
    #
    #            fcur = f0 # Current known good position
    #            rzold = [ri0,zi0]
    #            rzcur = rzold
    #            while True:
    #                fbndry = (fcur + fmax) / 2
    #     # Try to go half-way to fmax
    #                #CATCH, theError
    #                if theError != 0 :
    #        # Crossed boundary. Change fmax
    #                    #CATCH, /cancel
    #                    fmax = fbndry
    #                    ibndry = inds2[0] # refined boundary index
    #                else:
    #                    solver=lsode(radial_differential, f0, rzold)
    #                    rznew=solver.integrate(fbndry - f0)
    #
    #       # Didn't cross. Make this the new current location
    #                    fcur = fbndry
    #                    rzcur = rznew
    #
    #                    #CATCH, /cancel
    #
    #                if numpy.abs(fmax - fcur) < 0.01*numpy.abs(ftarget - f0):
    #                    break
    #            ri = rzcur[0]
    #            zi = rzcur[1]
    #            fbndry = fcur
    #
    #    #PRINT, "Hit boundary", ri, zi, " f =", fbndry
    #            status = 2
    #            return Bunch(status=status, rinext=ri, zinext=zi, fbndry=fbndry, ibndry=ibndry)

    #print "follow_gradient"
    #print ri, zi
    return Bunch(status = 0, rinext=ri, zinext=zi, fbndry=fbndry, ibndry=ibndry)

#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# Line crossing detection
#
# (r1, z1) and (r2, z2) are two lines
# period1 and period2 determine whether the lines are periodic
#
# Returns a list of locations where they intersect

def line_crossings( r1, z1, period1, r2, z2, period2, ncross=0,
                         inds1=0, inds2=0 ):
    n1 = numpy.size(r1)
    n2 = numpy.size(r2)

    result = 0
    ncross = 0

    for i in range (n1) :
        ip = i + 1
        if i == n1-1 :
            if period1 :
                ip = 0
            else:
                break


        for j in range (n2-1) :
            jp = j+1
            if j == n2-1 :
                if period2 :
                    jp = 0
                else :
                    break


      # Test if line (i to ip) and (j to jp) intersects
      # cast as a 2x2 matrix

            a = r1[ip] - r1[i]
            b = r2[j] - r2[jp]
            c = z1[ip] - z1[i]
            d = z2[j] - z2[jp]

            dr = r2[j] - r1[i]
            dz = z2[j] - z1[i]

            det = a*d - b*c

      # Get location along the line segments
            if numpy.abs(det) > 1.e-6 :
                alpha = old_div((d*dr - b*dz),det)
                beta =  old_div((a*dz - c*dr),det)
            else:
                alpha = -1.
                beta = -1.


            if (alpha >= 0.0) and (alpha <= 1.0) and (beta >= 0.0) and (beta <= 1.0) :
        # Intersection

                r = r1[i] + alpha * a
                z = z1[i] + alpha * c

                if ncross == 0 :
                    result = numpy.zeros((2,1))
                    result[0,0] = r
                    result[1,0] = z

                    inds1 = [numpy.float(i)+alpha]
                    inds2 = [numpy.float(j)+beta]
                else :
                    rold = result
                    result = numpy.zeros((2, ncross+1))
                    result[:,0:(ncross-1)] = rold
                    result[0,ncross] = r
                    result[1,ncross] = z

                    inds1 = [inds1, numpy.float(i)+alpha]
                    inds2 = [inds2, numpy.float(j)+beta]

                ncross = ncross + 1



    return result, ncross, inds2
