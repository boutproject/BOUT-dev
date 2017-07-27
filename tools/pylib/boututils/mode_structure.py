from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import range
from past.utils import old_div
import numpy as numpy
import sys
from pylab import plot,xlabel,ylim,savefig,gca, xlim, show, clf, draw, title
from boututils.fft_integrate import fft_integrate
from .ask import query_yes_no

#; Calculates mode structure from BOUT++ output
#; for comparison to ELITE
#;
#; April 2009 - Added ERGOS flag. This is intended
#;              for producing plots similar to the ERGOS
#;              vacuum RMP code

# interpolates a 1D periodic function
def zinterp( v, zind):
    
  v = numpy.ravel(v)

  nz = numpy.size(v)
  z0 = numpy.round(zind)

  p = zind - float(z0)          # between -0.5 and 0.5

  if p < 0.0 :
      z0 = z0 - 1
      p = p + 1.0
  

  z0 = ((z0 % (nz-1)) + (nz-1)) % (nz-1)

  # for now 3-point interpolation

  zp = (z0 + 1) % (nz - 1)
  zm = (z0 - 1 + (nz-1)) % (nz - 1)


  result = 0.5*p*(p-1.0)*v[zm.astype(int)] \
    + (1.0 - p*p)*v[z0.astype(int)] \
    + 0.5*p*(p+1.0)*v[zp.astype(int)]

  return result
 
    
def mode_structure( var_in, grid_in, period=1, 
                    zangle=0.0, n=None, addq=None, output=None, 
                    xq=None, xpsi=None, slow=None, subset=None, 
                    filter=None, famp=None, quiet=None, 
                    ergos=None, ftitle=None, 
                    xrange=None, yrange=None, rational=None, pmodes=None, 
                    _extra=None):


  #ON_ERROR, 2
  #
  # period = 1 ; default = full torus

    if n is None :
        if filter is not None :
            n = filter*period
        else: n = period
  

  #  if (grid_in.JYSEPS1_1 GE 0) OR (grid_in.JYSEPS1_2 NE grid_in.JYSEPS2_1) OR (grid_in.JYSEPS2_2 NE grid_in.ny-1) THEN BEGIN
  #  PRINT, "Mesh contains branch-cuts. Keeping only core"
  #  
  #  grid = core_mesh(grid_in)
  #  var = core_mesh(var_in, grid_in)
  #ENDIF ELSE BEGIN
    grid = grid_in
    vr = var_in
  #ENDELSE
  

  #IF KEYWORD_SET(filter) THEN BEGIN
  #  var = zfilter(var, filter)
  #ENDIF

    nx = grid.get('nx')
    ny = grid.get('ny')
    
    s = numpy.shape(vr)
    if numpy.size(s) != 3 :
        print("Error: Variable must be 3 dimensional")
        return
  
    if (s[0] != nx) or (s[1] != ny) :
      print("Error: Size of variable doesn't match grid")
      
      return
  
    nz = s[2]

    dz = 2.0*numpy.pi / numpy.float(period*(nz-1))
  
  # GET THE TOROIDAL SHIFT
    tn = list(grid.keys())
    tn = numpy.char.upper(tn)
    count = numpy.where(tn == "QINTY")
    if numpy.size(count) > 0 :
        print("Using qinty as toroidal shift angle")
        zShift = grid.get('qinty')
    else:
        count = numpy.where(tn == "ZSHIFT")
        if numpy.size(count) > 0 :
           print("Using zShift as toroidal shift angle")
           zShift = grid.get('zShift')
        else:
           print("ERROR: Can't find qinty or zShift variable")
           return

    zshift=grid.get('zShift')
    
    rxy=grid.get('Rxy')
    zxy=grid.get('Zxy')
    Btxy=grid.get('Btxy')
    Bpxy=grid.get('Bpxy')
    shiftangle=grid.get('ShiftAngle')
    psixy=grid.get('psixy')
    psi_axis=grid.get('psi_axis')
    psi_bndry=grid.get('psi_bndry')

    np = 4*ny

    nf = old_div((np - 2), 2)
    famp = numpy.zeros((nx, nf))

    for x in range (nx):
      #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      # transform data into fixed poloidal angle
      
      # get number of poloidal points
        nskip = numpy.zeros(ny-1)
        for y in range (ny-1):
            yp = y + 1
            nskip[y] = old_div(numpy.abs(zshift[x,yp] - zshift[x,y]), dz) - 1
      
    
        nskip =numpy.int_(numpy.round(nskip))
        nskip=numpy.where(nskip > 0, nskip, 0) 
     
       
        ny2 = numpy.int_(ny + numpy.sum(nskip)) # number of poloidal points

       # IF NOT KEYWORD_SET(quiet) THEN PRINT, x, ny2

        f = numpy.zeros(ny2)      # array for values
        R = numpy.zeros(ny2)      # Rxy
        Z = numpy.zeros(ny2)      # Zxy
        BtBp = numpy.zeros(ny2)   # Bt / Bp
      
      # interpolate values onto points
        
        ypos = 0
        for y in range(ny-1):
          # original points
            zind = old_div((zangle - zshift[x,y]),dz)
         
          
            if numpy.size(zind) != 1 : sys.exit()
            f[ypos] = zinterp(vr[x,y,:], zind)
            R[ypos] = rxy[x,y]
            Z[ypos] = zxy[x,y]
            BtBp[ypos] = old_div(Btxy[x,y], Bpxy[x,y])

            ypos = ypos + 1

          # add the extra points
          
            zi0 = old_div((zangle - zshift[x,y]),dz)
            zip1 = old_div((zangle - zshift[x,y+1]),dz)

            dzi = old_div((zip1 - zi0), (nskip[y] + 1))

            for i in range (nskip[y]):
                zi = zi0 + numpy.float(i+1)*dzi # zindex 
                w = old_div(numpy.float(i+1),numpy.float(nskip[y]+1)) # weighting
              
                f[ypos+i] = w*zinterp(vr[x,y+1,:], zi) + (1.0-w)*zinterp(vr[x,y,:], zi)
              
                R[ypos+i] = w*rxy[x,y+1] + (1.0-w)*rxy[x,y]
                Z[ypos+i] = w*zxy[x,y+1] + (1.0-w)*zxy[x,y]
                BtBp[ypos+i] = old_div((w*Btxy[x,y+1] + (1.0-w)*Btxy[x,y]), (w*Bpxy[x,y+1] + (1.0-w)*Bpxy[x,y]))
             
            ypos = ypos + nskip[y]
            
            # final point

            zind = old_div((zangle - zShift[x,ny-1]),dz)

            f[ypos] = zinterp(vr[x,ny-1,:], zind)
            R[ypos] = rxy[x,ny-1]
            Z[ypos] = zxy[x,ny-1]
            BtBp[ypos] = old_div(Btxy[x,ny-1], Bpxy[x,ny-1])
         

      #STOP

      #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      #; calculate poloidal angle

      
        drxy = numpy.gradient(R)
        dzxy = numpy.gradient(Z)
        dl = numpy.sqrt(drxy*drxy + dzxy*dzxy)
      
        nu = dl * BtBp / R # field-line pitch
        theta = old_div(numpy.real(fft_integrate(nu)), shiftangle[x])
      
        if numpy.max(theta) > 1.0 :
          # mis-match between q and nu (integration error?)
            if quiet is None : print("Mismatch  ", x, numpy.max(theta))
            theta = old_div(theta, (numpy.max(theta) + numpy.abs(theta[1] - theta[0])))
       
      
        theta = 2.0*numpy.pi * theta

      #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      #; take Fourier transform in theta angle

        tarr = 2.0*numpy.pi*numpy.arange(np) / numpy.float(np) # regular array in theta

        farr = numpy.interp(tarr, theta, f)

      #STOP

        ff = old_div(numpy.fft.fft(farr),numpy.size(farr))

        for i in range (nf):
            famp[x, i] = 2.0*numpy.abs(ff[i+1])
         
      
   
  
  # sort modes by maximum size

    fmax = numpy.zeros(nf)
    for i in range(nf):
        fmax[i] = numpy.max(famp[:,i])
     

    inds = numpy.argsort(fmax)[::-1]
 

    if pmodes is None : pmodes = 10

    qprof = old_div(numpy.abs(shiftangle), (2.0*numpy.pi))

    xarr = numpy.arange(nx)
    xtitle="Radial index"
    if xq is not None :
        # show as a function of q*n
        xarr = qprof*numpy.float(n)

        xtitle="q * n"
    elif xpsi is not None :
        # show as a function of psi. Should be normalised psi
        xarr = psixy[:,0]

        # Check if the grid includes psi axis and boundary
        count1 = numpy.where(tn == "PSI_AXIS")
        count2 = numpy.where(tn == "PSI_BNDRY")
      
        if (numpy.size(count1) > 0) and (numpy.size(count2) > 0) :
            xarr = old_div((xarr - psi_axis), (psi_bndry - psi_axis))
        
        else:
            # Use hard-wired values
            print("WARNING: Using hard-wired psi normalisation")
            # for circular case
            #xarr = (xarr + 0.1937) / (0.25044 + 0.1937)
            # for ellipse case
            #xarr = xarr / 0.74156
        
            # cbm18_dens8
            xarr = old_div((xarr + 0.854856), (0.854856 + 0.0760856))
         
      
        xtitle="Psi normalised"
     

  
    if slow is not None :
        # plot modes slowly for examination
        #safe_colors, /first
#        ax = fig.add_subplot(111)
        # go through and plot each mode
        for i in range(nf):
            if numpy.max(famp[:,i]) > 0.05*numpy.max(famp):
                print("Mode m = ", i+1, " of ", nf)
                plot(xarr, famp[:,i], 'k')
                ylim(0,numpy.max(famp))
                xlim(xrange)
                xlabel(xtitle)
                show(block=False)
                
                q = old_div(numpy.float(i+1), numpy.float(n))
        
                pos = numpy.interp(q, qprof, xarr)
                
                plot( [pos, pos],[0, 2.*numpy.max(fmax)], 'k--')
                draw()
                
                ans=query_yes_no('next mode')
                if ans:
                    clf()
            
         
    
    elif ergos is not None :
        # ERGOS - style output
    
        if output is not None and slow is None :
            savefig('output.png')
            
        
#
#        contour2, famp, xarr, indgen(nf)+1, $
#              xlabel=xtitle, xrange=xrange, yrange=yrange, _extra=_extra
#
#        ; overplot the q profile
#
#        oplot, xarr, qprof * n, color=1, thick=2
#  
#        IF KEYWORD_SET(rational) THEN BEGIN
#            maxm = FIX(MAX(qprof)) * n
#
#            qreson = (FINDGEN(maxm)+1) / FLOAT(n)
#    
#            ; get x location for each of these resonances
#            qloc = INTERPOL(xarr, qprof, qreson)
#
#            oplot, qloc, findgen(maxm)+1., psym=4, color=1
#        ENDIF
#    
#        IF KEYWORD_SET(output) THEN BEGIN
#            ; output data to save file
#            SAVE, xarr, qprof, famp, file=output+".idl"
#      
#            DEVICE, /close
#            SET_PLOT, 'X'
#        ENDIF

    else:
        if output is not None and slow is None :
            savefig('output.png')
         #  savefig('output.ps')
        
  #  
  #  
        if subset is not None :
      
            # get number of modes larger than 5% of the maximum
            count = numpy.size(numpy.where(fmax > 0.10*numpy.max(fmax)))
            
            minind = numpy.min(inds[0:count])
            maxind = numpy.max(inds[0:count])
      
            print("Mode number range: ", minind, maxind)
      
            plot( xarr, famp[:,0], 'k', visible=False)
            ylim(0,numpy.max(famp))
            xlabel(xtitle)
            xlim(xrange)
            title(ftitle)
            
            gca().set_color_cycle(['red', 'red', 'black', 'black'])
                          
            for i in range(minind, maxind+1, subset):
                plot( xarr, famp[:,i])
        
                q = old_div(numpy.float(i+1), numpy.float(n))
                pos = numpy.interp(q, qprof, xarr)
        
                plot( [pos, pos], [0, 2.*numpy.max(fmax)], '--')
        

  #    
        else:
            # default - just plot everything
             gca().set_color_cycle(['black', 'red'])
            
             plot(xarr, famp[:,0])
             ylim(0,numpy.max(famp)) #, color=1, 
             xlabel(xtitle) #, chars=1.5, xrange=xrange,title=title, _extra=_extra
             xlim(xrange)   
             for i in range (nf):
                plot( xarr, famp[:,i])
  
              
  #        
  #          IF KEYWORD_SET(addq) THEN BEGIN
  #      
  #              FOR i=0, pmodes-1 DO BEGIN
  #                  PRINT, "m = "+STRTRIM(STRING(inds[i]+1), 2)+" amp = "+STRTRIM(STRING(fmax[inds[i]]),2)
  #                  q = FLOAT(inds[i]+1) / FLOAT(n)
  #          
  #                  pos = INTERPOL(xarr, qprof, q)
  #          
  #                  oplot, [pos, pos], [0, 2.*MAX(fmax)], lines=2, color=1
  #              ENDFOR
  #          ENDIF
  #    
  #      ENDELSE
  #          IF KEYWORD_SET(output) THEN BEGIN
  #              ; output data to save file
  #              SAVE, xarr, qprof, famp, file=output+".idl"
  #    
  #              DEVICE, /close
  #              SET_PLOT, 'X'
  #          ENDIF
  #      ENDELSE
  #
