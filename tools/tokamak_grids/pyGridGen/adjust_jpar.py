from __future__ import print_function
from __future__ import division
from past.utils import old_div
# Adjusts Jpar0 to get force balance in reduced MHD models
# i.e. the equilibrium terms in the vorticity equation: 
#
# B^2 Grad_par(Jpar0/B) + 2*b0xk dot Grad(P0) = 0
#
# First finds the Pfirsch-Schluter current, then matches
# the current with the input at the inboard midplane
# 
# Usage:
# 
# IDL> adjust_jpar, "some.grid.nc"
# 
# will calculate the new Jpar, then offer to write it to the file
# (overwriting the old value of Jpar)
#
# IDL> g = file_import("some.grid.nc")
# IDL> adjust_jpar, g, jpar=jpar
#
# Given a grid structure, just calculates the new Jpar
#
# B.Dudson, May 2011
#
# 
#
import numpy
from gen_surface import gen_surface
from ddy import DDY
from ddx import  DDX
from int_y import  int_y
import sys

def grad_par( var, mesh):
  dtheta = 2.*numpy.pi / numpy.float(numpy.sum(mesh.npol))
  return (old_div(mesh.Bpxy, (mesh.Bxy * mesh.hthe))) * DDY(var, mesh)*dtheta / mesh.dy
 

def adjust_jpar( grid, smoothp=None, jpar=None, noplot=None):
  
  #type = numpy.type(grid)
  #if type == 'str' :
  #  #; Input is a string. Read in the data
  #  data = file_import(grid)
  #elif type == 'bunch.Bunch' :
  #  #; A structure, hopefully containing the grid data
    data = grid
  #else:
  #  print "ERROR: Not sure what to do with this type of grid input"
  #  return
   
  
  #; Find the inboard midplane. Use inboard since this is maximum B
  #; Matching here rather than outboard produces more realistic results
  #; (current doesn't reverse direction at edge)
    mid_ind = -1
    status = gen_surface(mesh=data) # Start generator
    while True:
        period, yi, xi, last = gen_surface(last=None, xi=None, period=None)
  
        if period :
            mid_ind = numpy.argmin(data.Rxy[xi, yi])
            out_mid = numpy.argmax(data.Rxy[xi, yi])
        break
     
        if last==1 : break
  
    if mid_ind < 0 :
        print("ERROR: No closed flux surfaces?")
        return
   
  #; Calculate 2*b0xk dot Grad P
  
    kp = 2.*data.bxcvx*DDX(data.psixy, data.pressure)
    
  #; Calculate B^2 Grad_par(Jpar0)
  
    gj = data.Bxy**2 * grad_par(old_div(data.jpar0,data.Bxy), data)

  
  #; Generate Jpar0 by integrating kp (Pfirsch-Schluter current)
  #; Grad_par = (Bp / (B*hthe))*d/dy
  
    gparj = -kp * data.hthe / (data.Bxy * data.Bpxy)
        
    
    ps = data.Bxy * int_y(gparj, data, nosmooth='nosmooth') * data.dy


  #; In core region add divergence-free parallel current to match input at
  #; inboard midplane. Using inboard as if the outboard is matched then
  #; unphysical overshoots in jpar can result on the inboard side
  #
  #; Need to make sure this bootstrap current is always in the same
  #; direction 
    dj = data.jpar0[:,mid_ind] - ps[:,mid_ind]
    ind = numpy.argmax(numpy.abs(dj))
    s = numpy.sign(dj[ind])

    w = numpy.where(dj * s < 0.0)[0] # find where contribution reverses
    if w.size > 0 : dj[w] = 0.0 # just zero in this region


    jpar = ps
    status = gen_surface(mesh=data) # Start generator
    while True:
        period, yi, xi, last = gen_surface(period=period, last=last, xi=xi)
    
        if period == None :
      # Due to multi-point differencing, dp/dx can be non-zero outside separatrix
            ps[xi,yi] = 0.0
            jpar[xi,yi] = 0.0
     

        w = numpy.size(numpy.where(yi == mid_ind))

        if (w != 0) and period != None :
      # Crosses midplane
      
            dj_b = old_div(dj[xi], data.Bxy[xi,mid_ind])
            jpar[xi,yi] = jpar[xi,yi] + dj_b * data.Bxy[xi,yi]
     
        if last==1 : break
  
  

        
 # if noplot!=None :
  #  WINDOW, xsize=800, ysize=800
  #  !P.multi=[0,2,2,0,0]
  #  SURFACE, data.jpar0, tit="Input Jpar0", chars=2
  #  SURFACE, jpar, tit="New Jpar0", chars=2
  #  PLOT, data.jpar0[0,*], tit="jpar at x=0. Solid=input", yr=[MIN([data.jpar0[0,*],jpar[0,*]]), $
  #                                                             MAX([data.jpar0[0,*],jpar[0,*]])]
  #  OPLOT, jpar[0,*], psym=1
  #
  #  #;x = data.ixseps1-1
  #  #;PLOT, data.jpar0[x,*], tit="Jpar at x="+STR(x)+" Solid=input", $
  #  #;  yr=[MIN([data.jpar0[x,*],jpar[x,*]]), $
  #  #;      MAX([data.jpar0[x,*],jpar[x,*]])]
  #  #;OPLOT, jpar[x,*], psym=1
  #  
  #  y = out_mid
  #  PLOT, data.jpar0[*,y], tit="Jpar at y="+STR(y)+" Solid=input", $
  #    yr=[MIN([data.jpar0[*,y],jpar[*,y]]), $
  #        MAX([data.jpar0[*,y],jpar[*,y]])]
  #  OPLOT, jpar[*,y], psym=1
    
    return jpar
   
  
  #if type == 'str' :
  # # Ask if user wants to write this new Jpar to file
  #  
  #  if query_yes_no("Write new Jpar0 to file?") :
  #    
  #    f = file_open(grid, /write) # Read/write mode
  #    
  #    status = file_write(f, "Jpar0", jpar)
  #    
  #    if status :
  #      print "ERROR writing Jpar0 to file '"+ grid+"'"
  #     
  #    
  #    file_close, f
     
   
