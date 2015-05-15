from __future__ import division
from builtins import range
from past.utils import old_div
#; Monotone interpolation using Hermite splines
#; 
#; x, y - input points
#; u - points where result is needed
import numpy



def h00( t ):
  return (1. + 2.*t)*(1. - t)^2
 

def h10( t ):
  return t*(1.-t)^2
 

def h01( t ):
  return t**2*(3. - 2.*t)
 

def h11( t ):
  return t**2*(t - 1.)
 

def spline_mono( x, y, u, yp0=None, ypn_1=None):
  n = numpy.size(x)

  if n < 3 :
    # Just linear interpolate
    return numpy.interp(u, x, y)
   

  # Calculate delta
  D = old_div((y[1::] - y[0:(n-1)]), (x[1::] - x[0:(n-1)]))
  
  if numpy.size(yp0) == 0 : yp0 = D[0]
  if numpy.size(ypn_1) == 0 : ypn_1 = D[n-2]
  
  m = numpy.zeros(n)
  m[1:(n-2)] = 0.5*(D[0:(n-3)] + D[1:(n-2)])
  m[0] = yp0
  
  for i in range (n-1) :
    if numpy.abs(D[i]) < 1.e-6 :
      m[i] = 0.
      m[i+1 < (n-1)] = 0.
    else:
      a = old_div(m[i], D[i])
      b = old_div(m[i+1], D[i])
      
      c = numpy.sqrt(a**2 + b**2)
      if c > 3. :
        t = old_div(3., c)
        m[i] = t * a * D[i]
        m[i+1] = t* b * D[i]
  
  nout = numpy.size(u)
  
  result = numpy.zeros(nout)
  
  for i in range (nout) :
    count=numpy.sum(numpy.where(x > u[i]))
    xup = numpy.min(count)
    
    
    if count < 0 :
      result[i] = y[n-1]
       
     
    if xup == 0 :
      result[i] = y[0]
       
     
    xlow = xup - 1
    
    h = numpy.float(x[xup] - x[xlow])
    t = old_div((numpy.float(u[i]) - numpy.float(x[xlow])), h)
    
    result[i] = (
                y[xlow] * h00(t) +  
                h*m[xlow]*h10(t) +  
                y[xup]*h01(t) +  
                h*m[xup]*h11(t)
                )
  
  return result
