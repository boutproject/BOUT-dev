from __future__ import division
from past.utils import old_div
import numpy
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
    if include_in is None :
        x = x + 0.5
        m = m + 0.5
  
    
    if include_out is None:
        m = m + 0.5
 
    x = old_div(x, m) 
  

    if in_dp is None and out_dp is None :
    # Neither inner or outer gradients set. Just return equal spacing
        return pin + (pout - pin)*x
  
    
    norm = (x[1] - x[0])*(pout - pin)

    if in_dp is not None and out_dp is not None :
    # Fit to dist = a*i^3 + b*i^2 + c*i
        c = old_div(in_dp,norm)
        b = 3.*(1. - c) - old_div(out_dp,norm) + c
        a = 1. - c - b
    elif in_dp is not None :
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
