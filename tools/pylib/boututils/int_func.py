from __future__ import division
from builtins import range
from past.utils import old_div
import numpy
from scipy.integrate import simps
import copy

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

def int_func( xin, fin=None, simple=None):
    if fin is None :
        f = copy.deepcopy(xin)
        x = numpy.arange(numpy.size(f)).astype(float)
    else:
        f = copy.deepcopy(fin)
        x = copy.deepcopy(xin)
    
    n = numpy.size(f)

    g = numpy.zeros(n)

    if simple is not None :
     # Just use trapezium rule
     
        g[0] = 0.0
        for i in range (1, n) :
            g[i] = g[i-1] + 0.5*(x[i] - x[i-1])*(f[i] + f[i-1])
         
    else:
     
        n2 = numpy.int(old_div(n,2))
     
        g[0] = 0.0
        for i in range (n2, n) :
            g[i] = simps( f[0:i+1], x[0:i+1])


            
        for i in range (1, n2) :
            g[i] = g[n-1] - simps( f[i::], x[i::])
             
    return g
 
 
