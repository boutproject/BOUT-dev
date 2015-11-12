from builtins import range
import numpy
from boututils.calculus import deriv

# calculates x (psi) derivative for 2D variable
def DDX( psi, var):
  s = numpy.shape(var)
  nx = s[0]
  ny = s[1] 
  
  
  dv = numpy.zeros((nx, ny))
 
  for i in range (ny) :
      
      dv[:,i] = deriv( psi[:,i], var[:,i] )

   
  return dv
