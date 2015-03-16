from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
#; Fast 2D Discrete Cosine Transform
#;
#; http://fourier.eng.hmc.edu/e161/lectures/dct/node2.html
#; 
#; Ben Dudson, University of York, Feb 2010
#; 
#; NOTE: SOMETHING NOT QUITE RIGHT HERE

import numpy
  
def DCT2D ( sig, inverse=None):
    s = numpy.shape(sig)
    if numpy.size(s) != 2 :
        print("ERROR: input to DCT2Dfast must be 2D")
        return 0
   
    nx = s[0]
    ny = s[1]
  
    result = numpy.size((nx, ny))
    
    for i in range (ny) :
        result[:,i] = DCT(sig[:,i], inverse=inverse)
   
    for i in range (nx) :
        result[i,:] = DCT(result[i,:], inverse=inverse)
    
  
    if inverse == None :
        result = result * 2. * numpy.sqrt(nx*ny)
    else:
        result = old_div(result, (2.* numpy.sqrt(nx*ny)))
   
      
    return result
    

