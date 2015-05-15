from __future__ import division
from builtins import range
from past.utils import old_div
#; 2D Discrete Cosine Transform 
#; 
#; Author: Maxim Umansky, LLNL
#;
#; 18-02-2010 Ben Dudson <bd512@york.ac.uk>
#;      * Modified to accept non-square arrays
#;      * Speed up using matrix operations instead of FOR loops
#;
#; 2011/06/02 Ilon Joseph
#;       * Double precision for forward and inverse transformation
#;       * Matrix operations for inverse
#;       * Forward # Inverse = 1 with errors on the order of 1.0d-14 of maximum

import numpy



def dct2dslow( sig, fsig, inverse=None):
#;
#;-calculate 2D discrete cosine transform
#;----------------------------------------

    if inverse==None : #---direct transform---

        s = numpy.size(sig)

        nx=s[0]
        ny=s[1]

        fsig=numpy.zeros((nx,ny))
  
        for iu in range (nx) :
            for jv in range (ny) :
      
                fsig[iu,jv] = numpy.sum( numpy.float_(sig) * 
               numpy.dot( numpy.cos(iu*numpy.float_(numpy.pi)*(2*numpy.arange(nx).astype(float)+1)/(2*nx)) , numpy.cos(jv*numpy.float_(numpy.pi)*(2*numpy.arange(ny).astype(float)+1)/(2*ny))) )
      

        fsig *= old_div(2,numpy.sqrt(numpy.float_(nx*ny)))
        fsig[0,:] *= numpy.sqrt(numpy.float_(0.5))
        fsig[:,0] *= numpy.sqrt(numpy.float_(0.5))
        
        return fsig

    else : #---inverse transform---
  
        s = numpy.shapeE(fsig)

        nx=s[0]
        ny=s[1]

        dsig=numpy.float_(fsig)
        dsig[0,:] *= numpy.sqrt(numpy.float_(0.5))
        dsig[:,0] *= numpy.sqrt(numpy.float_(0.5))

        sig=numpy.zeros((nx,ny))

        for ix in range (nx) :
            for jy in range (ny) :

                sig[ix,jy]= numpy.sum( dsig *  
                   numpy.dot( numpy.cos(numpy.arange(nx).astype(float)*numpy.float_(numpy.pi)*(2*ix+1)/(2*nx)) , numpy.cos(numpy.arange(ny).astype(float)*numpy.float_(numpy.pi)*(2*jy+1)/(2*ny)) ) )


        sig *= old_div(2,numpy.sqrt(numpy.float_(nx*ny)))


        return sig