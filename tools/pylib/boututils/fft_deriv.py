from __future__ import division
try:
    from builtins import range
except:
    pass
try:
    from past.utils import old_div
except:
    import sys
    if sys.version_info[0]==2:
        def old_div(a,b):
            return a/b
    else:
        raise
import numpy

#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# FFT_DERIV: Calculates the derivative of a variable on a         ;
# periodic domain.                                                ;
#                                                                 ;
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
def fft_deriv ( var ):
    #on_error, 2
                 
    n = numpy.size(var)

    F = old_div(numpy.fft.fft(var),n)  #different definition between IDL - python
                 
    imag = numpy.complex(0.0, 1.0)
    imag = numpy.complex_(imag)

    F[0] = 0.0
          
                
    if (n % 2) == 0 :
      # even number
        for i in range (1, old_div(n,2)) :
          a = imag*2.0*numpy.pi*numpy.float(i)/numpy.float(n)
          F[i] = F[i] * a         # positive frequencies
          F[n-i] = - F[n-i] * a   # negative frequencies
           
        F[old_div(n,2)] = F[old_div(n,2)] * (imag*numpy.pi)
    else:
      # odd number
        for i in range (1, old_div((n-1),2)+1) :
          a = imag*2.0*numpy.pi*numpy.float(i)/numpy.float(n)
          F[i] = F[i] * a
          F[n-i] = - F[n-i] * a 
        
    

    result = numpy.fft.ifft(F)*n  #different definition between IDL - python
    
      
    return result 
    
