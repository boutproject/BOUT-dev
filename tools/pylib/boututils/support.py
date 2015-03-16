from __future__ import division
from builtins import range
from past.utils import old_div
import numpy
from scipy.integrate import simps
import copy


def deriv(y, x=None):
    
    n=numpy.size(y)
    dy = numpy.zeros(n)
    
    if x==None :
  
        dy[0] = old_div((-3*y[0] + 4*y[1] - y[2]), 2)
    
        for i in range (1,n-1):        
            dy[i] = old_div((y[i+1] - y[i-1]), 2) #where i = 1...N-2

        
        dy[n-1] = old_div((3*y[n-1] - 4*y[n-2] + y[n-3]), 2)


    else:
        
        x12 = -x[1] + x[0]   #x1 - x2
        x01 = -x[2] + x[1]    #x0 - x1
        x02 = -x[2] + x[0] #x0 - x2

        dy[0] = (
                y[0] * (x01+x02)/(x01*x02) -  #First point
                y[1] * x02/(x01*x12) +   
                y[2] * x01/(x02*x12)
                )
                
        for i in range (1,n-1):  
             
             
            x12 = x[i] - x[i-1]   #x1 - x2
            x01 = x[i+1] - x[i]    #x0 - x1
            x02 = x[i+1] - x[i-1] #x0 - x2
            
            
            dy[i] = (
                    y[i+1] * (old_div(x12, (x01*x02))) +   #Middle points
                    y[i] * (old_div(1.,x12) - old_div(1.,x01)) -  
                    y[i-1] * (old_div(x01, (x02 * x12)))
                    )
                    
                    
        x12 = -x[n-2] + x[n-3]   #x1 - x2
        x01 = -x[n-1] + x[n-2]    #x0 - x1
        x02 = -x[n-1] + x[n-3] #x0 - x2
            
            
            
        dy[n-1] = (
                    -y[n-3] * x12/(x01*x02) + #Last point
                    y[n-2] * x02/(x01*x12) -  
                    y[n-1] * (x02+x12) / (x02*x12)
                    )

    return dy

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
    if fin == None :
        f = copy.deepcopy(xin)
        x = numpy.arange(numpy.size(f)).astype(float)
    else:
        f = copy.deepcopy(fin)
        x = copy.deepcopy(xin)
    
    n = numpy.size(f)

    g = numpy.zeros(n)

    if simple != None :
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
    
