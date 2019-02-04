"""
Derivatives and integrals of periodic and non-periodic functions


B.Dudson, University of York, Nov 2009
"""
from __future__ import print_function
from __future__ import division

from builtins import range

try:
    from past.utils import old_div
except ImportError:
    def old_div(a, b):
        return a / b

from numpy import zeros, pi, array, transpose, sum, where, arange, multiply
from numpy.fft import rfft, irfft

def deriv(*args, **kwargs):
    """Take derivative of 1D array

    result = deriv(y)
    result = deriv(x, y)

    keywords

    periodic = False    Domain is periodic
    """

    nargs = len(args)
    if nargs == 1:
        var = args[0]
        x = arange(var.size)
    elif nargs == 2:
        x = args[0]
        var = args[1]
    else:
        raise RuntimeError("deriv must be given 1 or 2 arguments")

    try:
        periodic = kwargs['periodic']
    except:
        periodic = False

    n = var.size
    if periodic:
        # Use FFTs to take derivatives
        f = rfft(var)
        f[0] = 0.0 # Zero constant term
        if n % 2 == 0:
            # Even n
            for i in arange(1,old_div(n,2)):
                f[i] *= 2.0j * pi * float(i)/float(n)
            f[-1] = 0.0 # Nothing from Nyquist frequency
        else:
            # Odd n
            for i in arange(1,old_div((n-1),2) + 1):
                f[i] *= 2.0j * pi * float(i)/float(n)
        return irfft(f)
    else:
        # Non-periodic function
        result = zeros(n) # Create empty array
        if n > 2:
            for i in arange(1, n-1):
                # 2nd-order central difference in the middle of the domain
                result[i] = old_div((var[i+1] - var[i-1]), (x[i+1] - x[i-1]))
            # Use left,right-biased stencils on edges (2nd order)
            result[0]   = old_div((-1.5*var[0]   + 2.*var[1]   - 0.5*var[2]), (x[1] - x[0]))
            result[n-1] =  old_div((1.5*var[n-1] - 2.*var[n-2] + 0.5*var[n-3]), (x[n-1] - x[n-2]))
        elif n == 2:
            # Just 1st-order difference for both points
            result[0] = result[1] = old_div((var[1] - var[0]),(x[1] - x[0]))
        elif n == 1:
            result[0] = 0.0
        return result

def deriv2D(data,axis=-1,dx=1.0,noise_suppression=True):
  """ Takes 1D or 2D Derivative of 2D array using convolution
	
	result = deriv2D(data)
	result = deriv2D(data, dx)
	
	output is 2D (if only one axis specified)
	output is 3D if no axis specified [nx,ny,2] with the third dimension being [dfdx, dfdy]
	
	keywords:
	axis = 0/1  If no axis specified 2D derivative will be returned
	dx = 1.0    axis spacing, must be 2D if 2D deriv is taken - default is [1.0,1.0]
	noise_suppression = True   noise suppressing coefficients used to take derivative - default = True
  """

  from scipy.signal import convolve
  
  s = data.shape
  if axis > len(s)-1:
    raise RuntimeError("ERROR: axis out of bounds for derivative")

  if noise_suppression:
    if s[axis] < 11:
      raise RuntimeError("Data too small to use 11th order method")
    tmp = array([old_div(-1.0,512.0),old_div(-8.0,512.0),old_div(-27.0,512.0),old_div(-48.0,512.0),old_div(-42.0,512.0),0.0,old_div(42.0,512.0),old_div(48.0,512.0),old_div(27.0,512.0),old_div(8.0,512.0),old_div(1.0,512.0)])
  else:
    if s[axis] < 9:
      raise RuntimeError("Data too small to use 9th order method")
    tmp = array([old_div(1.0,280.0),old_div(-4.0,105.0),old_div(1.0,5.0),old_div(-4.0,5.0),0.0,old_div(4.0,5.0),old_div(-1.0,5.0),old_div(4.0,105.0),old_div(-1.0,280.0)])
    
  N = int((tmp.size-1)/2)
  if axis==1:
    W = transpose(tmp[:,None])
    data_deriv = convolve(data,W,mode='same')/dx*-1.0
    for i in range(s[0]):
      data_deriv[i,0:N-1] = old_div(deriv(data[i,0:N-1]),dx)
      data_deriv[i,s[1]-N:] = old_div(deriv(data[i,s[1]-N:]),dx)

  elif axis==0:
    W = tmp[:,None]
    data_deriv = convolve(data,W,mode='same')/dx*-1.0
    for i in range(s[1]):
      data_deriv[0:N-1,i] = old_div(deriv(data[0:N-1,i]),dx)
      data_deriv[s[0]-N:,i] = old_div(deriv(data[s[0]-N:,i]),dx)
  else:
    data_deriv = zeros((s[0],s[1],2))
    if (not hasattr(dx, '__len__')) or len(dx)==1:
      dx = array([dx,dx])

    W = tmp[:,None]#transpose(multiply(tmp,ones((s[1],tmp.size))))
    data_deriv[:,:,0] = convolve(data,W,mode='same')/dx[0]*-1.0
    for i in range(s[1]):
      data_deriv[0:N-1,i,0]  =  old_div(deriv(data[0:N-1,i]),dx[0])
      data_deriv[s[0]-N:s[0]+1,i,0] = old_div(deriv(data[s[0]-N:s[0]+1,i]),dx[0])

    W = transpose(tmp[:,None])#multiply(tmp,ones((s[0],tmp.size)))
    data_deriv[:,:,1] = convolve(data,W,mode='same')/dx[1]*-1.0
    for i in range(s[0]):
      data_deriv[i,0:N-1,1] = old_div(deriv(data[i,0:N-1]),dx[1])
      data_deriv[i,s[1]-N:s[1]+1,1] = old_div(deriv(data[i,s[1]-N:s[1]+1]),dx[1])

  return data_deriv

def integrate(var, periodic=False):
    """Integrate a 1D array

    Return array is the same size as the input
    """
    if periodic:
        # Use FFT
        f = rfft(var)
        n = var.size
        # Zero frequency term
        result = f[0].real*arange(n, dtype=float)
        f[0] = 0.
        if n % 2 == 0:
            # Even n
            for i in arange(1,old_div(n,2)):
                f[i] /= 2.0j * pi * float(i)/float(n)
            f[-1] = 0.0 # Nothing from Nyquist frequency
        else:
            # Odd n
            for i in arange(1,old_div((n-1),2) + 1):
                f[i] /= 2.0j * pi * float(i)/float(n)
        return result + irfft(f)
    else:
        # Non-periodic function
        def int_total(f):
            """Integrate over a set of points"""
            n = f.size
            if n > 7:
                # Need to split into several segments
                # Use one 5-point, leaving at least 4-points
                return int_total(f[0:5]) + int_total(f[4:])
            elif (n == 7) or (n == 6):
                # Try to keep 4th-order
                # Split into 4+4 or 4+3
                return int_total(f[0:4]) + int_total(f[3:])
            elif n == 5:
                # 6th-order Bool's rule
                return 4.*(7.*f[0] + 32.*f[1] + 12.*f[2] + 32.*f[3] + 7.*f[4])/90.
            elif n == 4:
                # 4th-order Simpson's 3/8ths rule
                return 3.*(f[0] + 3.*f[1] + 3.*f[2] + f[3])/8.
            elif n == 3:
                # 4th-order Simpson's rule
                return (f[0] + 4.*f[1] + f[2])/3.
            elif n == 2:
                # 2nd-order Trapezium rule
                return 0.5*(f[0] + f[1])
            else:
                print("WARNING: Integrating a single point")
                return 0.0
        # Integrate using maximum number of grid-points
        n = var.size
        n2 = int(old_div(n,2))
        result = zeros(n)
        for i in arange(n2, n):
            result[i] = int_total(var[0:(i+1)])
        for i in arange(1, n2):
            result[i] = result[-1] - int_total(var[i:])
        return result

def simpson_integrate(data,dx,dy,kernel=0.0,weight=1.0):
  """ Integrates 2D data to one value using the simpson method and matrix convolution

      result = simpson_integrate(data,dx,dy)

      keywords:

      kernel - can be supplied if the simpson matrix is calculated ahead of time
	    - if not supplied, is calculated within this function
	    - if you need to integrate the same shape data over and over, calculated
		it ahead of time using:
		  kernel = simpson_matrix(Nx,Ny,dx,dy)

      weight - can be used to scale data if single number
	    - can be used to mask data if weight is array (same size as data)
  """
  s = data.shape
  Nx = s[0]
  Ny = s[1]

  if len(kernel)==1:
    kernel = simpson_matrix(Nx,Ny,dx,dy)

  return sum(multiply(multiply(weight,kernel),data))/sum(multiply(weight,kernel))


def simpson_matrix(Nx,Ny,dx,dy):
  """
      Creates a 2D matrix of coefficients for the simpson_integrate function

      Call ahead of time if you need to perform integration of the same size data with the
	same dx and dy
      
      Otherwise, simpson_integrate will automatically call this

  """
  Wx = arange(Nx) + 2
  Wx[where(arange(Nx) % 2 == 1)] = 4
  Wx[0] = 1
  Wx[Nx-1] = 1

  Wy = arange(Ny) + 2
  Wy[where(arange(Ny) % 2 == 1)] = 4
  Wy[0] = 1
  Wy[Ny-1] = 1

  W = Wy[None,:] * Wx[:,None]

  A = dx*dy/9.0

  return W*A
