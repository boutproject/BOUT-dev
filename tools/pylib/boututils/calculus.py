"""
Derivatives and integrals of periodic and non-periodic functions


B.Dudson, University of York, Nov 2009
"""

try:
    from numpy import zeros, arange, pi
    from numpy.fft import rfft, irfft
except ImportError:
    print "ERROR: NumPy module not available"
    raise


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
            for i in arange(1,n/2):
                f[i] *= 2.0j * pi * float(i)/float(n)
            f[-1] = 0.0 # Nothing from Nyquist frequency
        else:
            # Odd n
            for i in arange(1,(n-1)/2 + 1):
                f[i] *= 2.0j * pi * float(i)/float(n)
        return irfft(f)
    else:
        # Non-periodic function
        result = zeros(n) # Create empty array
        if n > 2:
            for i in arange(1, n-1):
                # 2nd-order central difference in the middle of the domain
                result[i] = (var[i+1] - var[i-1]) / (x[i+1] - x[i-1])
            # Use left,right-biased stencils on edges (2nd order)
            result[0]   = (-1.5*var[0]   + 2.*var[1]   - 0.5*var[2]) / (x[1] - x[0])
            result[n-1] =  (1.5*var[n-1] - 2.*var[n-2] + 0.5*var[n-3]) / (x[n-1] - x[n-2])
        elif n == 2:
            # Just 1st-order difference for both points
            result[0] = result[1] = (var[1] - var[0])/(x[1] - x[0])
        elif n == 1:
            result[0] = 0.0
        return result


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
            for i in arange(1,n/2):
                f[i] /= 2.0j * pi * float(i)/float(n)
            f[-1] = 0.0 # Nothing from Nyquist frequency
        else:
            # Odd n
            for i in arange(1,(n-1)/2 + 1):
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
                print "WARNING: Integrating a single point"
                return 0.0
        # Integrate using maximum number of grid-points
        n = var.size
        n2 = int(n/2)
        result = zeros(n)
        for i in arange(n2, n):
            result[i] = int_total(var[0:(i+1)])
        for i in arange(1, n2):
            result[i] = result[-1] - int_total(var[i:])
        return result

