"""

Transform data for input to BOUT++

"""
from builtins import range

from numpy.fft import rfft
from numpy import ndarray

def transform3D(arr):
    """
    Transforms a 3D array. BOUT++ expects 3D inputs
    to be Fourier transformed in the Z direction.

    Inputs
    ======
    
    arr   - 3D array [x,y,z]


    Returns
    =======

    A 3D array [x,y,kz]

    where kz is organised in the standard FFT order,
    with constant (DC, kz=0) component first, followed by
    real/imaginary pairs.

    kz = [0, (real, imag), (real, imag), ...]
    """
    
    if len(arr.shape) != 3:
        raise ValueError("Input array must be 3D")

    # Take FFT over z (last index), returning a complex array
    fa = rfft(arr, axis=-1)

    nmodes = fa.shape[-1]

    # scipy fft normalises to N, but fftw doesn't
    fa /= arr.shape[-1]
    # Unpack complex array into a real array

    shape = list(arr.shape)
    shape[-1] = 1 + (nmodes-1)*2 # One for DC + 2 for other modes

    result = ndarray(shape)

    # kz = 0 (DC) component only has real part
    result[:,:,0] = fa[:,:,0].real

    # All other components have both real and imaginary parts
    for k in range(1,nmodes):
        result[:,:,2*k-1] = fa[:,:,k].real
        result[:,:,2*k] = fa[:,:,k].imag

    return result

