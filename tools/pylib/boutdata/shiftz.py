from numpy import ndarray, pi, cos, sin
from numpy import fft


def shiftz(var, zangle, zperiod=1.0):
    """Shift a variable in Z, changing between field-aligned and
    orthogonal X-Z coordinates. This mainly used for tokamak
    simulations in field-aligned coordinates.

    Parameters
    ----------
    var : array_like
        Data to be shifted
            4D [t,x,y,z]
            3D [x,y,z] or [t,x,z]
            2D [x,z]
    zangle : array_like
        The shift angle
            2D [x,y]  (if var is 4D or 3D [x,y,z])
            1D [x]    (if var is 3D [t,x,z] or 2D)
    zperiod : float, optional
        The fraction of 2pi covered by the variable in Z. This
        corresponds to the ZPERIOD variable in BOUT.inp and multiplies
        the kz wavenumbers by this factor.

    Returns
    -------
    ndarray
        A numpy array of the same size and shape as var

    Examples
    --------

    >>> from boutdata import collect
    >>> from boututils.datafile import DataFile
    >>> from boutdata.shiftz import shiftz
    >>> n = collect("Ne")  # Read 4D variable [t,x,y,z]
    >>> d = DataFile("grid.nc")    # Read the grid file
    >>> nxz = shiftz(n, d["zShift"], zperiod=4)

    nxz is now in orthogonal X-Z coordinates (X is psi).

    Note that in older grid files "qinty" is used rather
    than "zShift".

    """
    
    if len(var.shape) == 4:
        # 4D variable [t,x,y,z]
        result = ndarray(var.shape)
        for t in range(var.shape[0]):
            # Shift each time slice separately
            result[t,:,:,:] = shiftz(var[t,:,:,:], zangle, zperiod=zperiod)
        return result
    elif len(var.shape) == 3:
        if len(zangle.shape) == 2:
            # 3D variable [x,y,z], array [x,y]
            result = ndarray(var.shape)
            for y in range(var.shape[1]):
                result[:,y,:] = shiftz(var[:,y,:], zangle[:,y], zperiod=zperiod)
            return result
        elif len(zangle.shape) == 1:
            # 3D variable [t,x,z], array [x]
            result = ndarray(var.shape)
            for t in range(var.shape[0]):
                result[t,:,:] = shiftz(var[t,:,:], zangle, zperiod=zperiod)
            return result
        else:
            raise ValueError("Expecting zangle to be 1 or 2D")
    elif len(var.shape) == 2:
        if len(zangle.shape) != 1:
            raise ValueError("Expecting zangle to be 1D")
        
        ################################
        # Main algorithm here
        # var is [x,z]
        # zangle is [x]
        
        # Take FFT in Z direction
        f = fft.rfft(var, axis=1)
        
        zlength = 2.*pi/zperiod
        
        for z in range(1, f.shape[1]):
            kwave=z*2.0*pi/zlength
            f[:,z] *= cos(kwave * zangle) - 1j*sin(kwave*zangle)
        return fft.irfft(f, var.shape[1], axis=1)
        
    else:
        raise ValueError("Don't know how to handle 1D variable")
    
