#A collection of functions related to the shifting of
#fields for both twist-shift and shifted-x-derivs
#
#Author : David Dickinson
#Date : 25/03/15

try:
    from numpy import cos, sin, zeros, product, complex
except:
    print("Error: Couldn't import from numpy.")
    raise ImportError

try:
    from scipy.fftpack import fft, ifft
except:
    try:
        from numpy.fft import fft, ifft
    except:
        print("Error: Couldn't import fft routines from either scipy or numpy.")
        raise ImportError

def shiftZSlice(zSlice=None,zAngle=None,fac=None,forward=True,inPlace=False,**kwargs):
    """ 
    Shifts a single z slice using zangle.
    
    Inputs:
             zSlice  = 1d numpy array slice of Field to shift
             zAngle  = Phase angle to use in shift
             fac     = Wave-number factor (quanta) = 2*pi/zlength which is usually just zPeriod.
             forward = Determines direction of transform
             inPlace = If true then returned shifted slice also stored in zSlice
      
    Output:
        returns the shifted slice and optionally stores in zSlice (see inPlace).
    """

    #Handle inputs
    if zSlice is None:
        print("Error: Must pass zSlice.")
        return

    if zAngle is None:
        print("Error: Must pass zAngle.")
        return

    if fac is None:
        print("Error: Must pass fac.")

    if forward:
        dirMult=-1
    else:
        dirMult=1

    #Get parameters
    nz = zSlice.shape[-1] #Should check this is even
    nkz = 1+(nz/2)
    nrest = 1

    #Backup shape
    origShape=zSlice.shape
        
    #Check if we have more than one dimension, we should be able to cope with multi-dimensional fields
    #but this hasn't really been tested so print a warning.
    if len(zSlice.shape) != 1:
        print("Warning: Expecting zSlice to be 1D, proceeding anyway (applying shift to last axis) but be careful not well tested.")

        #Update number of other elements
        nrest = product(origShape[:-1])

    #/Flatten array
    fldIn = zSlice.reshape((nrest,nz))
    
    #/Forward transform
    fft_coef = fft(fldIn,n=nz,axis=-1)

    #/Shift : Note we skip jz=0 (phase==1) and jz=nkz-1 handled separately 
    for jother in xrange(nrest):
        for jz in xrange(1,nkz-1):
            kwave=jz*fac
            phase =complex(cos(kwave*zAngle),dirMult*sin(kwave*zAngle))
            fft_coef[jother][jz]  *= phase
            fft_coef[jother][-jz] *= phase.conjugate()

        #-- Handle last wavenumber
        jz=nkz-1
        kwave=jz*fac
        phase =complex(cos(kwave*zAngle),dirMult*sin(kwave*zAngle))
        fft_coef[jother][jz]  *= phase

    #/Now transform back
    shifted = ifft(fft_coef,n=nz,axis=-1)

    #/Restore the correct shape
    shifted = shifted.reshape(origShape)

    #/Update zSlice if requested
    if inPlace:
        zSlice = shifted
        
    #Return
    return shifted

def shiftZ(fld=None,zShift=None,zPeriod=None,inPlace=False,**kwargs):
    """ 
    Shifts a field using passed shift angle array.
    
    Inputs:
             fld     = 4d/3d numpy array representing field ([t,x,y,z]/[x,y,z])
             zShift  = 2d/1d numpy array of phase angle to use in shift ([x,y],[x])
             inPlace = If true then shifted field overwrites input fld.
    Output:
        returns the shifted field and optionally stores in fld (see inPlace).
    """

    #Handle inputs
    if fld is None:
        print("Error: Must pass fld.")
        raise ValueError
    if zShift is None:
        print("Error: Must pass zShift.")
        raise ValueError
    if zPeriod is None:
        print("Error: Must pass zPeriod.")
        raise ValueError

    #Now lets try to work out what we've been passed and what we want to achieve
    ndimFld=len(fld.shape)
    ndimSft=len(zShift.shape)

    #/Debug output -- can be removed following testing, though probably want to retain error checking
    print("Provided field has "+str(ndimFld)+" dimensions.")
    if ndimFld == 4:
        print(" -- Assuming dims are [t,x,y,z]")
    elif ndimFld == 3:
        print(" -- Assuming dims are [x,y,z]")
    else:
        print("Error: Currently require field to be either 4D or 3D.")
        raise IndexError #Not really index error, but could lead to it!
    print("Provided shift has "+str(ndimSft)+" dimensions.")        
    if ndimSft == 2:
        print(" -- Assuming dims are [x,y] --> This is usually for ShiftXDerivs.")
    elif ndimSft == 1:
        print(" -- Assuming dims are [x]   --> This is usually for twist-shift.")
    else:
        print("Error: Currently require shift to be either 2D or 1D.")
        raise IndexError #Not really index error, but could lead to it!


    #Should really check that nx_fld == nx_shift, for example.

    #Now we'd like to generalise the rest of the code so we reshape things
    #/First backup field shape
    origShape=fld.shape

    #/Now setup sizes
    if ndimFld == 4:
        nt = origShape[0]
        start = 1
    else:
        nt = 1
        start = 0
    
    nx = origShape[start] ; ny = origShape[start+1] ; nz = origShape[start+2]
    fldIn = fld.reshape([nt,nx,ny,nz])

    #/Now ensure zShift 2d
    if ndimSft == 2:
        zShiftIn = zShift.copy()
    else:
        zShiftIn = zeros([nx,ny])
        for jy in xrange(ny):
            zShiftIn[:,jy] = zShift[:]

    #/Now we are definitely operating on 4D field with 2D shift
    #--> Make output field
    fldOut = fldIn.copy()

    #--> Now actually do the shift
    for jt in xrange(nt):
        for jx in xrange(nx):
            for jy in xrange(ny):
                fldOut[jt,jx,jy,:] = shiftZSlice(fldIn[jt,jx,jy,:],zShiftIn[jx,jy],
                                                 fac=zPeriod,inPlace=inPlace,**kwargs)

    #/Now restore field shape
    fldOut = fldOut.reshape(origShape)

    #/Copy field if we want an in place transform
    if inPlace:
        fld = fldOut

    #/Return
    return fldOut
