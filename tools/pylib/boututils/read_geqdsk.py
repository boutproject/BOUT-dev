from __future__ import print_function
from builtins import range
import numpy
from geqdsk import Geqdsk
from boututils.bunch import Bunch

def read_geqdsk (file):
    
    data=Geqdsk()
    
    data.openFile(file)

    nxefit =data.get('nw')
    nyefit =data.get('nh')
    xdim   =data.get('rdim')
    zdim   =data.get('zdim')
    rcentr =data.get('rcentr')
    rgrid1 =data.get('rleft')
    zmid   =data.get('zmid')

    rmagx  =data.get('rmaxis')
    zmagx  =data.get('zmaxis')
    simagx =data.get('simag')
    sibdry =data.get('sibry')
    bcentr =data.get('bcentr')
    
    cpasma =data.get('current')
    #simagx =data.get('simag')
    #xdum   =data.get()
    #rmagx  =data.get('rmaxis')
    #xdum   =data.get()

    #zmagx  =data.get('zmaxis')
    #xdum   =data.get()
    #sibdry =data.get('sibry')
    #xdum   =data.get()
    #xdum   =data.get()

# Read arrays

    fpol=data.get('fpol')
    pres=data.get('pres')

    f=data.get('psirz')
    qpsi=data.get('qpsi')

    nbdry=data.get('nbbbs')
    nlim=data.get('limitr')
    
    if(nlim != 0) :
        xlim=data.get('rlim')
        ylim=data.get('zlim')
    else:
        xlim=[0]
        ylim=[0]
        
    rbdry=data.get('rbbbs')
    zbdry=data.get('zbbbs')


    # Reconstruct the (R,Z) mesh
    r=numpy.zeros((nxefit, nyefit), numpy.float64)
    z=numpy.zeros((nxefit, nyefit), numpy.float64)


    for i in range(0,nxefit):
        for j in range(0,nyefit):
            r[i,j] = rgrid1 + xdim*i/(nxefit-1)
            z[i,j] = (zmid-0.5*zdim) + zdim*j/(nyefit-1)

    f=f.T
    
    print('nxefit =  ', nxefit, '  nyefit=  ', nyefit)

    return Bunch(nx=nxefit, ny=nyefit, # Number of horizontal and vertical points
            r=r, z=z,   # Location of the grid-points
            xdim=xdim, zdim=zdim, # Size of the domain in meters
            rcentr=rcentr, bcentr=bcentr, # Reference vacuum toroidal field (m, T)
            rgrid1=rgrid1, # R of left side of domain
            zmid=zmid, # Z at the middle of the domain
            rmagx=rmagx, zmagx=zmagx, # Location of magnetic axis
            simagx=simagx, # Poloidal flux at the axis (Weber / rad)
            sibdry=sibdry, # Poloidal flux at plasma boundary (Weber / rad)
            cpasma=cpasma, #
            psi=f, # Poloidal flux in Weber/rad on grid points
            fpol=fpol, # Poloidal current function on uniform flux grid
            pres=pres, # Plasma pressure in nt/m^2 on uniform flux grid
            qpsi=qpsi, # q values on uniform flux grid
            nbdry=nbdry, rbdry=rbdry, zbdry=zbdry, # Plasma boundary
            nlim=nlim, xlim=xlim, ylim=ylim) # Wall boundary
