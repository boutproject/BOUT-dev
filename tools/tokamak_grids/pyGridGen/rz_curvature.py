from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
from boututils.bunch import Bunch
import numpy
from boututils.calculus import deriv
from create_grid import local_gradient
from scipy import interpolate

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

def dct2dslow( fun, inverse=None):
#;
#;-calculate 2D discrete cosine transform
#;----------------------------------------

    if inverse==None : #---direct transform---
        
        sig=fun

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
        
        fsig=fun
  
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



# Calculate curvature on R-Z mesh
#
# NOTE: theta component needs to be divided by hthe

def pdiff ( nr, nz, r, z, f):
    print("Calculating DCT...")
    dctf= dct2dslow( f )
    print("Finished DCT")

    drdi = deriv(r)
    dzdi = deriv(z)

    # Get field components
    dfdR = numpy.zeros((nr, nz))
    dfdZ = numpy.zeros((nr, nz))
    for i in range (nr) :
        for j in range (nz) :
            g = local_gradient(dctf, i, j, status=None)
            status=g.status
            dfdr=g.dfdr[0][0]
            dfdz=g.dfdz[0][0]

            # dfd* are derivatives wrt the indices. Need to divide by dr/di etc
            dfdR[i,j] = old_div(dfdr,drdi[i])
            dfdZ[i,j] = old_div(dfdz,dzdi[j])
    
    
    return Bunch(r=dfdR, z=dfdZ, phi=0.0)
 

def pdiff_xy ( nr, nz, r, z, f ):
    # Get field components
    dfdR = numpy.zeros((nr, nz))
    dfdZ = numpy.zeros((nr, nz))
  
    for i in range (nz) :
        dfdR[:,i] = deriv(f[:,i], r)
     
  
    for i in range (nr) :
        dfdZ[i,:] = deriv(f[i,:], z)
     
  
    return Bunch(r=dfdR, z=dfdZ, phi=0.0)
 

def curlcyl ( vecR, vecV, gradVr, gradVphi, gradVz ):
#
# Calculate curl of a axisymmetric vector field V
# in cylindrical coords
#
# Inputs: 
#        vecR - location vector in cylindrical components {r:r,z:z}
#        vecV - vector V in cylindrical components {r:Vr,phi:Vphi,z:Vz} 
#        gradVr - gradient of the r-component,     {dVr/dr,dVr/dz}
#        gradVphi - gradient of the phi-component, {dVphi/dr,dVphi/dz}
#        gradVz - gradient of the z-component,     {dVz/dr,dVz/dz}
#
# Output: curl in cylindrical coordinates
#-------------------------------------------------


    curl=Bunch(r=-gradVphi.z, phi=gradVr.z-gradVz.r, z=old_div(vecV.phi,vecR.r)+gradVphi.r)
#
#
#
    return curl
 

def xprod ( v1, v2, minus=None ):
#
# Calculate cross-product of two vectors
# in cylindrical coordinates
#
# Inputs:
#        v1={r,phi,z}
#        v2={r,phi,z}
#
# Output: v1xv2 {r,phi,z}
#---------------------------------------


    r = v1.phi*v2.z-v1.z*v2.phi
    phi = v1.z*v2.r-v1.r*v2.z
    z = v1.r*v2.phi-v1.phi*v2.r

#
    if minus != None :
        res=Bunch(r=-r,phi=-phi,z=-z) 
    else:
        res=Bunch(r=r,phi=phi,z=z)
  

    return res
 

def dotprod ( v1, v2 ):
#
# Calculate dot-product of two vectors
# in cylindrical coordinates
#
# Inputs:
#        v1={r,phi,z}
#        v2={r,phi,z}
#
# Output: (v1,v2)
#---------------------------------------

    res=v1.r*v2.r + v1.phi*v2.phi + v1.z*v2.z

    return res
 

def rz_curvature( mesh, rixy=None, zixy=None ):
    nr = mesh.nr
    nz = mesh.nz
  
    grad_psi = pdiff_xy(nr, nz, mesh.R, mesh.Z, mesh.psi)
  
    R2D = numpy.zeros((nr, nz))
    Z2D = numpy.zeros((nr, nz))
    for i in range (nr) :
        R2D[i,:] = mesh.R[i]
        Z2D[i,:] = mesh.Z
     

    Br = old_div(grad_psi.Z, R2D)
    Bz = old_div(-grad_psi.R, R2D)

    Bphi = numpy.zeros((nr, nz))
    for i in range (nr) :
        for j in range (nz):
            psinorm = old_div((mesh.psi[i,j] - mesh.simagx), (mesh.sibdry - mesh.simagx))
            if psinorm > 1. :
                fpol = mesh.fpol[numpy.size(mesh.fpol)-1]
            else:
                #fpol = INTERPOL(mesh.fpol, mesh.npsigrid, psinorm, /spline)
                sol=interpolate.UnivariateSpline(mesh.npsigrid, mesh.fpol,s=1)
                fpol=sol(psinorm)

                #fpol = SPLINE(mesh.npsigrid, mesh.fpol, psinorm)
             
            Bphi[i,j] = old_div(fpol, mesh.R[i])
  
    # Total B field
    Bpol = numpy.sqrt(Br**2 + Bz**2)
    B = numpy.sqrt(Bphi**2 + Bpol**2)
  
    # DCT method produces very oscillatory solution
    grad_Br_unit   = pdiff_xy(nr, nz, mesh.R, mesh.Z, old_div(Br,B))
    grad_Bz_unit   = pdiff_xy(nr, nz, mesh.R, mesh.Z, old_div(Bz,B))
    grad_Bphi_unit = pdiff_xy(nr, nz, mesh.R, mesh.Z, old_div(Bphi,B))
  
    vecR=Bunch(r=R2D,z=Z2D)
    vecB_unit=Bunch(r=old_div(Br,B),z=old_div(Bz,B),phi=old_div(Bphi,B))
  
    Bpxy = Bpol
    Rxy = R2D
  
    # Get grad phi
    grad_Phi=Bunch(r=0.0,z=0.0,phi=old_div(1.,Rxy)) #-gradient of the toroidal angle
  
    # Curl of unit b vector
    curlb_unit = curlcyl(vecR, vecB_unit, grad_Br_unit, grad_Bphi_unit, grad_Bz_unit)
  
    # Cross product with b to get curvature vector
    curvec   = xprod(vecB_unit,curlb_unit, minus='minus')
    #-unit b cross curvature vector at cell center
    bxcurvec = xprod(vecB_unit,curvec)
  
    # grad Theta (without factor of 1/hthe)
    grad_Theta = xprod(grad_Phi, grad_psi)
    grad_Theta.r   = old_div(grad_Theta.r, Bpxy)
    grad_Theta.z   = old_div(grad_Theta.z, Bpxy)

    #-calculate bxcurvec dotted with grad_psi, grad_theta, and grad_phi
    bxcv = Bunch(psi=dotprod(bxcurvec,grad_psi),  
                theta=dotprod(bxcurvec,grad_Theta),  
                phi=dotprod(bxcurvec,grad_Phi))
  
    return bxcv
 
