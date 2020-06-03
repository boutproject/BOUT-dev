from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
#; Calculates curvature from GATO grid
#; adapted from M.Umansky's code
import numpy
from boututils.bunch import Bunch
from gen_surface import gen_surface
from boututils.calculus import deriv
from boututils.fft_deriv import fft_deriv
import sys


def pdiff_rz( rxy, zxy, fxy, i, j, jp, jm):
   s = numpy.shape(rxy)
   nx = s[0]
   ny = s[1]

   im = numpy.max([i-1 , 0])
   ip = numpy.min([i+1 , nx-1])


   r = [rxy[im,jm], rxy[ip,jm], rxy[ip,jp], rxy[im,jp]]
   z = [zxy[im,jm], zxy[ip,jm], zxy[ip,jp], zxy[im,jp]]
   f = [fxy[im,jm], fxy[ip,jm], fxy[ip,jp], fxy[im,jp]]

   #IF j EQ 0 THEN STOP

   A=numpy.transpose([numpy.zeros(4)+1,r-r[0],z-z[0]])

   res, _, _, _ = numpy.linalg.lstsq(A,f)


   pdiff=Bunch(r=res[1],z=res[2],phi=0.0)

   return pdiff


def curlcyl( vecR, vecV, gradVr, gradVphi, gradVz):
#;
#; Calculate curl of a axisymmetric vector field V
#; in cylindrical coords
#;
#; Inputs:
#;        vecR - location vector in cylindrical components {r:r,z:z}
#;        vecV - vector V in cylindrical components {r:Vr,phi:Vphi,z:Vz}
#;        gradVr - gradient of the r-component,     {dVr/dr,dVr/dz}
#;        gradVphi - gradient of the phi-component, {dVphi/dr,dVphi/dz}
#;        gradVz - gradient of the z-component,     {dVphi/dr,dVphi/dz}
#;
#; Output: curl in cylindrical coordinates
#;-------------------------------------------------


    curl=Bunch(r=-gradVphi.z, phi=gradVr.z-gradVz.r, z=old_div(vecV.phi,vecR.r)+gradVphi.r)

    return curl


def xprod( v1, v2, minus=None):
#;
#; Calculate cross-product of two vectors
#; in cylindrical coordinates
#;
#; Inputs:
#;        v1={r,phi,z}
#;        v2={r,phi,z}
#;
#; Output: v1xv2 {r,phi,z}
#;---------------------------------------


    r = v1.phi*v2.z-v1.z*v2.phi
    phi = v1.z*v2.r-v1.r*v2.z
    z = v1.r*v2.phi-v1.phi*v2.r


    if minus :
        res=Bunch(r=-r,phi=-phi,z=-z)
    else:
        res=Bunch(r=r,phi=phi,z=z)


    return res


def dotprod( v1, v2):
#;
#; Calculate dot-product of two vectors
#; in cylindrical coordinates
#;
#; Inputs:
#;        v1={r,phi,z}
#;        v2={r,phi,z}
#;
#; Output: (v1,v2)
#;---------------------------------------

    res=v1.r*v2.r + v1.phi*v2.phi + v1.z*v2.z

    return res


def curvature( nx, ny, Rxy, Zxy, BRxy, BZxy, BPHIxy, PSIxy, THETAxy, hthexy,
               CURLB=None, JXB=None, CURVEC=None, BXCURVEC=None, BXCV=None,
               DEBUG=None, mesh=None):
#;
#; Calculate the magnetic field curvature and other related quantities
#;--------------------------------------------------------------------

    print('Calculating curvature-related quantities...')

#;;-vector quantities are stored as 2D arrays of structures {r,phi,z}
    vec=Bunch( r=0.,phi=0.,z=0.)
    curlb=numpy.tile(vec,(nx,ny))
    jxb=numpy.tile(vec,(nx,ny))
    curvec=numpy.tile(vec,(nx,ny))
    bxcurvec=numpy.tile(vec,(nx,ny))

    bxcv=Bunch()
    bxcv.psi=numpy.zeros((nx,ny))
    bxcv.theta=numpy.zeros((nx,ny))
    bxcv.phi=numpy.zeros((nx,ny))


    status = gen_surface(mesh=mesh) # Start generator


    while True:
        period, yi, xi, last = gen_surface(period=None, last=None, xi=None)
        nys = numpy.size(yi)
        x=xi


     # Get vector along the surface
        if period ==1 :
            dr = fft_deriv(Rxy[x,yi])
            dz = fft_deriv(Zxy[x,yi])
        else:
            dr = deriv(Rxy[x,yi])
            dz = deriv(Zxy[x,yi])

        dl = numpy.sqrt(dr**2 + dz**2)

        dr = old_div(dr, dl)
        dz = old_div(dz, dl)



        for j in range (nys) :
            y = yi[j]

            if period :
                yp = yi[ (j+1)     % nys ]
                ym = yi[ (j-1+nys) % nys ]
            else:
                yp = yi[ numpy.min([j+1 , nys-1]) ]
                ym = yi[ numpy.max([j-1 , 0]) ]


            grad_Br   = pdiff_rz(Rxy, Zxy, BRxy, x, y, yp, ym)
            grad_Bz   = pdiff_rz(Rxy, Zxy, BZxy, x, y, yp, ym)
            grad_Bphi = pdiff_rz(Rxy, Zxy, BPHIxy, x, y, yp, ym)



            grad_Psi  = pdiff_rz(Rxy, Zxy, PSIxy, x, y, yp, ym)


            #grad_Theta = pdiff_rz(Rxy, Zxy, THETAxy, x, y, yp, ym)
            grad_Theta = Bunch( r=old_div(dr[j],hthexy[x,y]), z=old_div(dz[j],hthexy[x,y]), phi=0.0 )

            grad_Phi=Bunch( r=0.0,z=0.0,phi=old_div(1.,Rxy[x,y]) ) #-gradient of the toroidal angle

            vecR=Bunch(  r=Rxy[x,y],z=Zxy[x,y] )
            vecB=Bunch( r=BRxy[x,y],z=BZxy[x,y],phi=BPHIxy[x,y] )


            curlb[x,y]=curlcyl(vecR, vecB, grad_Br, grad_Bphi, grad_Bz)


            jxb[x,y]=xprod(curlb[x,y], vecB)


            #-magnitude of B at 5 locations in cell
            bstrength = numpy.sqrt(BRxy**2 + BZxy**2 + BPHIxy**2)

            #-unit B vector at cell center
            vecB_unit=Bunch( r=old_div(BRxy[x,y],bstrength[x,y]),
                  z=old_div(BZxy[x,y],bstrength[x,y]),
                  phi=old_div(BPHIxy[x,y],bstrength[x,y]) )

            #-components of gradient of unit B vector at 5 locations in cell
            grad_Br_unit = pdiff_rz(Rxy, Zxy, old_div(BRxy,bstrength), x, y, yp, ym)

            grad_Bz_unit = pdiff_rz(Rxy, Zxy, old_div(BZxy,bstrength), x, y, yp, ym)

            grad_Bphi_unit = pdiff_rz(Rxy, Zxy, old_div(BPHIxy,bstrength), x, y, yp, ym)

            #-curl of unit B vector at cell center
            curlb_unit=curlcyl(vecR, vecB_unit, grad_Br_unit, grad_Bphi_unit, grad_Bz_unit)

            #-curvature vector at cell center
            curvec[x,y]=xprod(vecB_unit,curlb_unit,minus='MINUS')

            #-unit b cross curvature vector at cell center
            bxcurvec[x,y]=xprod(vecB_unit,curvec[x,y])

            #-calculate bxcurvec dotted with grad_psi, grad_theta, and grad_phi
            bxcv.psi[x,y]=dotprod(bxcurvec[x,y],grad_Psi)
            bxcv.theta[x,y]=numpy.real(dotprod(bxcurvec[x,y],grad_Theta))
            bxcv.phi[x,y]=dotprod(bxcurvec[x,y],grad_Phi)


        if last==1 : break

#   if DEBUG : sys.exit()

    print('...done')

    return bxcv
