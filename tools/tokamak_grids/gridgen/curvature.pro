; Calculates curvature from GATO grid
; adapted from M.Umansky's code

FUNCTION pdiff_rz, rxy, zxy, fxy, i, j, jp, jm
   s = SIZE(rxy, /dim)
   nx = s[0]
   ny = s[1]

   im = (i-1) > 0
   ip = (i+1) < (nx-1)

   r = [rxy[im,jm], rxy[ip,jm], rxy[ip,jp], rxy[im,jp]]
   z = [zxy[im,jm], zxy[ip,jm], zxy[ip,jp], zxy[im,jp]]
   f = [fxy[im,jm], fxy[ip,jm], fxy[ip,jp], fxy[im,jp]]
   
   ;IF j EQ 0 THEN STOP

   A=TRANSPOSE([[DBLARR(4)+1],[r-r(0)],[z-z(0)]])

   SVDC, A,W,U,V

   res=SVSOL(U,W,V,f)

   pdiff={r:res[1],z:res[2],phi:0.0D}

   RETURN, pdiff
END

function curlcyl, vecR, vecV, gradVr, gradVphi, gradVz
;
; Calculate curl of a axisymmetric vector field V
; in cylindrical coords
;
; Inputs: 
;        vecR - location vector in cylindrical components {r:r,z:z}
;        vecV - vector V in cylindrical components {r:Vr,phi:Vphi,z:Vz} 
;        gradVr - gradient of the r-component,     {dVr/dr,dVr/dz}
;        gradVphi - gradient of the phi-component, {dVphi/dr,dVphi/dz}
;        gradVz - gradient of the z-component,     {dVphi/dr,dVphi/dz}
;
; Output: curl in cylindrical coordinates
;-------------------------------------------------


  curl={r:-gradVphi.z, phi:gradVr.z-gradVz.r, z:vecV.phi/vecR.r+gradVphi.r}
;
;
;
return, curl
end

function xprod, v1, v2, minus=minus
;
; Calculate cross-product of two vectors
; in cylindrical coordinates
;
; Inputs:
;        v1={r,phi,z}
;        v2={r,phi,z}
;
; Output: v1xv2 {r,phi,z}
;---------------------------------------


    r = v1.phi*v2.z-v1.z*v2.phi
    phi = v1.z*v2.r-v1.r*v2.z
    z = v1.r*v2.phi-v1.phi*v2.r

;
 if keyword_set(MINUS) then begin
   res={r:-r,phi:-phi,z:-z} 
 endif else begin
   res={r:r,phi:phi,z:z}
 endelse

return, res
end

function dotprod, v1, v2
;
; Calculate dot-product of two vectors
; in cylindrical coordinates
;
; Inputs:
;        v1={r,phi,z}
;        v2={r,phi,z}
;
; Output: (v1,v2)
;---------------------------------------

    res=v1.r*v2.r + v1.phi*v2.phi + v1.z*v2.z

return, res
end

PRO curvature, nx, ny, Rxy, Zxy, BRxy, BZxy, BPHIxy, PSIxy, THETAxy, HTHExy, $
               CURLB=CURLB, JXB=JXB, CURVEC=CURVEC, BXCURVEC=BXCURVEC, BXCV=BXCV,$
               DEBUG=DEBUG, mesh=mesh
;
; Calculate the magnetic field curvature and other related quantities
;--------------------------------------------------------------------

   PRINT, 'Calculating curvature-related quantities...'
   
;;-vector quantities are stored as 2D arrays of structures {r,phi,z}
   vec={r:0.D,phi:0.D,z:0.D}
   curlb=REPLICATE(vec,nx,ny) 
   jxb=REPLICATE(vec,nx,ny) 
   curvec=REPLICATE(vec,nx,ny) 
   bxcurvec=REPLICATE(vec,nx,ny)

   vec2={psi:0.D,theta:0.D,phi:0.D}
   bxcv=REPLICATE(vec2,nx,ny)

   status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
   REPEAT BEGIN
     yi = gen_surface_hypnotoad(last=last, xi=x, period=period)
     nys = N_ELEMENTS(yi)

     ; Get vector along the surface
     IF period THEN BEGIN
       dr = fft_deriv(Rxy[x,yi])
       dz = fft_deriv(Zxy[x,yi])
     ENDIF ELSE BEGIN
       dr = DERIV(Rxy[x,yi])
       dz = DERIV(Zxy[x,yi])
     ENDELSE
     dl = SQRT(dr^2 + dz^2)
     dr = dr / dl
     dz = dz / dl
     
     FOR j=0, nys-1 DO BEGIN
       y = yi[j]
       
       IF period THEN BEGIN
         yp = yi[ (j+1)     MOD nys ]
         ym = yi[ (j-1+nys) MOD nys ]
       ENDIF ELSE BEGIN
         yp = yi[ (j+1) < (nys-1) ]
         ym = yi[ (j-1) > 0 ]
       ENDELSE
       
       grad_Br   = pdiff_rz(Rxy, Zxy, BRxy, x, y, yp, ym)
       grad_Bz   = pdiff_rz(Rxy, Zxy, BZxy, x, y, yp, ym)
       grad_Bphi = pdiff_rz(Rxy, Zxy, BPHIxy, x, y, yp, ym)
       
       grad_Psi  = pdiff_rz(Rxy, Zxy, PSIxy, x, y, yp, ym)
       
       ;grad_Theta = pdiff_rz(Rxy, Zxy, THETAxy, x, y, yp, ym)
       grad_Theta = {r:dr[j]/hthexy[x,y], z:dz[j]/hthexy[x,y], phi:0.0D}

       grad_Phi={r:0.0D,z:0.0D,phi:1.D/Rxy[x,y]} ;-gradient of the toroidal angle

       vecR={r:Rxy[x,y],z:Zxy[x,y]}
       vecB={r:BRxy[x,y],z:BZxy[x,y],phi:BPHIxy[x,y]}
       
       curlb[x,y]=CurlCyl(vecR, vecB, grad_Br, grad_Bphi, grad_Bz)
       jxb[x,y]=Xprod(curlb[x,y], vecB)
       
       ;-magnitude of B at 5 locations in cell
       bstrength = SQRT(BRxy^2 + BZxy^2 + BPHIxy^2)
       
       ;-unit B vector at cell center
       vecB_unit={r:BRxy[x,y]/bStrength[x,y], $
                  z:BZxy[x,y]/bStrength[x,y], $
                  phi:BPHIxy[x,y]/bStrength[x,y]}
       
       ;-components of gradient of unit B vector at 5 locations in cell
       grad_Br_unit = pdiff_rz(Rxy, Zxy, BRxy/bStrength, x, y, yp, ym)
       
       grad_Bz_unit = pdiff_rz(Rxy, Zxy, BZxy/bStrength, x, y, yp, ym)
       
       grad_Bphi_unit = pdiff_rz(Rxy, Zxy, BPHIxy/bStrength, x, y, yp, ym)

       ;-curl of unit B vector at cell center
       curlb_unit=CurlCyl(vecR, vecB_unit, grad_Br_unit, grad_Bphi_unit, grad_Bz_unit)

       ;-curvature vector at cell center
       curvec[x,y]=Xprod(vecB_unit,curlb_unit,/MINUS)

       ;-unit b cross curvature vector at cell center
       bxcurvec[x,y]=Xprod(vecB_unit,curvec[x,y])

       
       ;-calculate bxcurvec dotted with grad_psi, grad_theta, and grad_phi
       bxcv[x,y].psi=Dotprod(bxcurvec[x,y],grad_Psi)
       bxcv[x,y].theta=Dotprod(bxcurvec[x,y],grad_Theta)
       bxcv[x,y].phi=Dotprod(bxcurvec[x,y],grad_Phi)
       
     ENDFOR
   ENDREP UNTIL last
   
   IF KEYWORD_SET(DEBUG) THEN STOP
   PRINT, '...done'

END

