; Calculates curvature from GATO grid
; adapted from M.Umansky's code

FUNCTION pdiff_rz, rxy, zxy, fxy, i, j
   s = SIZE(rxy, /dim)
   nx = s[0]
   ny = s[1]

   jp = (j+1) MOD ny
   jm = (j-1 + ny) MOD ny

   r = [rxy[i-1,jm], rxy[i+1,jm], rxy[i+1,jp], rxy[i-1,jp]]
   z = [zxy[i-1,jm], zxy[i+1,jm], zxy[i+1,jp], zxy[i-1,jp]]
   f = [fxy[i-1,jm], fxy[i+1,jm], fxy[i+1,jp], fxy[i-1,jp]]
   


   ;IF j EQ 0 THEN STOP

   A=TRANSPOSE([[fltarr(4)+1],[r-r(0)],[z-z(0)]])

   SVDC, A,W,U,V

   res=SVSOL(U,W,V,f)

   pdiff={r:res[1],z:res[2],phi:0.0}

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

PRO curvature, nx, ny, Rxy, Zxy, BRxy, BZxy, BPHIxy, PSIxy, THETAxy,$
               CURLB=CURLB, JXB=JXB, CURVEC=CURVEC, BXCURVEC=BXCURVEC, BXCV=BXCV,$
               DEBUG=DEBUG
;
; Calculate the magnetic field curvature and other related quantities
;--------------------------------------------------------------------

   PRINT, 'Calculating curvature-related quantities...'

;;-vector quantities are stored as 2D arrays of structures {r,phi,z}
   vec={r:0.,phi:0.,z:0.}
   curlb=REPLICATE(vec,nx,ny) 
   jxb=REPLICATE(vec,nx,ny) 
   curvec=REPLICATE(vec,nx,ny) 
   bxcurvec=REPLICATE(vec,nx,ny)

   vec2={psi:0.,theta:0.,phi:0.}
   bxcv=REPLICATE(vec2,nx,ny)

   FOR i=1,nx-2 DO BEGIN
       FOR j=0,ny-1 DO BEGIN
           x = i
           y = j
           
           grad_Br   = pdiff_rz(Rxy, Zxy, BRxy, x, y)
           grad_Bz   = pdiff_rz(Rxy, Zxy, BZxy, x, y)
           grad_Bphi = pdiff_rz(Rxy, Zxy, BPHIxy, x, y)

           grad_Psi  = pdiff_rz(Rxy, Zxy, PSIxy, x, y)
           IF (j EQ 0) OR (j EQ ny-1) THEN BEGIN
               grad_Theta= pdiff_rz(Rxy, Zxy, transpose([transpose(THETAxy[*,ny/2:*]), transpose(THETAxy[*,0:ny/2-1])]), x, y)
           ENDIF ELSE grad_Theta= pdiff_rz(Rxy, Zxy, THETAxy, x, y)
           
           grad_Phi={r:0.0,z:0.0,phi:1./Rxy[x,y]} ;-gradient of the toroidal angle

           vecR={r:Rxy[x,y],z:Zxy[x,y]}
           vecB={r:BRxy[x,y],z:BZxy[x,y],phi:BPHIxy[x,y]}

           curlb[i,j]=CurlCyl(vecR, vecB, grad_Br, grad_Bphi, grad_Bz)
           jxb[i,j]=Xprod(curlb[i,j], vecB)

           ;-magnitude of B at 5 locations in cell
           bstrength = SQRT(BRxy^2 + BZxy^2 + BPHIxy^2)

           ;-unit B vector at cell center
           vecB_unit={r:BRxy[x,y]/bStrength[x,y], $
                      z:BZxy[x,y]/bStrength[x,y], $
                      phi:BPHIxy[x,y]/bStrength[x,y]}

           ;-components of gradient of unit B vector at 5 locations in cell
           grad_Br_unit = pdiff_rz(Rxy, Zxy, BRxy/bStrength, x, y)

           grad_Bz_unit = pdiff_rz(Rxy, Zxy, BZxy/bStrength, x, y)
           
           grad_Bphi_unit = pdiff_rz(Rxy, Zxy, BPHIxy/bStrength, x, y)

           ;-curl of unit B vector at cell center
           curlb_unit=CurlCyl(vecR, vecB_unit, grad_Br_unit, grad_Bphi_unit, grad_Bz_unit)

           ;-curvature vector at cell center
           curvec[i,j]=Xprod(vecB_unit,curlb_unit,/MINUS)

           ;-unit b cross curvature vector at cell center
           bxcurvec[i,j]=Xprod(vecB_unit,curvec[i,j])


           ;-calculate bxcurvec dotted with grad_psi, grad_theta, and grad_phi
           bxcv[i,j].psi=Dotprod(bxcurvec[i,j],grad_Psi)
           bxcv[i,j].theta=Dotprod(bxcurvec[i,j],grad_Theta)
           bxcv[i,j].phi=Dotprod(bxcurvec[i,j],grad_Phi)

       ENDFOR
   ENDFOR
   
   IF KEYWORD_SET(DEBUG) THEN STOP
   PRINT, '...done'

END

