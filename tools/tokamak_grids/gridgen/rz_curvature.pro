; Calculate curvature on R-Z mesh
;
; NOTE: theta component needs to be divided by hthe

FUNCTION pdiff, nr, nz, r, z, f
  PRINT, "Calculating DCT..."
  DCT2Dslow, f, dctf
  PRINT, "Finished DCT"

  drdi = DERIV(R)
  dzdi = DERIV(Z)

  ; Get field components
  dfdR = DBLARR(nr, nz)
  dfdZ = dfdR
  FOR i=0,nr-1 DO BEGIN
    FOR j=0,nz-1 DO BEGIN
      g = local_gradient(dctf, i, j, status=status)
      ; dfd* are derivatives wrt the indices. Need to divide by dr/di etc
      dfdR[i,j] = g.dfdr/drdi[i]
      dfdZ[i,j] = g.dfdz/dzdi[j]
    ENDFOR
  ENDFOR
  
  RETURN, {r:dfdR, z:dfdZ, phi:0.0D}
END

FUNCTION pdiff_xy, nr, nz, r, z, f
  ; Get field components
  dfdR = DBLARR(nr, nz)
  dfdZ = dfdR
  
  FOR i=0,nz-1 DO BEGIN
    dfdR[*,i] = DERIV(r, f[*,i])
  ENDFOR
  FOR i=0, nr-1 DO BEGIN
    dfdZ[i,*] = DERIV(z, f[i,*])
  ENDFOR
  
  RETURN, {r:dfdR, z:dfdZ, phi:0.0D}
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
;        gradVz - gradient of the z-component,     {dVz/dr,dVz/dz}
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

FUNCTION rz_curvature, mesh, rixy=rixy, zixy=zixy
  nr = mesh.nr
  nz = mesh.nz
  
  grad_Psi = pdiff_xy(nr, nz, mesh.R, mesh.Z, mesh.psi)
  
  R2D = DBLARR(nr, nz)
  Z2D = DBLARR(nr, nz)
  FOR i=0,nr-1 DO BEGIN
    R2D[i,*] = mesh.R[i]
    Z2D[i,*] = mesh.Z
  ENDFOR

  Br = grad_psi.Z / R2D
  Bz = -grad_psi.R / R2D

  Bphi = DBLARR(nr, nz)
  FOR i=0,nr-1 DO BEGIN
    FOR j=0,nz-1 DO BEGIN
      psinorm = (mesh.psi[i,j] - mesh.simagx) / (mesh.sibdry - mesh.simagx)
      IF psinorm GT 1.D THEN BEGIN
        fpol = mesh.fpol[N_ELEMENTS(mesh.fpol)-1]
      ENDIF ELSE BEGIN
        ;fpol = INTERPOL(mesh.fpol, mesh.npsigrid, psinorm, /spline)
        fpol = SPLINE(mesh.npsigrid, mesh.fpol, psinorm)
      ENDELSE
      Bphi[i,j] = fpol / mesh.R[i]
    ENDFOR
  ENDFOR
  
  ; Total B field
  Bpol = SQRT(Br^2 + Bz^2)
  B = SQRT(Bphi^2 + Bpol^2)
  
  ; DCT method produces very oscillatory solution
  grad_Br_unit   = pdiff_xy(nr, nz, mesh.R, mesh.Z, Br/B)
  grad_Bz_unit   = pdiff_xy(nr, nz, mesh.R, mesh.Z, Bz/B)
  grad_Bphi_unit = pdiff_xy(nr, nz, mesh.R, mesh.Z, Bphi/B)
  
  vecR={r:R2D,z:Z2D}
  vecB_unit={r:Br/B,z:Bz/B,phi:Bphi/B}
  
  Bpxy = Bpol
  Rxy = R2D
  
  ; Get grad phi
  grad_Phi={r:0.0D,z:0.0D,phi:1.D/Rxy} ;-gradient of the toroidal angle
  
  ; Curl of unit b vector
  curlb_unit = CurlCyl(vecR, vecB_unit, grad_Br_unit, grad_Bphi_unit, grad_Bz_unit)
  
  ; Cross product with b to get curvature vector
  curvec   = Xprod(vecB_unit,curlb_unit, /MINUS)
  ;-unit b cross curvature vector at cell center
  bxcurvec = Xprod(vecB_unit,curvec)
  
  ; grad Theta (without factor of 1/hthe)
  grad_Theta = Xprod(grad_Phi, grad_Psi)
  grad_Theta.r   = grad_Theta.r   / Bpxy
  grad_Theta.z   = grad_Theta.z   / Bpxy

  ;-calculate bxcurvec dotted with grad_psi, grad_theta, and grad_phi
  bxcv = {psi:Dotprod(bxcurvec,grad_Psi), $
          theta:Dotprod(bxcurvec,grad_Theta), $
          phi:Dotprod(bxcurvec,grad_Phi)}
  
  RETURN, bxcv
END
