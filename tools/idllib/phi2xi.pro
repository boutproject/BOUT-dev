; Calculate the radial (psi) component of the MHD displacement
; vector from the electrostatic potential phi
; Actually returns the psi component of the ExB velocity
;
; NOTE: Does not take into account branch-cuts
;
; Input
;   u.dy, u.hthe
;   u.Btxy, u.Bpxy, u.Bxy
;   u.Rxy
;
; Changelog
; 16 Oct 2008: First version, B.Dudson

; Take derivative in y, taking into account branch-cuts
FUNCTION ddy, var, mesh, dy=dy
  f = var

  IF NOT KEYWORD_SET(dy) THEN dy = 2.*!PI / FLOAT(TOTAL(mesh.npol))

  status = gen_surface(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=period)
    IF period THEN BEGIN
       f[xi,yi] = fft_deriv(var[xi,yi])
    ENDIF ELSE f[xi,yi] = DERIV(var[xi,yi])
  ENDREP UNTIL last
  RETURN, f / dy
END

FUNCTION phi2xi_3d, phi, u, period=period
  ; this is the main function. Takes 3D (x,y,z) variable
  ON_ERROR, 2

  s = SIZE(phi, /dim)

  IF N_ELEMENTS(s) NE 3 THEN BEGIN
    PRINT, "Error: phi must be 3D"
    RETURN, 0
  ENDIF

  nx = s[0]
  ny = s[1]
  nz = s[2]

  IF NOT is_pow2(nz) THEN BEGIN
    nz = nz - 1  ; data must include the "extra" data point
    IF NOT is_pow2(nz) THEN BEGIN
      PRINT, "Error: All Z points expected"
    ENDIF
  ENDIF

  IF NOT KEYWORD_SET(period) THEN period = 10

  dz = 2.0 * !PI / FLOAT(period * nz)
  
  ; Take derivatives

  dphidy = FLTARR(nx,ny,nz)
  dphidz = dphidy

  ; d/dy

  FOR z=0, nz-1 DO BEGIN
    dphidy[*,*,z] = DDY( REFORM(phi[*,*,z]), u, dy=u.dy )
  ENDFOR

  ; d/dz
  
  FOR x=0, nx-1 DO BEGIN
    FOR y=0, ny-1 DO BEGIN
      dphidz[x,y,*] = fft_deriv(phi[x,y,*])/dz
    ENDFOR 
  ENDFOR

  ; calculate result

  result = FLTARR(nx, ny, nz)

  FOR z=0, nz-1 DO BEGIN
    result[*,*,z] = (u.Btxy * u.Rxy * dphidy[*,*,z] $
                     - u.Bxy*u.Bxy*u.hthe*dphidz[*,*,z] / u.Bpxy) / (u.Bxy*u.Bxy * u.hthe * u.Rxy)
  ENDFOR

  RETURN, result
END

FUNCTION phi2xi, phi, u, period=period
  ON_ERROR, 2
  
  s = SIZE(phi, /dim)
  
  IF N_ELEMENTS(s) EQ 3 THEN BEGIN
    RETURN, phi2xi_3d(phi, u)
  ENDIF ELSE IF N_ELEMENTS(s) EQ 4 THEN BEGIN
    result = phi

    FOR t=0, s[3]-1 DO BEGIN
      result[*,*,*,t] = phi2xi_3d(REFORM(phi[*,*,*,t]), u)
    ENDFOR
    
    RETURN, result
  ENDIF ELSE BEGIN
    PRINT, "Error in phi2xi: Variable must be 3 or 4D"
  ENDELSE
  RETURN, 0
END
