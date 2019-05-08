FUNCTION simple_psi, r, z, r0, z0, mag
  nx = N_ELEMENTS(r)
  ny = N_ELEMENTS(z)
  
  f = DBLARR(nx,ny)
  
  mind = 5*SQRT((r[1]-r[0])^2 + (z[1]-z[0])^2)
  FOR i=0, nx-1 DO BEGIN
    FOR j=0, ny-1 DO BEGIN
      d = SQRT( (R[i]-r0)^2 + (Z[j]-z0)^2 )
      F[i,j] = mag / ( d + mind )
    ENDFOR
  ENDFOR
  RETURN, f
END

PRO gridgen, debug=debug, settings=settings
  ; Generate some artificial data
  nx = 200
  ny = 200

  R = FINDGEN(nx) / DOUBLE(nx-1)
  Z = (FINDGEN(ny) / DOUBLE(ny-1)) - 0.5D

  f = simple_psi(R, Z, 0.5D, 0.0D, -1.0D) $
    + simple_psi(R, Z, 0.3D, -0.55D, -0.39D) $
    + simple_psi(R, Z, 0.8D, -0.55D, -0.39D) $
    + simple_psi(R, Z, 0.4D, 0.55D, -0.4D) 
  
  ;;;;;;;;;;;;;;;; Get input settings ;;;;;;;;;;;;;;;
  
  IF KEYWORD_SET(settings) THEN BEGIN
    ; Read all settings from the input file
    
  ENDIF
  
  REPEAT BEGIN
    psi_inner = get_float("Core normalised psi: ")
  ENDREP UNTIL (psi_inner GE 0.D) AND (psi_inner LT 1.0D)

  REPEAT BEGIN
    psi_outer = get_float("Outer normalised psi: ")
  ENDREP UNTIL (psi_outer GT psi_inner)

  REPEAT BEGIN
    nrad = get_integer("Enter number of radial points  : ")
  ENDREP UNTIL nrad GT 1

  REPEAT BEGIN
    npol = get_integer("Enter number of poloidal points: ")
  ENDREP UNTIL npol GT 1

  ; Get separatrix peaking factors
  rad_peaking = get_float("Enter radial peaking factor  : ")
  pol_peaking = get_float("Enter poloidal peaking factor: ")
  
  ;;;;;;;;;;;;;;;; Call grid generator ;;;;;;;;;;;;;;
  
  settings = {psi_inner:psi_inner, psi_outer:psi_outer, $  ; Normalised psi range
              nrad:nrad, npol:npol, $    ; Number of radial and poloidal points
              rad_peaking:rad_peaking, pol_peaking:pol_peaking}

  grid = create_grid(F, R, Z, settings, /simple)

END
