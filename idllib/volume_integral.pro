; Integrate over a volume

FUNCTION volume_integral, var, grid, xrange=xrange
  s = SIZE(var)

  IF s[0] EQ 4 THEN BEGIN
    ; 4D [x,y,z,t] - integrate for each t
    nx = s[1]
    ny = s[2]
    nt = s[4]
    
    result = FLTARR(nt)
    FOR t=0,nt-1 DO result[t] = volume_integral(var[*,*,*,t], grid, xrange=xrange)
    RETURN, result
  ENDIF ELSE IF s[0] EQ 3 THEN BEGIN
    ; 3D [x,y,z] - average in Z
    nx = s[1]
    ny = s[2]
    nz = s[3]
    
    zi = FLTARR(nx, ny)
    FOR x=0, nx-1 DO BEGIN
      FOR y=0, ny-1 DO BEGIN
        zi[x,y] = MEAN(var[x,y,*])
      ENDFOR
    ENDFOR
    RETURN, volume_integral(zi, grid, xrange=xrange)
  ENDIF ELSE IF s[0] NE 2 THEN BEGIN
    PRINT, "ERROR: volume_integral var must be 2, 3 or 4D"
  ENDIF
  
  ; 2D [x,y]
  nx = s[1]
  ny = s[2]
  
  IF NOT KEYWORD_SET(xrange) THEN xrange=[0,nx-1]
  
  result = 0.0d

  ;status = gen_surface(mesh=grid) ; Start generator
  xi = -1
  yi = INDGEN(ny)
  last = 0
  iy = FLTARR(nx)
  REPEAT BEGIN
    ;yi = gen_surface(last=last, xi=xi, period=periodic)
    xi = xi + 1
    IF xi EQ nx-1 THEN last = 1
    
    IF (xi GE MIN(xrange)) AND (xi LE MAX(xrange)) THEN BEGIN
      dtheta = 2.*!PI / FLOAT(ny)
      r = REFORM(grid.Rxy[xi,yi])
      z = REFORM(grid.Zxy[xi,yi])
      n = N_ELEMENTS(r)
      dl = SQRT( DERIV(r)^2 + DERIV(z)^2 ) / dtheta
      
      ; Area of flux-surface
      dA = (REFORM(grid.Bxy[xi,yi]/grid.Bpxy[xi,yi])*dl) * (r*2.*!PI)
      ; Volume
      IF xi EQ nx-1 THEN dpsi = REFORM(grid.psixy[xi,yi] - grid.psixy[xi-1,yi]) ELSE dpsi = REFORM(grid.psixy[xi+1,yi] - grid.psixy[xi,yi])
      
      dV = dA * dpsi / (r*REFORM(grid.Bpxy[xi,yi])) ; May need factor of 2pi
      dV = ABS(dV)
      
      result = result + TOTAL(var[xi,yi] * dV)
    ENDIF
  ENDREP UNTIL last
  
  RETURN, result
END
