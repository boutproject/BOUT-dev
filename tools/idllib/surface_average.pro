; Perform surface average
;
; var - [x,y,z] or [x,y,z,t]
;
; 
; area=area : Average by flux-surface area = (B/Bp)*dl * R*dz
;
; By default, averages over poloidal angle such that
; surface_average(nu) = q
;

FUNCTION surface_average, var, grid, area=area, simple=simple
  
  s = SIZE(var)
  
  IF s[0] EQ 4 THEN BEGIN
    nx = s[1]
    ny = s[2]
    nt = s[4]
    
    result = FLTARR(nx,nt)
    FOR t=0,nt-1 DO result[*,t] = surface_average(var[*,*,*,t])
    
  ENDIF ELSE IF s[0] NE 3 THEN BEGIN
    PRINT, "ERROR: surface_average var must be 3 or 4D"
    RETURN, 0
  ENDIF
  
  ; 3D [x,y,z]
  nx = s[1]
  ny = s[2]
  nz = s[3]
  
  ; Calculate poloidal angle from grid
  theta = FLTARR(nx,ny)
  
  ;status = gen_surface(mesh=grid) ; Start generator
  xi = -1
  yi = INDGEN(ny)
  last = 0
  REPEAT BEGIN
    ;yi = gen_surface(last=last, xi=xi, period=periodic)
    xi = xi + 1
    IF xi EQ nx-1 THEN last = 1
    
    dtheta = 2.*!PI / FLOAT(ny)
    r = REFORM(grid.Rxy[xi,yi])
    z = REFORM(grid.Zxy[xi,yi])
    n = N_ELEMENTS(r)
    dl = SQRT( DERIV(r)^2 + DERIV(z)^2 ) / dtheta
    IF KEYWORD_SET(area) THEN BEGIN
      dA = REFORM(grid.Bxy[xi,yi]/grid.Bpxy[xi,yi])*r*dl
      A = int_func(FINDGEN(n),dA,/simple=simple)
      theta[xi,yi] = 2.*!PI*A/A[n-1]
    ENDIF ELSE BEGIN
      nu = dl * REFORM(grid.Btxy[xi,yi]) / ( REFORM(grid.Bpxy[xi,yi]) * r )
      theta[xi,yi] = int_func(FINDGEN(n)*dtheta,nu,/simple=simple)
      theta[xi,yi] = 2.*!PI*theta[xi,yi] / theta[xi,yi[n-1]]
    ENDELSE
  ENDREP UNTIL last
  
  vy = FLTARR(ny)
  result = FLTARR(nx)
  IF KEYWORD_SET(simple) THEN BEGIN
    FOR x=0,nx-1 DO BEGIN
      FOR y=0,ny-1 DO vy[y] = MEAN(var[x,y,*])
      result[x] = INT_TRAPEZOID(REFORM(theta[x,*]), vy) / (2.*!PI)
    ENDFOR
  ENDIF ELSE BEGIN
    FOR x=0,nx-1 DO BEGIN
      FOR y=0,ny-1 DO vy[y] = MEAN(var[x,y,*])
      result[x] = INT_TABULATED(REFORM(theta[x,*]), vy) / (2.*!PI)
    ENDFOR
  ENDELSE
  
  RETURN, result
END
