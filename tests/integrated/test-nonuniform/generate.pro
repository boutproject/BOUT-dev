; Generate a grid file for the Delp2 test
; Since each y location is inverted separately, they can all have
; different mesh spacings
; 
; For each mesh spacing, the same set of data are used
; 
; Usage:
; =====
;
; IDL> generate, nx=nx
;

; Generate data for test cases
FUNCTION gen_data, xpos, nx, nz, ref=ref
  in = DBLARR(nx, 4, nz)
  ref = in

  in[*,0,0]  = xpos
  ref[*,0,0] = 0.d

  in[*,1,0]  = xpos^2
  ref[*,1,0] = 2.d

  in[*,2,0]  = COS(4.25*!PI*xpos)
  ref[*,2,0] = -(4.25*!PI)^2 * COS(4.25*!PI*xpos)

  in[*,3,0]  = ALOG(COSH(xpos*8. - 4.))
  ref[*,3,0] =64./( COSH(xpos*8 - 4)^2 )

  FOR z=1,nz-1 DO BEGIN
    in[*,0:3,z] = in[*,0:3,0]
    ref[*,0:3,z] = ref[*,0:3,0]
  ENDFOR

  RETURN, in
END

FUNCTION add_data, xpos, dx, nz, a=a
  nx = N_ELEMENTS(xpos)
  
  in = gen_data(xpos, nx, nz, ref=ref)
  ny = (SIZE(in, /dim))[2] ; Number of y values returned
  
  xpos2 = DBLARR(nx, ny)
  dx2   = xpos2
  FOR y=0, ny-1 DO BEGIN
    xpos2[*,y] = xpos
    dx2[*,y] = dx
  ENDFOR

  IF NOT KEYWORD_SET(a) THEN BEGIN
    d = {nx:nx, ny:ny, nz:nz, $
         xpos:xpos2, dx:dx2, $
         input:in, reference:ref}
  ENDIF ELSE BEGIN
    ny2 = a.ny+ny
    in2 = DBLARR(nx, ny2, nz)
    in2[*,0:(a.ny-1), *] = a.input
    in2[*,a.ny:*, *] = in
    
    ref2 = DBLARR(nx, ny2, nz)
    ref2[*,0:(a.ny-1), *] = a.reference
    ref2[*,a.ny:*,*] = ref
    
    d = {nx:nx, ny:ny2, nz:nz, $
         xpos:[[a.xpos], [xpos2]], dx:[[a.dx], [dx2]], $
         input:in2, reference:ref2 }
  ENDELSE

  RETURN, d
END

PRO generate, nx=nx
  IF NOT KEYWORD_SET(nx) THEN nx = 20
  
  nz = 4 
  
  ; First uniform mesh tests
  xpos = DINDGEN(nx)/DOUBLE(nx-1)
  dx   = 1.d / DOUBLE(nx-1)
  d = add_data(xpos, dx, nz)
  
  ; Constant change in mesh size
  xpos = 0.5*(DINDGEN(nx)/DOUBLE(nx-1))^2 + 0.5*(DINDGEN(nx)/DOUBLE(nx-1))
  dx   = DINDGEN(nx)/DOUBLE(nx-1)/DOUBLE(nx-1) + 0.5/DOUBLE(nx-1)
  d = add_data(xpos, dx, nz, a=d)

  ; Quadratic change in mesh size
  xpos = 0.5 + 0.25*(2.*DINDGEN(nx)/DOUBLE(nx-1) - 1.)^3 + 0.25*(2.*DINDGEN(nx)/DOUBLE(nx-1) - 1.)
  dx   = 0.75*(2.*DINDGEN(nx)/DOUBLE(nx-1) - 1.)^2 * 2./DOUBLE(nx-1) + 0.5/DOUBLE(nx-1)
  d = add_data(xpos, dx, nz, a=d)
  
  ; Jump in mesh size
  xi = FIX(nx/2)
  xpos[0:(xi-1)] = 0.25*DINDGEN(xi)/DOUBLE(xi)
  xpos[xi:*]     = 0.25 + 0.75*DINDGEN(nx-xi)/DOUBLE(nx-xi-1)
  dx[0:(xi-1)] = 0.25 / DOUBLE(xi)
  dx[xi:*]     = 0.75 / DOUBLE(nx-xi-1)  
  d = add_data(xpos, dx, nz, a=d)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Write to file
  
  f = file_open("test_delp2.grd.nc", /create)
  
  status = file_write(f, "nx", d.nx)
  status = file_write(f, "ny", d.ny)
  
  status = file_write(f, "dx", d.dx)

  status = file_write(f, "xpos", d.xpos)
 
  status = file_write(f, "input", bout3dvar(d.input))
  status = file_write(f, "reference", bout3dvar(d.reference))
  
  file_close, f
END
