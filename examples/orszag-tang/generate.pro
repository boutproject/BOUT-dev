; generates an input file for Orszag-Tang vortex problem

PRO generate, n=n, mxg=mxg, file=file
  IF NOT KEYWORD_SET(n) THEN n = 128
  IF NOT KEYWORD_SET(mxg) THEN mxg = 2
  IF NOT KEYWORD_SET(file) THEN file="data/otv.grd.nc"

  ny = n
  nx = ny + 2*FIX(mxg)

  dy = 1.0 / FLOAT(ny)
  dx = dy

  dxarr = FLTARR(nx, ny) + dx
  dyarr = FLTARR(nx, ny) + dy
  
  xarr = FINDGEN(nx) * dx
  yarr = FINDGEN(ny) * dy
  
  rho = FLTARR(nx, ny) + 25./(36.*!PI)
  P = FLTARR(nx, ny) + 5./(12.*!PI)
  
  v_x = FLTARR(nx, ny)
  Bx = FLTARR(nx, ny)

  FOR y=0, ny-1 DO BEGIN
      v_x[*,y] = -SIN(2.*!PI*yarr[y])
      Bx[*,y]  = -SIN(2.*!PI*yarr[y])
  ENDFOR
  Bx = Bx / SQRT(4.*!PI)

  v_y = FLTARR(nx, ny)
  By = FLTARR(nx, ny)

  FOR x=0, nx-1 DO BEGIN
      v_y[x,*] = SIN(2.*!PI*xarr[x])
      By[x,*] = SIN(4.*!PI*xarr[x])
  ENDFOR
  By = By / SQRT(4.*!PI)

  ; domain inside core (periodic)

  ixseps1 = nx  ; Index of separatrix 1 and 2
  ixseps2 = nx

  f = file_open(file, /create)
  status = file_write(f, "nx", nx)
  status = file_write(f, "ny", ny)
  status = file_write(f, "dx", dxarr)
  status = file_write(f, "dy", dyarr)

  status = file_write(f, "ixseps1", ixseps1)
  status = file_write(f, "ixseps2", ixseps2)
  
  status = file_write(f, "rho0", rho)
  status = file_write(f, "p0", P)
  status = file_write(f, "v0_x", v_x)
  status = file_write(f, "v0_y", v_y)
  status = file_write(f, "B0x", Bx)
  status = file_write(f, "B0y", By)
  
  file_close, f
END
