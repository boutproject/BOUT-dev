; create an input file for the Rayleigh-Taylor instability
; intended for use with the gas_compress.cpp physics module

PRO generate, ny=ny, mxg=mxg, file=file
  IF NOT KEYWORD_SET(ny) THEN ny = 128
  IF NOT KEYWORD_SET(mxg) THEN mxg = 2
  IF NOT KEYWORD_SET(file) THEN file="advect.grd.nc"

  nx = 1 + 2*FIX(mxg)

  density = FLTARR(nx, ny)
  pressure = FLTARR(nx, ny)
  vy = FLTARR(nx, ny)

  dy = 1.0 / FLOAT(ny)
  
  dxarr = FLTARR(nx, ny) + 1.0;
  dyarr = FLTARR(nx, ny) + dy;
  
  yarr = FINDGEN(ny) * dy
  
  FOR y=0, ny-1 DO BEGIN
     density[*,y] = 1.0
     pressure[*,y] = 3.0/5.0   ; so gamma*P = 1 => sound speed = 1
     vy[*,y] = 1.0 ; change in BOUT.inp 
 ENDFOR

  ; entire domain inside 'core' - periodic
  ixseps1 = nx
  ixseps2 = nx
  
  jyseps1_1 = -1
  jyseps2_2 = ny-1

  fp = file_open(file, /create)

  status = file_write(fp, "nx", nx)
  status = file_write(fp, "ny", ny)
  status = file_write(fp, "dx", dxarr)
  status = file_write(fp, "dy", dyarr)
  
  status = file_write(fp, "ixseps1", ixseps1)
  status = file_write(fp, "ixseps2", ixseps2)
  status = file_write(fp, "jyseps1_1", jyseps1_1)
  status = file_write(fp, "jyseps2_2", jyseps2_2)

  status = file_write(fp, "density", density)
  status = file_write(fp, "pressure", pressure)
  status = file_write(fp, "vy", vy)
  
  file_close, fp
END
