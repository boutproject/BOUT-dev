; create an input file for the Rayleigh-Taylor instability
; intended for use with the gas_compress.cpp physics module

PRO generate, ny=ny, mxg=mxg, file=file
  IF NOT KEYWORD_SET(ny) THEN ny = 128
  IF NOT KEYWORD_SET(mxg) THEN mxg = 2
  IF NOT KEYWORD_SET(file) THEN file="data/advect.grd.pdb"

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

  PD_write, file, "nx", nx
  PD_write, file, "ny", ny
  PD_write, file, "dx", dxarr
  PD_write, file, "dy", dyarr

  PD_write, file, "ixseps1", ixseps1
  PD_write, file, "ixseps2", ixseps2
  PD_write, file, "jyseps1_1", jyseps1_1
  PD_write, file, "jyseps2_2", jyseps2_2

  PD_write, file, "density", density
  PD_write, file, "pressure", pressure
  PD_write, file, "vy", vy
END
