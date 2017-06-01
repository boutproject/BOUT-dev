; create an input file for the Rayleigh-Taylor instability
; intended for use with the gas_compress.cpp physics module

PRO generate, ny=ny, mxg=mxg, file=file
  IF NOT KEYWORD_SET(ny) THEN ny = 128
  IF NOT KEYWORD_SET(mxg) THEN mxg = 2
  IF NOT KEYWORD_SET(file) THEN file="data/rt.grd.pdb"

  nx = 1 + 2*FIX(mxg)

  density = FLTARR(nx, ny)
  pressure = FLTARR(nx, ny)

  n2 = FIX(ny / 2)

  dy = 1.0 / FLOAT(ny - 1)
  
  dxarr = FLTARR(nx, ny) + 1.0;
  dyarr = FLTARR(nx, ny) + dy;
  
  yarr = FINDGEN(ny) * dy
  
  FOR y=0, n2-1 DO BEGIN
     density[*,y] = 2.0
     pressure[*,y] = 2.0*yarr[y] + 1.0
  ENDFOR

  FOR y=n2, ny-1 DO BEGIN
     density[*, y] = 1.0
     pressure[*,y] = yarr[y] + 1.5;
  ENDFOR
  
  gy = FLTARR(nx, ny) + 1.0;

  ; entire domain outside 'core': just a box
  ixseps1 = 0
  ixseps2 = 0
  
  PD_write, file, "nx", nx
  PD_write, file, "ny", ny
  PD_write, file, "dx", dxarr
  PD_write, file, "dy", dyarr

  PD_write, file, "ixseps1", ixseps1
  PD_write, file, "ixseps2", ixseps2

  PD_write, file, "density", density
  PD_write, file, "pressure", pressure
  PD_write, file, "gy", gy
END
