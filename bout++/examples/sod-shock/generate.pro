; create an input file for the Sod shock-tube problem
; intended for use with the gas_compress.cpp physics module

PRO generate, ny=ny, mxg=mxg, file=file
  IF NOT KEYWORD_SET(ny) THEN ny = 128
  IF NOT KEYWORD_SET(mxg) THEN mxg = 2
  IF NOT KEYWORD_SET(file) THEN file="sod.grd.nc"

  nx = 1 + 2*FIX(mxg)

  density = FLTARR(nx, ny)
  pressure = FLTARR(nx, ny)

  dy = 1.0 / FLOAT(ny)
  
  dxarr = FLTARR(nx, ny) + 1.0;
  dyarr = FLTARR(nx, ny) + dy;
  
  yarr = FINDGEN(ny) * dy
  
  hy = FIX(ny / 2)
  
  FOR y=0, hy-1 DO BEGIN ; High pressure region
    density[*,y] = 1.0
    pressure[*,y] = 1.0
  ENDFOR
  FOR y=hy, ny-1 DO BEGIN ; Low pressure region
    density[*,y] = 0.125
    pressure[*,y] = 0.1
  ENDFOR
   
  ; Domain outside 'core' so non-periodic
  ixseps1 = 0
  ixseps2 = 0
  
  jyseps1_1 = -1
  jyseps2_2 = ny-1
  
  handle = file_open(file, /CREATE)
  
  s = file_write(handle, "nx", nx)
  s = file_write(handle, "ny", ny)
  s = file_write(handle, "dx", dxarr)
  s = file_write(handle, "dy", dyarr)

  s = file_write(handle, "ixseps1", ixseps1)
  s = file_write(handle, "ixseps2", ixseps2)
  s = file_write(handle, "jyseps1_1", jyseps1_1)
  s = file_write(handle, "jyseps2_2", jyseps2_2)
  
  s = file_write(handle, "density", density)
  s = file_write(handle, "pressure", pressure)
  
  file_close, handle
  
  PRINT, "DONE"
END
