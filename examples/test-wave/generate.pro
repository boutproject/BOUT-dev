
PRO generate
  nx = 5
  ny = 64
  
  len = 1.
  
  dx = FLTARR(nx, ny) + 1.0
  dy = FLTARR(nx, ny) + len / FLOAT(ny)
  
  ixseps1 = LONG(5)
  ixseps2 = LONG(5)
  
  f = file_open('test_wave.grd.nc', /create)
  
  status = file_write(f, 'nx', nx)
  status = file_write(f, 'ny', ny)
  
  status = file_write(f, 'dx', dx)
  status = file_write(f, 'dy', dy)
  
  status = file_write(f, 'ixseps1', ixseps1)
  status = file_write(f, 'ixseps2', ixseps1)
  
  file_close, f
END
