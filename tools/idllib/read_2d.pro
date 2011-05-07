FUNCTION read_2d, nx, ny
  d = FLTARR(nx, ny)
  FOR i=0, ny-1 DO d[*,i] = read_1d(nx)

  RETURN, d
END
