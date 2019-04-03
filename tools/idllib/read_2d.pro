FUNCTION read_2d, nx, ny, fid=fid
  IF KEYWORD_SET(fid) THEN status = next_double(fid=fid)
  
  d = DBLARR(nx, ny)
  
  FOR i=0, ny-1 DO d[*,i] = read_1d(nx)

  RETURN, d
END
