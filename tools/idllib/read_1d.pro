FUNCTION read_1d, n, fid=fid
  IF KEYWORD_SET(fid) THEN status = next_double(fid=fid)
  
  d = DBLARR(n)
  
  FOR i=0, n-1 DO BEGIN
    d[i] = next_double()
  ENDFOR
  
  RETURN, d
END
