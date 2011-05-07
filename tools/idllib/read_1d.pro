FUNCTION read_1d, n
  d = DBLARR(n)
  
  FOR i=0, n-1 DO BEGIN
    d[i] = next_double()
  ENDFOR
  
  RETURN, d
END
