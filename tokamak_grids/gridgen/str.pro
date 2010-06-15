
FUNCTION STR, val
  n = N_ELEMENTS(val)
  IF n GT 1 THEN BEGIN
    s = "(" + STRTRIM(STRING(val[0]),2)
    FOR i=1, n-1 DO s = s + ", " + STRTRIM(STRING(val[i]),2)
    s = s + ")"
    RETURN, s
  ENDIF ELSE RETURN, STRTRIM(STRING(val),2)
END

