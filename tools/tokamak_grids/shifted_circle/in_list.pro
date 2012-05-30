FUNCTION in_list, list, val, need=need
  ON_ERROR, 2

  IF SIZE(val, /TYPE) EQ 7 THEN BEGIN
      ; string
      list = STRUPCASE(list)
      val = STRUPCASE(val)
  ENDIF

  w = WHERE(list EQ val, count)
  IF count GT 0 THEN RETURN, 1

  IF KEYWORD_SET(need) THEN BEGIN
      PRINT, "ERROR: Need variable ", val
      STOP
  ENDIF
  RETURN, 0
END
