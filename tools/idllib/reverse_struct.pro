;******************************
; reverse_struct.pro
;******************************
; Created: 09/18/2012 I. Joseph
;******************************

FUNCTION reverse_struct, s
  IF size(s,/tname) eq 'STRUCT' THEN BEGIN
    tags = TAG_NAMES(s)
    n_tags = N_TAGS(s)
    sout = CREATE_STRUCT(tags(n_tags-1),s.(n_tags-1))
    FOR i=n_tags-2,0,-1 DO BEGIN
      sout = CREATE_STRUCT(sout,tags(i),s.(i))
    ENDFOR
    RETURN, sout
  ENDIF ELSE BEGIN 
    PRINT, 'ERROR: REVERSE_STR requires a STRUCT as input'
    RETURN, -1
  END
END
