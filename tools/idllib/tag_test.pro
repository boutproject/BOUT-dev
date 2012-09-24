;******************************
; tag_test.pro
;******************************
; Created: 03/26/2012 I. Joseph
;******************************

FUNCTION tag_test, s, tag
  tagname=STRTRIM(STR(tag))
  IF size(s,/tname) eq 'STRUCT' THEN BEGIN
    tags = TAG_NAMES(s)
    FOR i=1,N_TAGS(s) DO BEGIN
	IF tags(i-1) EQ tagname THEN RETURN, i
    ENDFOR
    PRINT, 'TAG_TEST: '+tagname+' not found'
    RETURN, -2
  ENDIF ELSE BEGIN 
    PRINT, 'ERROR: TAG_TEST requires a STRUCT as input'
    RETURN, -1
  END
END
