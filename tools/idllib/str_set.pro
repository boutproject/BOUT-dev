; Set a value in a structure
PRO str_set, str, name, value
  tag = TAG_NAMES(str)
  IF MAX(STRCMP(tag, name, /FOLD_CASE), ind) EQ 0 THEN BEGIN
    PRINT, " Adding "+name+" to structure"
    str=CREATE_STRUCT(str, name, value)
  ENDIF ELSE BEGIN
    ; Set the value
    str.(ind) = value
  ENDELSE
END
