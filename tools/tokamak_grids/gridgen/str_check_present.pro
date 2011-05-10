; Check if a given field is in a structure, and if not then add it
; with a default value
PRO str_check_present, str, name, default
  tag = TAG_NAMES(str)
  IF MAX(STRCMP(tag, name, /FOLD_CASE)) EQ 0 THEN BEGIN
    PRINT, "  "+name+" not specified. Setting to "+STR(default)
    str=CREATE_STRUCT(str, name, default)
  ENDIF
END
