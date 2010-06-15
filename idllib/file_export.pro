; export a structure to a file
; return non-zero on error

FUNCTION file_export, file, data, lowercase=lowercase, uppercase=uppercase
  
  ON_ERROR, 2
  IF SIZE(file, /TYPE) EQ 7 THEN BEGIN ; a string, so filename
    handle = file_open(file, /create)
  ENDIF ELSE handle = file
  
  t = SIZE(data, /structure)
  
  IF (SIZE(handle, /TYPE) NE 8) OR (t.type NE 8) THEN BEGIN
    PRINT, "Useage: file_export(file, data)"
    RETURN, 1
  ENDIF

  names = TAG_NAMES(data)

  IF KEYWORD_SET(lowercase) THEN names = STRLOWCASE(names)
  IF KEYWORD_SET(uppercase) THEN names = STRUPCASE(names)

  FOR i=0, N_TAGS(data)-1 DO BEGIN
    var = data.(i)
    status = file_write(handle, names[i], var)
  ENDFOR

  IF SIZE(file, /TYPE) EQ 7 THEN file_close, handle

  RETURN, 0
END
