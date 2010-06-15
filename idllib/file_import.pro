; Import an entire data file
;
; file   - Either a file name, or a file handle

FUNCTION file_import, file, _extra=_extra
  ON_ERROR, 2
  
  IF SIZE(file, /TYPE) EQ 7 THEN BEGIN ; a string, so filename
    handle = file_open(file, _extra=_extra)
  ENDIF ELSE handle = file

  IF SIZE(handle, /TYPE) NE 8 THEN BEGIN
    PRINT, "Cannot import"
    RETURN, 0
  ENDIF

  vars = file_list(handle)
  nvars = N_ELEMENTS(vars)
  FOR i=0, nvars-1 DO BEGIN
    data = file_read(handle, vars[i])

    ; Need to convert to a legal name
    label = IDL_VALIDNAME(vars[i], /CONVERT_ALL)
    
    IF i EQ 0 THEN BEGIN
      str = CREATE_STRUCT(label, data)
    ENDIF ELSE str = CREATE_STRUCT(label, data, str)
    
  ENDFOR
  
  IF SIZE(file, /TYPE) EQ 7 THEN file_close, handle

  RETURN, str
END
