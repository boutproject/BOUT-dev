

PRO scale_restarts, factor, path=path, format=format
  IF NOT KEYWORD_SET(path) THEN path = "./"
  
  SPAWN, "\ls "+path+"/BOUT.restart.*", result, exit_status=status

  IF status NE 0 THEN BEGIN
     PRINT, "ERROR: No data found"
     RETURN
  ENDIF

  ; Get the file extension 
  i = STRPOS(result[0], '.', /REVERSE_SEARCH)
  fext = STRMID(result[0], i+1, STRLEN(result[0])-(i+1))

  ; Select again only the ones with this extension
  SPAWN, "\ls "+path+"/BOUT.restart.*."+fext, result, exit_status=status

  nfiles = N_ELEMENTS(result)
  
  PRINT, "Number of files found: ", nfiles

  FOR f=0, nfiles-1 DO BEGIN
    PRINT, "Changing "+result[f]
    ; open the restart file in read-write mode
    handle = file_open(result[f], /write, format=format)
    ; get list of variables
    var_list = file_list(handle)
    
    nvars = N_ELEMENTS(var_list)
    FOR i=0, nvars-1 DO BEGIN
      var = var_list[i]
      nd = file_ndims(handle, var)
      
      ; find 3D variables
      IF nd EQ 3 THEN BEGIN
        PRINT, "   Scaling "+var
        data = file_read(handle, var)
        
        data = data * factor
        
        status = file_write(handle, var, data)
      ENDIF
    ENDFOR
    
    file_close, handle
  ENDFOR
END
