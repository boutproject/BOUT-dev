; List all variables in a file
FUNCTION file_list, file
  ON_ERROR, 2
  
  IF SIZE(file, /TYPE) EQ 7 THEN BEGIN ; a string, so filename
    handle = file_open(file)
  ENDIF ELSE handle = file

  IF handle.type EQ 1 THEN BEGIN
    status = CALL_FUNCTION('pdb_open', handle.name)
    
    var_list = CALL_FUNCTION('pdb_get_list')
    
    CALL_PROCEDURE, 'pdb_close'

    RETURN, var_list
  ENDIF

  ; a netCDF file
    
  desc = CALL_FUNCTION('NCDF_INQUIRE', handle.id )
  
  FOR i=0, desc.nvars-1 DO BEGIN
    info = CALL_FUNCTION('NCDF_VARINQ', handle.id, i)
    
    IF i EQ 0 THEN var_list = [info.name] ELSE var_list = [var_list, info.name]
  ENDFOR
  
  IF SIZE(file, /TYPE) EQ 7 THEN file_close, handle

  RETURN, var_list
END
