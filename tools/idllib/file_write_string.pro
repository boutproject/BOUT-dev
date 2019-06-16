; Write a string to a file as an attribute
FUNCTION file_write_string, handle, varname, data

  CATCH, errcode
  
  IF errcode NE 0 THEN BEGIN
    ; error occurred
    PRINT, "Error occurred in file_write: "+!ERROR_STATE.MSG
    PRINT, "=> Couldn't write variable '"+varname+"'"
    RETURN, 1
  ENDIF

  IF handle.type EQ 1 THEN BEGIN
    ; PDB
    status = CALL_FUNCTION('pdb_open', handle.name, /write)
    
    status = CALL_FUNCTION('pdb_write_var', varname, data)

    CALL_PROCEDURE, 'pdb_close'
    
    CATCH, /CANCEL
    RETURN, status
  ENDIF

  ; NetCDF

  NCDF_CONTROL, handle.id, /REDEF

  NCDF_ATTPUT, handle.id, /GLOBAL, varname, data, /CHAR
    
  NCDF_CONTROL, handle.id, /ENDEF
  
  CATCH, /CANCEL
  
  RETURN, 0
END
