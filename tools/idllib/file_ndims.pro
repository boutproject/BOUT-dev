; Number of dimensions
FUNCTION file_ndims, file, varname
  ON_ERROR, 2
  
  IF SIZE(file, /TYPE) EQ 7 THEN BEGIN ; a string, so filename
    handle = file_open(file)
  ENDIF ELSE handle = file
  
  IF handle.type EQ 1 THEN BEGIN
    ; PDB
    
    status = CALL_FUNCTION('pdb_open', handle.name)
    query = CALL_FUNCTION('pdb_query_var', varname)
    CALL_PROCEDURE, 'pdb_close'
    
    RETURN, query[1]
  ENDIF
  
  ; NetCDF

  id = CALL_FUNCTION('NCDF_VARID', handle.id, varname) ; Get variable ID
  
  IF id EQ -1 THEN BEGIN
    PRINT, "Variable '"+varname+"' not found"
    RETURN, 0
  ENDIF

  info = CALL_FUNCTION('NCDF_VARINQ', handle.id, id)

  IF SIZE(file, /TYPE) EQ 7 THEN file_close, handle
  
  RETURN, info.ndims
END
