; Return array of dimension sizes

FUNCTION file_size, file, varname
  ON_ERROR, 2

  IF SIZE(file, /TYPE) EQ 7 THEN BEGIN ; a string, so filename
    handle = file_open(file)
  ENDIF ELSE handle = file

  IF handle.type EQ 1 THEN BEGIN
    ; PDB
    
    status = CALL_FUNCTION('pdb_open', handle.name)
    query = CALL_FUNCTION('pdb_query_var', varname)
    CALL_PROCEDURE, 'pdb_close'
    
    nd = query[1]
    IF nd EQ 0 THEN RETURN, 1
    dims = LONARR(nd)
    FOR i=0, nd-1 DO BEGIN
      dims[i] = query[3 + 2*i] - query[2 + 2*i] + 1
    ENDFOR

    RETURN, dims
  ENDIF

  ; NetCDF

  id = CALL_FUNCTION('NCDF_VARID', handle.id, varname) ; Get variable ID
  
  IF id EQ -1 THEN BEGIN
    PRINT, "Variable '"+varname+"' not found"
    RETURN, 0
  ENDIF

  info = CALL_FUNCTION('NCDF_VARINQ', handle.id, id)

  IF info.ndims EQ 0 THEN RETURN, 1

  dims = LONARR(info.ndims)
  FOR i=0, info.ndims-1 DO BEGIN
    CALL_PROCEDURE, 'NCDF_DIMINQ', handle.id, info.dim[i], name, len
    dims[i] = len
  ENDFOR

  IF SIZE(file, /TYPE) EQ 7 THEN file_close, handle

  RETURN, dims
END
