; Close the file
PRO file_close, handle
  ON_ERROR, 2
  IF handle.type EQ 1 THEN RETURN ; a PDB file 
  
  ; a netCDF file
    
  CALL_PROCEDURE, 'NCDF_CLOSE', handle.id
END
