
FUNCTION file_nc_find_dim, handle, name, length
  ; Either find the given dimension, or create a new one

  COMPILE_OPT hidden
  ON_ERROR, 2
  
  did = CALL_FUNCTION('NCDF_DIMID', handle.id, name)
  
  IF did NE -1 THEN BEGIN
    ; Check the dimension is the right size
    
    CALL_PROCEDURE, 'NCDF_DIMINQ', handle.id, did, dname, size
    
    IF size NE length THEN BEGIN
      ; Find a dimension with the correct size
      did = -1
      
      prop = CALL_FUNCTION('NCDF_INQUIRE', handle.id)
      
      FOR i=0, prop.ndims-1 DO BEGIN
        CALL_PROCEDURE, 'NCDF_DIMINQ', handle.id, i, dname, size
        IF size EQ length THEN BEGIN
          did = i
          BREAK
        ENDIF
      ENDFOR
    ENDIF

    IF did EQ -1 THEN BEGIN
      ; Dimension <name> wrong size, and no other matching. Define a
      ; new one
      
      i = 2
      REPEAT BEGIN
        ; Append a number to the end of name. Check if it exists
        dname = name + STRTRIM(STRING(i),2)
        id = CALL_FUNCTION('NCDF_DIMID', handle.id, dname)
        i = i + 1
      ENDREP UNTIL id EQ -1
      ; define this dimension
      did = CALL_FUNCTION('NCDF_DIMDEF', handle.id, dname, length)
    ENDIF
  ENDIF ELSE BEGIN
    ; Define the dimension
    did = CALL_FUNCTION('NCDF_DIMDEF', handle.id, name, length)
  ENDELSE

  RETURN, did
END

; Write a variable to a file
FUNCTION file_write, handle, varname, data

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

  ; See if variable is already defined
  vid = CALL_FUNCTION('NCDF_VARID', handle.id, varname)

  IF vid EQ -1 THEN BEGIN
    ; Need to define the variable
    
    CALL_PROCEDURE, 'NCDF_CONTROL', handle.id, /REDEF

    nd = SIZE(data, /N_DIM)
    
    IF nd EQ 0 THEN BEGIN
      ; Just a scalar
      vartype = SIZE(data, /TYPE)
      CASE vartype OF
        1: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, /BYTE)
        2: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, /SHORT)
        3: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, /LONG)
        4: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, /FLOAT)
        5: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, /DOUBLE)
        7: BEGIN ; String
           dim = file_nc_find_dim(handle, "strlen", STRLEN(data))
           vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, dim, /CHAR)
        END
        12: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, /LONG)
        13: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, /LONG)
        14: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, /LONG)
        15: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, /LONG)
        ELSE: BEGIN
          PRINT, "Sorry, file_write doesn't support this type: ", vartype
          PRINT, "=> Couldn't write variable '"+varname+"'"
          RETURN, 1
        ENDELSE
      ENDCASE
      IF vid LT 0 THEN STOP
    ENDIF ELSE BEGIN
      dimsize = REVERSE(SIZE(data, /DIMEN)) ; Reverse indices
      
      ; Need to match up dimension IDs
      
      ; 2D - yx
      ; 3D - zyx
      ; 4D - zyxt
      
      did = INTARR(nd) ; Dimension IDs

      CASE nd OF
        1: BEGIN
          did[0] = file_nc_find_dim(handle, "x", dimsize[0])
        END
        2: BEGIN
          did[0] = file_nc_find_dim(handle, "y", dimsize[0])
          did[1] = file_nc_find_dim(handle, "x", dimsize[1])
        END
        3: BEGIN
          did[0] = file_nc_find_dim(handle, "z", dimsize[0])
          did[1] = file_nc_find_dim(handle, "y", dimsize[1])
          did[2] = file_nc_find_dim(handle, "x", dimsize[2])
        END
        4: BEGIN
          did[0] = file_nc_find_dim(handle, "z", dimsize[0])
          did[1] = file_nc_find_dim(handle, "y", dimsize[1])
          did[2] = file_nc_find_dim(handle, "x", dimsize[2])
          did[3] = file_nc_find_dim(handle, "t", dimsize[3])
        END
        ELSE: BEGIN
          PRINT, "Sorry, file_write doesn't support this many dimensions: ", nd
          PRINT, "=> Couldn't write variable '"+varname+"'"
          RETURN, 1
        END
      ENDCASE

      ; Define variable
      
      vartype = SIZE(data, /TYPE)
      CASE vartype OF
        1: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, did, /BYTE)
        2: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, did, /SHORT)
        3: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, did, /LONG)
        4: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, did, /FLOAT)
        5: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, did, /DOUBLE)
        7: BEGIN ; Array of strings
           PRINT, "Sorry, arrays of strings not supported"
           RETURN, 1
        END
        12: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, /LONG)
        13: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, /LONG)
        14: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, /LONG)
        15: vid = CALL_FUNCTION('NCDF_VARDEF', handle.id, varname, /LONG)
        ELSE: BEGIN
          PRINT, "Sorry, file_write doesn't support this type: ", vartype
          PRINT, "=> Couldn't write variable '"+varname+"'"
          RETURN, 1
        ENDELSE
      ENDCASE

    ENDELSE
    
  ENDIF
  
  CALL_PROCEDURE, 'NCDF_CONTROL', handle.id, /ENDEF
  
  ; Reverse indices of data and write to the file

  vartype = SIZE(data, /TYPE)
  IF (vartype GE 12) AND (vartype LE 15) THEN BEGIN
    ; unsigned long and 64 integer types
    data2 = LONG(data)
    CALL_PROCEDURE, 'NCDF_VARPUT', handle.id, vid, REVERSE_INDS(data2)
  ENDIF ELSE BEGIN
    CALL_PROCEDURE, 'NCDF_VARPUT', handle.id, vid, REVERSE_INDS(data)
  ENDELSE
  
  CATCH, /CANCEL
  
  RETURN, 0
END
