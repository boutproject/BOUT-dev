;
; Wrapper functions for I/O to different file formats
;
; Ben Dudson, University of York, July 2009
; 
; FUNCTION file_ndims, handle, varname
; FUNCTION file_read, handle, varname
; FUNCTION file_list, handle
; PRO file_close, handle
; FUNCTION file_open, name, write=write, create=create, format=format
;
; NOTES:
;  - file_open is the last function, so everything
;    is compiled when this is called.
;
;  - Currently uses the PDB2IDL library which can't handle more than one
;    open file at a time. To handle this, 
;    PDB files are opened and closed often.
;  - All functions need to be called through CALL_FUNCTION
;    otherwise won't compile if NCDF not supported,
;    or PDB2IDL not compiled
;
; Types: 1 = PDB, 2 = NetCDF
;
; 

; Open a file, returning a handle. Default is read-only.
; /write  -> Read/write mode
; /create -> Start a new file
; format=""  Specify a file format
FUNCTION file_open, name, write=write, create=create, format=format
  ON_ERROR, 2
  
  IF NOT KEYWORD_SET(format) THEN format = name ; use extension to guess format

  IF KEYWORD_SET(create) THEN write = 1

  IF STREGEX(format, '(nc|cdl|cdf|ncdf|netcdf)$', /BOOL) THEN BEGIN
    ;PRINT, "Opening NetCDF file ", name
    
    IF NOT CALL_FUNCTION('NCDF_EXISTS') THEN BEGIN
      PRINT, "  -> NetCDF not supported. Add to IDL_DLM_PATH."
      RETURN, 0
    ENDIF

    IF KEYWORD_SET(create) THEN BEGIN
      fp = CALL_FUNCTION('NCDF_CREATE', name, /CLOBBER)
    ENDIF ELSE BEGIN
      fp = CALL_FUNCTION('NCDF_OPEN', name, WRITE=write)
    ENDELSE
    
    ; turn off error messages
    CALL_PROCEDURE, 'NCDF_CONTROL', fp, /NOVERBOSE 

    RETURN, {type:2, name:name, id:fp}
    
  ENDIF ELSE BEGIN
    ; Assume it's a PDB file
    
    IF KEYWORD_SET(create) THEN BEGIN
      ; delete the file if it already exists
      SPAWN, "\ls "+name, result, exit_status=status
      IF status EQ 0 THEN BEGIN
        SPAWN, "rm "+ name
      ENDIF
      
    ENDIF

    ;PRINT, "Opening PDB file ", name
    status = CALL_FUNCTION('pdb_open', name, write=write) ; Opening as a check
    IF status EQ 0 THEN BEGIN
      PRINT, "  -> File not found"
      RETURN, 0
    ENDIF
    
    CALL_PROCEDURE, 'pdb_close'   ; Close file again
    
    IF NOT KEYWORD_SET(write) THEN write = 0

    RETURN, {type:1, write:write, name:name} ; PDB2IDL can't have lots of open files
  ENDELSE
  
  RETURN, 0
END
