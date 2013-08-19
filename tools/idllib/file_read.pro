; Read a variable from file
; 
; 
; file    Can be a file handle or name
; 
; KEYWORDS
;   inds = [xmin, xmax, ymin, ymax, ...]

FUNCTION file_read, file, varname, inds=inds
  ON_ERROR, 2
  
  IF SIZE(file, /TYPE) EQ 7 THEN BEGIN ; a string, so filename
    handle = file_open(file)
  ENDIF ELSE handle = file

  IF handle.type EQ 1 THEN BEGIN
    ; PDB file
    
    IF NOT KEYWORD_SET(inds) THEN BEGIN
      ; just read the whole variable
      data = CALL_FUNCTION('pd_read', handle.name, varname)
      IF SIZE(file, /TYPE) EQ 7 THEN file_close, handle
      RETURN, data
    ENDIF ELSE BEGIN

      status = CALL_FUNCTION('pdb_open', handle.name)
      
      ; first query
      query = CALL_FUNCTION('pdb_query_var', varname)

      nd = query[1]
      
      ; Change query ranges
      nm = MIN([FIX(N_ELEMENTS(inds) / 2), nd])
      query[2:(1+nm*2)] = inds[0:(2*nm - 1)]
      
      ; Read data
      data = CALL_FUNCTION('pdb_read_var', varname, query)
      
      CALL_PROCEDURE, 'pdb_close'
      
      IF SIZE(file, /TYPE) EQ 7 THEN file_close, handle
      RETURN, data
    ENDELSE
  ENDIF

  ; NetCDF
  
  id = CALL_FUNCTION('NCDF_VARID', handle.id, varname) ; Get variable ID
  
  IF id EQ -1 THEN BEGIN
    PRINT, "Variable '"+varname+"' not found"
    RETURN, 0
  ENDIF

  IF NOT KEYWORD_SET(inds) THEN BEGIN
    CALL_PROCEDURE, 'NCDF_VARGET', handle.id, id, value
  ENDIF ELSE BEGIN
    nm = FIX(N_ELEMENTS(inds)/2)
    offset = LONARR(nm)
    count = LONARR(nm)
    
    FOR i=0, nm-1 DO BEGIN
      offset[i] = inds[2*i]
      count[i] = inds[2*i + 1] - offset[i] + 1
    ENDFOR
    
    ; NOTE: NetCDF reverses order of indices.

    offset = REVERSE(offset)
    count = REVERSE(count)

    CALL_PROCEDURE, 'NCDF_VARGET', handle.id, id, value, $
      OFFSET=offset, COUNT=count

    ; ncdf_varget truncates last dimension if size 1
    value = REFORM(value, count)
  ENDELSE
  
  ; Check if variable is a string
  varinfo = CALL_FUNCTION("NCDF_VarInq", handle.id, id)
  IF varinfo.datatype EQ 'CHAR' THEN value = STRING(value)

  ; Reverse data. Not very efficient.
  value = reverse_inds(value)

  IF SIZE(file, /TYPE) EQ 7 THEN file_close, handle
  RETURN, value
END
