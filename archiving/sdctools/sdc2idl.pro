
FUNCTION sdc_read, filename, range=range

  ; open the file
  status = CALL_EXTERNAL('sdc2idl.so', 'sdc2idl_open', filename, /i_value)
  IF status NE 0 THEN RETURN, 0

  ; get the name of the variable
  name = CALL_EXTERNAL('sdc2idl.so', 'sdc2idl_get_name', /s_value)

  ; get number of dimensions
  nd = CALL_EXTERNAL('sdc2idl.so', 'sdc2idl_get_ndim', /i_value)

  dims = LONARR(2*nd)
  
  status = CALL_EXTERNAL('sdc2idl.so', 'sdc2idl_get_dims', dims, /i_value)
  IF status NE 0 THEN RETURN, 0

  PRINT, "VARIABLE NAME: ", name
  PRINT, "NUMBER OF DIMENSIONS: ", nd

  
  ds = LONARR(nd)
  size = 1
  FOR i=0,nd-2 DO BEGIN
      size = size * (dims[2*i+1] - dims[2*i] + 1)
      ds[i] = dims[2*i+1] - dims[2*i] + 1
  ENDFOR
  nt = dims[2*(nd-1)+1] - dims[2*(nd-1)] + 1
  ds[nd-1] = nt
  PRINT, "Number of elements per time-slice: ", size

  mint = 0L
  maxt = LONG(nt-1)

  IF KEYWORD_SET(range) THEN BEGIN
      range = LONG(range)
      IF range[0] GT mint THEN mint = range[0]
      IF range[1] LT maxt THEN maxt = range[1]
  ENDIF

  nt = maxt - mint + 1
  var = FLTARR(size*nt)

  PRINT, "Time range: ", mint, " to ", maxt

  status = CALL_EXTERNAL('sdc2idl.so', 'sdc2idl_read', mint, maxt, var, /i_value)
  IF status NE 0 THEN RETURN, 0
  

  status = CALL_EXTERNAL('sdc2idl.so', 'sdc2idl_close', /i_value)
  IF status NE 0 THEN RETURN, 0
  
  IF nd EQ 1 THEN BEGIN
      ; no change
      output = var
  ENDIF ELSE IF nd EQ 2 THEN BEGIN
      output = reform(var, size, nt)
  ENDIF ELSE IF nd EQ 3 THEN BEGIN
      nx = ds[0]
      ny = ds[1]
      output = FLTARR(nx, ny, nt)
      
      var = reform(var, size, nt)
      FOR t=0, nt-1 DO BEGIN
          p = 0
          FOR i=0, nx-1 DO BEGIN
              FOR j=0, ny-1 DO BEGIN
                  output[i,j,t] = var[p,nt]
                  p = p + 1
              ENDFOR
          ENDFOR
      ENDFOR

  ENDIF ELSE IF nd EQ 4 THEN BEGIN
      nx = ds[0]
      ny = ds[1]
      nz = ds[2]
      output = FLTARR(nx, ny, nz, nt)
      
      var = reform(var, size, nt)
      FOR t=0, nt-1 DO BEGIN
          p = 0L
          FOR i=0, nx-1 DO BEGIN
              FOR j=0, ny-1 DO BEGIN
                  FOR k=0,nz-1 DO BEGIN
                      output[i,j,k,t] = var[p,t]
                      p = p + 1L
                  ENDFOR
              ENDFOR
          ENDFOR
      ENDFOR

  ENDIF ELSE BEGIN
      PRINT, "Cannot handle this number of dimensions yet"
      RETURN,0
  ENDELSE
  
  
  RETURN, output
END
