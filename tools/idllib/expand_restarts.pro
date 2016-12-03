;; Increases the number of Z points in restart files
;; NOTE: Can't over-write 
;
; path     Input path 
; output   Output path
; format   File extension of output

PRO expand_restarts, newz, path=path, output=output, format=format
  IF NOT KEYWORD_SET(path) THEN path = "data"
  IF NOT KEYWORD_SET(output) THEN output="./"

  IF STRCMP(path, output) THEN BEGIN
    PRINT, "Error: Can't overwrite restart files when resizing"
    RETURN
  ENDIF
  
  IF NOT is_pow2(newz-1) THEN BEGIN
    PRINT, "Error: New z size must be a power of 2 + 1"
    RETURN
  ENDIF
  
  SPAWN, "\ls "+path+"/BOUT.restart.*", result, exit_status=status

  IF status NE 0 THEN BEGIN
     PRINT, "ERROR: No data found"
     RETURN
  ENDIF

  ; Get the file extension 
  i = STRPOS(result[0], '.', /REVERSE_SEARCH)
  fext = STRMID(result[0], i+1, STRLEN(result[0])-(i+1))

  IF NOT KEYWORD_SET(format) THEN format = fext ; same as input format

  ; Select again only the ones with this extension
  SPAWN, "\ls "+path+"/BOUT.restart.*."+fext, result, exit_status=status

  nfiles = N_ELEMENTS(result)
  
  PRINT, "Number of files found: ", nfiles

  FOR fi=0, nfiles-1 DO BEGIN
    oldfile = path+"/BOUT.restart."+STRTRIM(STRING(fi),2)+"."+fext
    newfile = output+"/BOUT.restart."+STRTRIM(STRING(fi),2)+"."+format

    PRINT, "Changing "+oldfile + " => "+newfile
    ; open the restart file in read mode
    inhandle = file_open(oldfile)
    ; Create output file 
    outhandle = file_open(newfile, /create)

    ; get list of variables
    var_list = file_list(inhandle)
    
    nvars = N_ELEMENTS(var_list)
    FOR i=0, nvars-1 DO BEGIN
      
      var = var_list[i]
      
      ; Get number of dimensions
      nd = file_ndims(inhandle, var)
      ; read the data
      data = file_read(inhandle, var)
      
      ; find 3D variables
      IF nd EQ 3 THEN BEGIN
        PRINT, "   Resizing "+var
        
        s = SIZE(data, /dim)
        nx = s[0]
        ny = s[1]
        nz = s[2]

        newdata = FLTARR(nx, ny, newz)

        FOR x=0, nx-1 DO BEGIN
          FOR y=0, ny-1 DO BEGIN
            f = FFT(data[x,y,0:(nz-2)])

            ; Number of points in f is power of 2
            
            f2 = COMPLEXARR(newz-1)
            ; copy coefficients across (ignoring nyquist)
            f2[0] = f[0] ; DC
            minnz = MIN([nz, newz])
            FOR m=1, (minnz-1)/2 - 1 DO BEGIN
              f2[m] = f[m]       ; +ve frequencies
              f2[newz-1-m] = f[nz-1-m] ; -ve frequencies
            ENDFOR
            
            ; invert FFT
            
            newdata[x,y,0:(newz-2)] = REAL_PART(FFT(f2, /INV))
            newdata[x,y,newz-1] = newdata[x,y,0]

          ENDFOR
        ENDFOR
      ENDIF ELSE BEGIN
        PRINT, "   Copying "+var
        newdata = data
      ENDELSE

      status = file_write(outhandle, var, newdata)
    ENDFOR

    ; Close the files
    file_close, inhandle
    file_close, outhandle

  ENDFOR
  

END
