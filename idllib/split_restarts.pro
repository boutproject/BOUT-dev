; split the grid between a different number of processors
; Useage:
;   split_restarts, [nxpe], nype
; keywords
;   path = "path/to/input/"
;   output = "path/to/output"
;   format = "fext"  file extension of output ("pdb" or "nc")

PRO split_restarts, nxpe, nype, path=path, output=output, format=format
  IF NOT KEYWORD_SET(path) THEN path = "data"
  IF NOT KEYWORD_SET(output) THEN output="./"

  MXG = 2
  MYG = 2

  np = N_PARAMS()  ; number of parameters
  IF np EQ 0 THEN BEGIN
    PRINT, "Error: Need to specify number of output processors"
    RETURN
  ENDIF ELSE IF np EQ 1 THEN BEGIN
    ; specified only nype
    nype = nxpe
    nxpe = 1
  ENDIF

  npes = nxpe * nype

  IF npes LE 0 THEN BEGIN
    PRINT, "ERROR: Rediculous nxpe or nype specified!"
    RETURN
  ENDIF
  
  IF STRCMP(path, output) THEN BEGIN
    PRINT, "Error: Can't overwrite restart files when resizing"
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

  ; Select again only the ones with this extension
  SPAWN, "\ls "+path+"/BOUT.restart.*."+fext, result, exit_status=status

  nfiles = N_ELEMENTS(result)
  
  PRINT, "Number of files found: ", nfiles

  IF NOT KEYWORD_SET(format) THEN format = fext

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Read old processor layout from restart 0

  file = path+"/BOUT.restart.0."+fext
    
  ; open in read-only mode
  handle = file_open(file)
  
  ; get list of variables
  var_list = file_list(handle)
  
  ; Get size of the grid, and processor layout
  old_npes = file_read(handle, "NPES")
  old_nxpe = file_read(handle, "NXPE")

  file_close, handle

  IF nfiles NE old_npes THEN BEGIN
    PRINT, "WARNING: number of restart files inconsistent with NPES"
    PRINT, "Setting nfiles = ", old_npes
    
    nfiles = old_npes
  ENDIF

  IF old_npes MOD old_nxpe NE 0 THEN BEGIN
    PRINT, "ERROR: old NPES is not a multiple of old NXPE. This makes no sense."
    RETURN
  ENDIF

  old_nype = old_npes / old_nxpe
  
  IF nype MOD old_nype NE 0 THEN BEGIN
    PRINT, "SORRY: Currently can only split restart files in Y"
    RETURN
  ENDIF

  IF nxpe MOD old_nxpe NE 0 THEN BEGIN
    PRINT, "SORRY: Currently can only split restart files in X"
    RETURN
  ENDIF
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Get variable details

  nx = 0  ; Don't know the X and Y resolution yet

  n3d = 0 ; number of 3D variables
  n2d = 0 ; number of 2D variables
  nother = 0 ; number of other variables
  
  handle = file_open(file)

  nvars = N_ELEMENTS(var_list)
  FOR i=0, nvars-1 DO BEGIN
    ; query the variable type
    ndims = file_ndims(handle, var_list[i])
    
    IF ndims EQ 3 THEN BEGIN
      ; 3D variable
      
      IF n3d EQ 0 THEN var_3d = var_list[i] ELSE var_3d = [var_3d, var_list[i]]
      n3d = n3d + 1
      
      IF nx EQ 0 THEN BEGIN
        ; need to find out how big the grid is
        ; Simple way: read in the data
        
        data = file_read(handle, var_list[i])
        s = SIZE(data, /dim)
        
        old_MXSUB = s[0] - 2*MXG
        old_MYSUB = s[1] - 2*MYG
        MZ = s[2]
        
        ; calculate total size of the grid
        NX = old_MXSUB * old_NXPE
        NY = old_MYSUB * old_NYPE
        PRINT, "Grid size: ", NX, NY, MZ
        
        ; get new subgrid size
        IF (NX MOD NXPE NE 0) OR (NY MOD NYPE NE 0) THEN BEGIN
          PRINT, "ERROR: Cannot divide grid equally across processors"
          RETURN
        ENDIF
        MXSUB = NX / NXPE
        MYSUB = NY / NYPE
        
      ENDIF

    ENDIF ELSE IF ndims EQ 2 THEN BEGIN
      IF n2d EQ 0 THEN var_2d = var_list[i] ELSE var_2d = [var_2d, var_list[i]]
      n2d = n2d + 1
    ENDIF ELSE IF (STRCMP(var_list[i], "NPES") EQ 0) AND (STRCMP(var_list[i], "NXPE") EQ 0) THEN BEGIN
      ; other variables, ignoring NPES and NXPE
      
      IF nother EQ 0 THEN var_other = var_list[i] ELSE var_other = [var_other, var_list[i]]
      nother = nother + 1
      
    ENDIF
  ENDFOR

  file_close, handle

  PRINT, " 3D variables: ", n3d
  PRINT, " 2D variables: ", n2d
  PRINT, " Others      : ", nother

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Create new restart files

  xs = NXPE / old_NXPE ; number of x slices
  ys = NYPE / old_NYPE ; number of y slices

  FOR MYPE = 0, NPES-1 DO BEGIN
    ; calculate X and Y processor numbers
    PEX = MYPE MOD NXPE
    PEY = FIX(MYPE / NXPE)
    
    old_PEX = FIX(PEX / xs)
    old_PEY = FIX(PEY / ys)
    old_x = PEX MOD xs
    old_y = PEY MOD ys

    old_MYPE = old_NXPE * old_PEY + old_PEX ; get old restart file number

    ; calculate indices in old restart file
    xmin = old_x*MXSUB
    xmax = xmin + MXSUB - 1 + 2*MXG
    ymin = old_y*MYSUB
    ymax = ymin + MYSUB - 1 + 2*MYG
    
    PRINT, MYPE, PEX, PEY
    PRINT, "=>", old_MYPE, old_PEX, old_PEY, ":", old_x, old_y

    ; get filenames
    file = output+"/BOUT.restart."+STRTRIM(STRING(MYPE),2)+"."+format
    oldfile = path+"/BOUT.restart."+STRTRIM(STRING(old_MYPE),2)+"."+fext

    inhandle = file_open(oldfile)
    outhandle = file_open(file, /create)

    ; go through 3D variables
    FOR i=0, n3d-1 DO BEGIN
      ; read the original data
      data = file_read(inhandle, var_3d[i])
      
      ; chop out the required range
      data = data[xmin:xmax, ymin:ymax, *]
      
      ; write out the new data
      status = file_write(outhandle, var_3d[i], data)
    ENDFOR

    ; go through 2D variables
    FOR i=0, n2d-1 DO BEGIN
      ; read the original data
      data = file_read(inhandle, var_2d[i])
      
      ; chop out the required range
      data = data[xmin:xmax, ymin:ymax]
      
      ; write out the new data
      status = file_write(outhandle, var_2d[i], data)
    ENDFOR

    ; copy the other variables across
    FOR i=0, nother-1 DO BEGIN
      data = file_read(inhandle, var_other[i])
      
      status = file_write(outhandle, var_other[i], data)
    ENDFOR

    ; write NPES and NXPE
    status = file_write(outhandle, "NPES", NPES)
    status = file_write(outhandle, "NXPE", NXPE)

    file_close, inhandle
    file_close, outhandle
  ENDFOR
END
