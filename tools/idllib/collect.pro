; Collects a set of BOUT++ outputs. Converts into BOUT's XYZT format
;
; Jan 2014: Updated variable name handling.
;         o Tries to match abbreviations as well as case variation
;         o Allows variable name to be passed as a parameter
;           or a keyword
; 
; August 2008: Updated to read BOUT++ v0.3 output, which is split in X
;              and Y

FUNCTION in_arr, arr, val

  w = WHERE(arr EQ val, count)
  IF count GT 0 THEN RETURN, 1

  RETURN, 0
END

FUNCTION collect, arg, xind=xind, yind=yind, zind=zind, tind=tind, $
             path=path, var=var, t_array=t_array, use=use, old=old,  $
                  quiet=quiet, debug=debug, prefix=prefix

  MESSAGE, "This is currently broken for BOUT++ > v4.0.0. See issue #394"

  IF NOT KEYWORD_SET(prefix) THEN prefix="BOUT.dmp"

  IF NOT KEYWORD_SET(debug) THEN BEGIN
    on_error,2  ; If an error occurs, return to caller

    ; catch any error which occurs
    CATCH, errcode
  
    IF errcode NE 0 THEN BEGIN
      ; error occurred
      PRINT, "Error occurred in collect: "+!ERROR_STATE.MSG
      
      RETURN, 0
    ENDIF
  ENDIF

  ; first find out how many files there are

  IF NOT KEYWORD_SET(use) THEN use = 0

  IF NOT KEYWORD_SET(path) THEN path = "./"

  IF KEYWORD_SET(quiet) THEN quiet = 2 ELSE quiet = 1
  IF KEYWORD_SET(debug) THEN quiet = 0

  IF N_PARAMS() GT 0 THEN var = arg

  IF NOT KEYWORD_SET(var) THEN BEGIN
    var = "Ni"
    PRINT, "No variable specified. Trying to read "+var
  ENDIF

  ;; Allow user to specify a single number for the index

  IF KEYWORD_SET(xind) THEN BEGIN
      IF N_ELEMENTS(xind) EQ 1 THEN xind = [xind, xind]
  ENDIF
  IF KEYWORD_SET(yind) THEN BEGIN
      IF N_ELEMENTS(yind) EQ 1 THEN yind = [yind, yind]
  ENDIF
  IF KEYWORD_SET(zind) THEN BEGIN
      IF N_ELEMENTS(zind) EQ 1 THEN zind = [zind, zind]
  ENDIF
  IF KEYWORD_SET(tind) THEN BEGIN
      IF N_ELEMENTS(tind) EQ 1 THEN tind = [tind, tind]
  ENDIF

  SPAWN, "\ls "+path+"/"+prefix+".*", result, exit_status=status

  IF status NE 0 THEN BEGIN
      PRINT, "ERROR: No data found"
      RETURN, 0
  ENDIF

  nfiles = N_ELEMENTS(result)
  
  IF quiet EQ 0 THEN PRINT, "Number of files found: ", nfiles
  
  mainfile = result[use] ; File to get most data

  ; Get the file extension 
  i = STRPOS(mainfile, '.', /REVERSE_SEARCH)
  fext = STRMID(mainfile, i+1, STRLEN(mainfile)-(i+1))

  ; Select again only the ones with this extension
  SPAWN, "\ls "+path+"/"+prefix+".*."+fext, result, exit_status=status

  nfo = nfiles
  nfiles = N_ELEMENTS(result)
  
  IF (quiet EQ 0) AND (nfo NE nfiles) THEN PRINT, "Number of files with same type: ", nfiles

  ; get list of variables
  handle = file_open(mainfile)
  var_list = file_list(handle)
  
  ; Check if the variable requested is available

  IF MAX(STRCMP(var_list, var)) EQ 0 THEN BEGIN
    PRINT, "Variable '"+var+"' not found"
    ; Check if the case is different
    res = STRCMP(var_list, var, /FOLD_CASE)
    w = WHERE(res EQ 1, count)
    IF COUNT GT 0 THEN BEGIN
      var = var_list[w[0]]
      PRINT, "-> Variables are case-sensitive: Using '"+var+"'"
    ENDIF ELSE BEGIN
      
      ; Check if it's an abbreviation
      success = 0
      Nmax = STRLEN(var)
      FOR N=Nmax, 1,-1 DO BEGIN
        res = STRCMP(var_list, var, N, /FOLD_CASE)
        w = WHERE(res EQ 1, count)
       
        IF count GE 2 THEN BEGIN
          PRINT, "Variable name is ambiguous. Could match :", var_list[w]
          file_close, handle
          RETURN, 0
        ENDIF ELSE IF count EQ 1 THEN BEGIN
          ; Matched
          var = var_list[w[0]]
          PRINT, "-> Using '"+var+"'"
          success = 1
          BREAK
        ENDIF
      ENDFOR
      
      IF success EQ 0 THEN BEGIN
        file_close, handle
        RETURN, 0
      ENDIF
    ENDELSE
  ENDIF

  var_list = STRUPCASE(var_list)  
  ; get the version of BOUT++

  IF in_arr(var_list, "BOUT_VERSION") THEN version = file_read(handle, "BOUT_VERSION") ELSE VERSION= 0.2

  IF quiet EQ 0 THEN PRINT, "BOUT++ VERSION: ", version

  ; read parameters
  mxsub = file_read(handle, "MXSUB")
  mysub = file_read(handle, "MYSUB")
  MZ    = file_read(handle, "MZ")
  myg   = file_read(handle, "MYG")
  t_array = file_read(handle, "t_array")
  nt = N_ELEMENTS(t_array)

  IF version GT 0.25 THEN BEGIN
    ; 2D decomposition - need to read NXPE and MXG
    
    NXPE = file_read(handle, "NXPE")
    MXG  = file_read(handle, "MXG")
    NYPE = file_read(handle, "NYPE")
    NPE  = NXPE*NYPE

    IF NPE LT nfiles THEN BEGIN
        PRINT, "WARNING: More files than expected"
        nfiles = NPE
    ENDIF ELSE IF NPE GT nfiles THEN BEGIN
        PRINT, "WARNING: Some files missing"
    ENDIF

    IF quiet EQ 0 THEN PRINT, "MXG = ", mxg

    NX = NXPE * MXSUB + 2*MXG
  ENDIF ELSE BEGIN
    NX = MXSUB
    MXG = 0
    NXPE = 1
    NYPE = nfiles
  ENDELSE

  IF version LT 3.5 THEN BEGIN
    nz = MZ-1  ; Remove extra point
  ENDIF ELSE nz = MZ

  NY = MYSUB * NYPE

  IF quiet EQ 0 THEN BEGIN
      PRINT, "Size of the grid: ", NX, NY, MZ
      PRINT, "In each file: ", MXSUB, MYSUB, MZ
  ENDIF

  ; check that var is 4-D

  ndims = file_ndims(handle, var)
  dimsize = file_size(handle, var)

  file_close, handle ; Close the result[use] file

  IF NOT KEYWORD_SET(xind) THEN xind = [0, nx-1]
  IF NOT KEYWORD_SET(yind) THEN yind = [0, ny-1]
  IF NOT KEYWORD_SET(zind) THEN zind = [0, nz-1]
  IF NOT KEYWORD_SET(tind) THEN tind = [0, nt-1]
  
  ; check ranges
  
  IF xind[0] LT 0 THEN xind[0] = 0
  IF xind[0] GT nx-1 THEN xind[0] = nx-1
  IF xind[1] LT 0 THEN xind[1] = 0
  IF xind[1] GT nx-1 THEN xind[1] = nx-1
  IF xind[0] GT xind[1] THEN BEGIN
    tmp = xind[0]
    xind[0] = xind[1]
    xind[1] = tmp
  ENDIF
    
  IF yind[0] LT 0 THEN yind[0] = 0
  IF yind[0] GT ny-1 THEN yind[0] = ny-1
  IF yind[1] LT 0 THEN yind[1] = 0
  IF yind[1] GT ny-1 THEN yind[1] = ny-1
  IF yind[0] GT yind[1] THEN BEGIN
    tmp = yind[0]
    yind[0] = yind[1]
    yind[1] = tmp
  ENDIF
  
  IF zind[0] LT 0 THEN zind[0] = 0
  IF zind[0] GT nz-1 THEN zind[0] = nz-1
  IF zind[1] LT 0 THEN zind[1] = 0
  IF zind[1] GT nz-1 THEN zind[1] = nz-1
  IF zind[0] GT zind[1] THEN BEGIN
    tmp = zind[0]
    zind[0] = zind[1]
    zind[1] = tmp
  ENDIF
  
  IF tind[0] LT 0 THEN tind[0] = 0
  IF tind[0] GT nt-1 THEN tind[0] = nt-1
  IF tind[1] LT 0 THEN tind[1] = 0
  IF tind[1] GT nt-1 THEN tind[1] = nt-1
  IF tind[0] GT tind[1] THEN BEGIN
    tmp = tind[0]
    tind[0] = tind[1]
    tind[1] = tmp
  ENDIF
  
  xsize = xind[1] - xind[0] + 1
  ysize = yind[1] - yind[0] + 1
  zsize = zind[1] - zind[0] + 1
  tsize = tind[1] - tind[0] + 1
  
  IF ndims EQ 4 THEN BEGIN
    ; print ranges

    IF quiet EQ 0 THEN BEGIN
        PRINT, "X indices: ", xind[0], " to ", xind[1]
        PRINT, "Y indices: ", yind[0], " to ", yind[1]
        PRINT, "Z indices: ", zind[0], " to ", zind[1]
        PRINT, "T indices: ", tind[0], " to ", tind[1]
    ENDIF
    
    data = FLTARR(xsize, ysize, zsize, tsize)
    
    FOR i=0, nfiles-1 DO BEGIN
      ; get X and Y processor indices
      pe_yind = FIX(i / NXPE)
      pe_xind = i MOD NXPE
      
      ; get local y range
      ymin = yind[0] - pe_yind*MYSUB + MYG
      ymax = yind[1] - pe_yind*MYSUB + MYG
      
      xmin = xind[0] - pe_xind*MXSUB
      xmax = xind[1] - pe_xind*MXSUB
      
      inrange = 1
      
      IF (ymin GE (MYSUB + MYG)) OR (ymax LT MYG) THEN BEGIN
        inrange = 0 ; out of Y range
      ENDIF
      
      IF ymin LT MYG THEN ymin = MYG
      IF ymax GE MYSUB + MYG THEN ymax = MYG + MYSUB - 1

      ; check lower x boundary
      IF PE_XIND EQ 0 THEN BEGIN
        ; keeping inner boundary
        IF xmax LT 0 THEN inrange = 0
        IF xmin LT 0 THEN xmin = 0
      ENDIF ELSE BEGIN
        IF xmax LT MXG THEN inrange = 0
        IF xmin LT MXG THEN xmin = MXG
      ENDELSE
      
      ; check upper boundary
      IF PE_XIND EQ (NXPE - 1) THEN BEGIN
        ; keeping outer boundary
        IF xmin GE (MXSUB+2*MXG) THEN inrange = 0
        IF xmax GE (MXSUB+2*MXG-1) THEN xmax = (MXSUB+2*MXG-1)
      ENDIF ELSE BEGIN
        IF xmin GE (MXSUB+MXG) THEN inrange = 0
        IF xmax GE (MXSUB+MXG) THEN xmax = (MXSUB+MXG-1)
      ENDELSE
      
      ; calculate global indices
      xgmin = xmin + PE_XIND * MXSUB
      xgmax = xmax + PE_XIND * MXSUB
      
      ygmin = ymin + PE_YIND*MYSUB - MYG
      ygmax = ymax + PE_YIND*MYSUB - MYG
      
      IF inrange THEN BEGIN
        filename = path+"/"+prefix+"."+STRTRIM(STRING(i),2)+"."+fext
        IF quiet EQ 0 THEN BEGIN
            PRINT, ""
            PRINT, "Reading from "+filename
            PRINT, "PE_X = " + STRTRIM(STRING(pe_xind),2) + $
              " PE_Y = " + STRTRIM(STRING(pe_yind),2)
            PRINT, "Local indices: ["+STRTRIM(STRING(xmin),2)+"-"+STRTRIM(STRING(xmax),2)+"][" + $
              STRTRIM(STRING(ymin),2)+"-"+STRTRIM(STRING(ymax),2)+"]"
            PRINT, "   =>  Global: ["+STRTRIM(STRING(xgmin),2)+"-"+STRTRIM(STRING(xgmax),2)+"][" + $
              STRTRIM(STRING(ygmin),2)+"-"+STRTRIM(STRING(ygmax),2)+"]"
        ENDIF ELSE IF quiet EQ 1 THEN BEGIN
            WRITEU, -1, 13, "Reading from "+filename+": ["+STRTRIM(STRING(xmin),2)+"-"+STRTRIM(STRING(xmax),2)+"][" + $
              STRTRIM(STRING(ymin),2)+"-"+STRTRIM(STRING(ymax),2)+"] -> ["+STRTRIM(STRING(xgmin),2)+"-"+STRTRIM(STRING(xgmax),2)+"][" + $
              STRTRIM(STRING(ygmin),2)+"-"+STRTRIM(STRING(ygmax),2)+"]"
        ENDIF
        
        
        ; set the indices
        inds = LONARR(8)
        inds[0] = LONG(tind[0])
        inds[1] = LONG(tind[1])
        inds[2] = LONG(xmin)
        inds[3] = LONG(xmax)
        inds[4] = LONG(ymin)
        inds[5] = LONG(ymax)
        inds[6] = LONG(zind[0])
        inds[7] = LONG(zind[1])
        
        ny_loc = ymax - ymin + 1
        nx_loc = xmax - xmin + 1
        
        handle = file_open(filename)
        d = file_read(handle, var, inds=inds)
        file_close, handle
        
        ; now copy into data array
        FOR jx=0, nx_loc-1 DO BEGIN
          FOR jy=0, ny_loc-1 DO BEGIN
            FOR t=0, tsize-1 DO BEGIN
              data[xgmin + jx - xind[0], ygmin + jy - yind[0], *, t] = d[t, jx, jy, *]
            ENDFOR
          ENDFOR
        ENDFOR
        
      ENDIF
    ENDFOR
  ENDIF ELSE IF ndims EQ 3 THEN BEGIN
    
    ; could be xyz or txy
    isxyz = 0
    IF dimsize[0] EQ nz THEN isxyz = 1

    ; print ranges
    IF quiet EQ 0 THEN BEGIN
      PRINT, "X indices: ", xind[0], " to ", xind[1]
      PRINT, "Y indices: ", yind[0], " to ", yind[1]
      IF isxyz THEN BEGIN
        PRINT, "Z indices: ", zind[0], " to ", zind[1]
      ENDIF ELSE BEGIN
        PRINT, "T indices: ", tind[0], " to ", tind[1]
      ENDELSE
    ENDIF
    
    IF isxyz THEN BEGIN
      data = FLTARR(xsize, ysize, zsize)
    ENDIF ELSE BEGIN
      data = FLTARR(xsize, ysize, tsize)
    ENDELSE
      
    FOR i=0, nfiles-1 DO BEGIN
      ; get X and Y processor indices
      pe_yind = FIX(i / NXPE)
      pe_xind = i MOD NXPE
      
      ; get local y range
      ymin = yind[0] - pe_yind*MYSUB + MYG
      ymax = yind[1] - pe_yind*MYSUB + MYG
        
      xmin = xind[0] - pe_xind*MXSUB
      xmax = xind[1] - pe_xind*MXSUB
      
      inrange = 1
      
      IF (ymin GE (MYSUB + MYG)) OR (ymax LT MYG) THEN BEGIN
        inrange = 0 ; out of Y range
      ENDIF
      
      IF ymin LT MYG THEN ymin = MYG
      IF ymax GE MYSUB + MYG THEN ymax = MYG + MYSUB - 1
      
      ; check lower x boundary
      IF PE_XIND EQ 0 THEN BEGIN
        ; keeping inner boundary
        IF xmax LT 0 THEN inrange = 0
        IF xmin LT 0 THEN xmin = 0
      ENDIF ELSE BEGIN
        IF xmax LT MXG THEN inrange = 0
        IF xmin LT MXG THEN xmin = MXG
      ENDELSE
      
      ; check upper boundary
      IF PE_XIND EQ (NXPE - 1) THEN BEGIN
        ; keeping outer boundary
        IF xmin GE (MXSUB+2*MXG) THEN inrange = 0
        IF xmax GE (MXSUB+2*MXG-1) THEN xmax = (MXSUB+2*MXG-1)
      ENDIF ELSE BEGIN
        IF xmin GE (MXSUB+MXG) THEN inrange = 0
        IF xmax GE (MXSUB+MXG) THEN xmax = (MXSUB+2*MXG-1)
      ENDELSE
      
      ; calculate global indices
      xgmin = xmin + PE_XIND * MXSUB
      xgmax = xmax + PE_XIND * MXSUB
      
      ygmin = ymin + PE_YIND*MYSUB - MYG
      ygmax = ymax + PE_YIND*MYSUB - MYG
      
      IF inrange THEN BEGIN
        
        filename = path+"/"+prefix+"."+STRTRIM(STRING(i),2)+"."+fext
        IF quiet EQ 0 THEN BEGIN
            PRINT, ""
            PRINT, "Reading from "+filename
            PRINT, "PE_X = " + STRTRIM(STRING(pe_xind),2) + $
              " PE_Y = " + STRTRIM(STRING(pe_yind),2)
            PRINT, "Local indices: ["+STRTRIM(STRING(xmin),2)+"-"+STRTRIM(STRING(xmax),2)+"][" + $
              STRTRIM(STRING(ymin),2)+"-"+STRTRIM(STRING(ymax),2)+"]"
            PRINT, "   =>  Global: ["+STRTRIM(STRING(xgmin),2)+"-"+STRTRIM(STRING(xgmax),2)+"][" + $
              STRTRIM(STRING(ygmin),2)+"-"+STRTRIM(STRING(ygmax),2)+"]"
        ENDIF  ELSE IF quiet EQ 1 THEN BEGIN
            WRITEU, -1, 13, "Reading from "+filename+": ["+STRTRIM(STRING(xmin),2)+"-"+STRTRIM(STRING(xmax),2)+"][" + $
              STRTRIM(STRING(ymin),2)+"-"+STRTRIM(STRING(ymax),2)+"] -> ["+STRTRIM(STRING(xgmin),2)+"-"+STRTRIM(STRING(xgmax),2)+"][" + $
              STRTRIM(STRING(ygmin),2)+"-"+STRTRIM(STRING(ygmax),2)+"]"
        ENDIF
        
        inds = LONARR(6)
        IF isxyz THEN BEGIN
          inds[0] = LONG(xmin)
          inds[1] = LONG(xmax)
          inds[2] = LONG(ymin)
          inds[3] = LONG(ymax)
          inds[4] = LONG(zind[0])
          inds[5] = LONG(zind[1])
        ENDIF ELSE BEGIN
          inds[0] = LONG(tind[0])
          inds[1] = LONG(tind[1])
          inds[2] = LONG(xmin)
          inds[3] = LONG(xmax)
          inds[4] = LONG(ymin)
          inds[5] = LONG(ymax)
        ENDELSE

        ny_loc = ymax - ymin + 1
        nx_loc = xmax - xmin + 1

        ; read the data
        handle = file_open(filename)
        d = file_read(handle, var, inds=inds)
        file_close, handle
        
        ; now copy into data array
        IF isxyz THEN BEGIN
          data[(xgmin - xind[0]):(xgmin - xind[0] + nx_loc - 1), (ygmin - yind[0]):(ygmin - yind[0] + ny_loc - 1), *] = d
        ENDIF ELSE BEGIN
          FOR jx=0, nx_loc-1 DO BEGIN
            FOR jy=0, ny_loc-1 DO BEGIN
              FOR t=0, tsize-1 DO BEGIN
                data[xgmin + jx - xind[0], ygmin + jy - yind[0], t] = d[t, jx, jy]
              ENDFOR
            ENDFOR
          ENDFOR
        ENDELSE
        
      ENDIF
    ENDFOR

  ENDIF ELSE IF ndims EQ 2 THEN BEGIN
    ; 2-D data, e.g. profiles. Split across processors the same as 3D
    
    IF NOT KEYWORD_SET(xind) THEN xind = [0, nx-1]
    IF NOT KEYWORD_SET(yind) THEN yind = [0, ny-1]
    
    ; check ranges
    
    IF xind[0] LT 0 THEN xind[0] = 0
    IF xind[0] GT nx-1 THEN xind[0] = nx-1
    IF xind[1] LT 0 THEN xind[1] = 0
    IF xind[1] GT nx-1 THEN xind[1] = nx-1
    IF xind[0] GT xind[1] THEN BEGIN
      tmp = xind[0]
      xind[0] = xind[1]
      xind[1] = tmp
    ENDIF
    
    IF yind[0] LT 0 THEN yind[0] = 0
    IF yind[0] GT ny-1 THEN yind[0] = ny-1
    IF yind[1] LT 0 THEN yind[1] = 0
    IF yind[1] GT ny-1 THEN yind[1] = ny-1
    IF yind[0] GT yind[1] THEN BEGIN
      tmp = yind[0]
      yind[0] = yind[1]
      yind[1] = tmp
    ENDIF
    
    ; print ranges

    IF quiet EQ 0 THEN BEGIN
        PRINT, "X indices: ", xind[0], " to ", xind[1]
        PRINT, "Y indices: ", yind[0], " to ", yind[1]
    ENDIF
    
    xsize = xind[1] - xind[0] + 1
    ysize = yind[1] - yind[0] + 1
    
    data = FLTARR(xsize, ysize)
    
    FOR i=0, nfiles-1 DO BEGIN
      ; get X and Y processor indices
      pe_yind = FIX(i / NXPE)
      pe_xind = i MOD NXPE
      
      ; get local y range
      ymin = yind[0] - pe_yind*MYSUB + MYG
      ymax = yind[1] - pe_yind*MYSUB + MYG
      
      xmin = xind[0] - pe_xind*MXSUB
      xmax = xind[1] - pe_xind*MXSUB
      
      inrange = 1
      
      IF (ymin GE (MYSUB + MYG)) OR (ymax LT MYG) THEN BEGIN
        inrange = 0 ; out of Y range
      ENDIF

      IF ymin LT MYG THEN ymin = MYG
      IF ymax GE MYSUB + MYG THEN ymax = MYG + MYSUB - 1
      
      ; check lower x boundary
      IF PE_XIND EQ 0 THEN BEGIN
        ; keeping inner boundary
        IF xmax LT 0 THEN inrange = 0
        IF xmin LT 0 THEN xmin = 0
      ENDIF ELSE BEGIN
        IF xmax LT MXG THEN inrange = 0
        IF xmin LT MXG THEN xmin = MXG
      ENDELSE
      

      ; check upper boundary
      IF PE_XIND EQ (NXPE - 1) THEN BEGIN
        ; keeping outer boundary
        IF xmin GE (MXSUB+2*MXG) THEN inrange = 0
        IF xmax GE (MXSUB+2*MXG-1) THEN xmax = (MXSUB+2*MXG-1)
      ENDIF ELSE BEGIN
        IF xmin GE (MXSUB+MXG) THEN inrange = 0
        IF xmax GE (MXSUB+MXG) THEN xmax = (MXSUB+2*MXG-1)
      ENDELSE
      
      ; calculate global indices
      xgmin = xmin + PE_XIND * MXSUB
      xgmax = xmax + PE_XIND * MXSUB
      
      ygmin = ymin + PE_YIND*MYSUB - MYG
      ygmax = ymax + PE_YIND*MYSUB - MYG
      
      IF inrange THEN BEGIN
        
        filename = path+"/"+prefix+"."+STRTRIM(STRING(i),2)+"."+fext
        IF quiet EQ 0 THEN BEGIN
            PRINT, ""
            PRINT, "Reading from "+filename
            PRINT, "PE_X = " + STRTRIM(STRING(pe_xind),2) + $
              " PE_Y = " + STRTRIM(STRING(pe_yind),2)
            PRINT, "Local indices: ["+STRTRIM(STRING(xmin),2)+"-"+STRTRIM(STRING(xmax),2)+"][" + $
              STRTRIM(STRING(ymin),2)+"-"+STRTRIM(STRING(ymax),2)+"]"
            PRINT, "   =>  Global: ["+STRTRIM(STRING(xgmin),2)+"-"+STRTRIM(STRING(xgmax),2)+"][" + $
              STRTRIM(STRING(ygmin),2)+"-"+STRTRIM(STRING(ygmax),2)+"]"
        ENDIF  ELSE IF quiet EQ 1 THEN BEGIN
            WRITEU, -1, 13, "Reading from "+filename+": ["+STRTRIM(STRING(xmin),2)+"-"+STRTRIM(STRING(xmax),2)+"][" + $
              STRTRIM(STRING(ymin),2)+"-"+STRTRIM(STRING(ymax),2)+"] -> ["+STRTRIM(STRING(xgmin),2)+"-"+STRTRIM(STRING(xgmax),2)+"][" + $
              STRTRIM(STRING(ygmin),2)+"-"+STRTRIM(STRING(ygmax),2)+"]"
        ENDIF

        inds = LONARR(4)
        inds[0] = LONG(xmin)
        inds[1] = LONG(xmax)
        inds[2] = LONG(ymin)
        inds[3] = LONG(ymax)

        ny_loc = ymax - ymin + 1
        nx_loc = xmax - xmin + 1

        ; read the data
        handle = file_open(filename)
        d = file_read(handle, var, inds=inds)
        file_close, handle
        
        ; now copy into data array
        FOR jx=0, nx_loc-1 DO BEGIN
          FOR jy=0, ny_loc-1 DO BEGIN
            data[xgmin + jx - xind[0], ygmin + jy - yind[0]] = d[jx, jy]
          ENDFOR
        ENDFOR
        
      ENDIF
    ENDFOR
    
  ENDIF ELSE BEGIN
    ; 0 or 1-D
    ; Just read the variable from main file
    IF quiet LT 2 THEN PRINT, "Reading from "+mainfile
    handle = file_open(mainfile)
    data = file_read(handle, var)
    file_close, handle
  ENDELSE

  IF quiet EQ 1 THEN PRINT, ""

  CATCH, /CANCEL

  RETURN, data
END
