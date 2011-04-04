; Convert a BOUT-06 collect file to a set of BOUT++ restart files
; 
; 1. Use the BOUT-06 collect_fluc_txyz to collect 
;    the evolving variables
; 2. Run b06_to_bpp,  passing the collected file (BOUT_FLUC_xyzt.pdb)
;    as cfile argument (either a string or a structure)
; 3. Enter the names of variables to copy across, "none" to skip
;
;
; Doesn't handle Y boundaries (target plates)
;

PRO b06_to_bpp, cfile, tind=tind, NPES=NPES, NXPE=NXPE, path=path
  IF NOT KEYWORD_SET(NPES) THEN NPES = 1
  IF NOT KEYWORD_SET(NXPE) THEN NXPE = 1
  
  IF NOT KEYWORD_SET(path) THEN path = "."
  
  IF (NPES MOD NXPE) NE 0 THEN BEGIN
    PRINT, "ERROR: NPES ("+STR(NPES)+") not a multiple of NXPE ("+STR(NXPE)+")"
    RETURN
  ENDIF

  NYPE = NPES / NXPE

  NPES = LONG(NPES)
  NXPE = LONG(NXPE)

  IF SIZE(cfile, /type) EQ 7 THEN BEGIN
    ; A file-name, so import
    data = pd_import(cfile)
  ENDIF ELSE data = cfile
  
  s = SIZE(data.ni_xyzt, /dimensions)
  nx = s[0]
  ny = s[1]
  nz = s[2]
  nt = s[3]
  
  IF NOT KEYWORD_SET(tind) THEN tind = nt-1
  IF tind GT (nt-1) THEN tind = nt-1
  IF tind LT 0 THEN tind = 0
  
  IF (ny MOD nype) NE 0 THEN BEGIN
    PRINT, "ERROR: ny ("+str(ny)+") doesn't divide equally between nype ("+str(nype)+")"
    RETURN
  ENDIF
  
  IF ((nx - 4) MOD nxpe) NE 0 THEN BEGIN
    PRINT, "ERROR: nx-4 ("+str(nx-4)+") doesn't divide equally between nxpe ("+str(nxpe)+")"
    RETURN
  ENDIF
  
  mysub = FIX(ny / nype)
  mxsub = FIX((nx-4) / nxpe)
  
  IF mysub LT 2 THEN BEGIN
    PRINT, "ERROR: MYSUB too small ("+STR(mysub)+")"
    RETURN
  ENDIF

  IF mxsub LT 2 THEN BEGIN
    PRINT, "ERROR: MXSUB too small ("+STR(mxsub)+")"
    RETURN
  ENDIF

  ; Get list of variables
  tn = TAG_NAMES(data)
  ; Need to ask the user to translate the variable names
  PRINT, "Translate BOUT-06 variable names to BOUT++"
  PRINT, "Enter 'none' to skip"
  
  ncopy = 0
  FOR i=0, N_ELEMENTS(tn)-1 DO BEGIN
    IF SIZE(data.(i), /n_dim) EQ 4 THEN BEGIN
      ; an evolving variable
      
      bvar = ""
      READ, bvar, prompt="Variable "+tn[i]+" -> "
      IF STRCMP(bvar, "none", /fold_case) EQ 0 THEN BEGIN
        IF ncopy EQ 0 THEN BEGIN
          copyfrom = [tn[i]]
          copyto = [bvar]
          copyind = [i]
        ENDIF ELSE BEGIN
          copyfrom = [copyfrom, tn[i]]
          copyto = [copyto, bvar]
          copyind = [copyind, i]
        ENDELSE
        ncopy = ncopy + 1
      ENDIF
      
    ENDIF
  ENDFOR
  
  hist_hi = LONG(data.hist_hi)
  tt = DOUBLE(0.)
  
  xmin = 0
  xmax = nx-1
  FOR yp = 0, nype-1 DO BEGIN
    ymin = yp*mysub
    ymax = ymin + mysub - 1
    
    FOR xp=0, nxpe-1 DO BEGIN
      xmin = 2+xp*mxsub
      xmax = 1+mxsub + xp*mxsub

      xi0 = 2
      xi1 = 1+mxsub

      IF xp EQ 0 THEN BEGIN 
        xmin = 0
        xi0 = 0
      ENDIF
      IF xp EQ nxpe-1 THEN BEGIN
        xmax = nx-1
        xi1 = 3+mxsub
      ENDIF
      
      pe = NXPE*yp + xp
    
      name = path+"/BOUT.restart."+STR(pe)+".nc"
      fp = file_open(name, /create)
      
      status = file_write(fp, "hist_hi", hist_hi)
      status = file_write(fp, "tt", tt)
      status = file_write(fp, "NPES", NPES)
      status = file_write(fp, "NXPE", NXPE)
    
      ; Copy across the data
      FOR i=0, ncopy-1 DO BEGIN
        new_nz = nz
        IF is_pow2(nz) THEN BEGIN
          ; Need to add an additional point
          new_nz = nz + 1
        ENDIF
      
        var = DBLARR(mxsub+4, mysub+4, new_nz)
        var[xi0:xi1,2:(mysub+1),0:(nz-1)] = REFORM((data.(copyind[i]))[xmin:xmax, ymin:ymax, *,tind])
        
        IF is_pow2(nz) THEN var[*,*,nz] = var[*,*,0]
        
        status = file_write(fp, copyto[i], var)
      ENDFOR
    
      file_close, fp
    ENDFOR
  ENDFOR
END
