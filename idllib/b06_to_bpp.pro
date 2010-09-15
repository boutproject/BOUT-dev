; Convert a BOUT-06 collect file to a set of BOUT++ restart files
; 
; 1. Use the BOUT-06 collect_fluc_txyz to collect 
;    the evolving variables
; 2. Run b06_to_bpp,  passing the collected file (BOUT_FLUC_xyzt.pdb)
;    as cfile argument (either a string or a structure)
; 3. Enter the names of variables to copy across, "none" to skip
;
; 
; Currently won't decompose in X, only in Y
;

PRO b06_to_bpp, cfile, tind=tind, NPES=NPES, path=path
  
  NXPE = 1
  IF NOT KEYWORD_SET(NPES) THEN NPES = 1
  
  IF NOT KEYWORD_SET(path) THEN path = "."
  
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
  
  IF ny MOD nype NE 0 THEN BEGIN
    PRINT, "ERROR: ny ("+str(ny)+") doesn't divide equally between nype ("+str(nype)+")"
    RETURN
  ENDIF
  
  mysub = FIX(ny / nype)

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
    
    pe = NXPE*yp
    
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
      
      var = DBLARR(nx, ny+4, new_nz)
      var[*,2:(ny+1),0:(nz-1)] = REFORM((data.(copyind[i]))[xmin:xmax, ymin:ymax, *,tind])
      
      IF is_pow2(nz) THEN var[*,*,nz] = var[*,*,0]
      
      status = file_write(fp, copyto[i], var)
    ENDFOR
    
    file_close, fp
  ENDFOR
END
