; Produce an ERGOS-style plot of mode structure
; 

PRO ergos_plot, var3d, grid, period=period, mode=mode
  IF NOT KEYWORD_SET(period) THEN period = 1 ; default = full torus
  IF NOT KEYWORD_SET(mode) THEN mode = 1 ; Lowest harmonic
  s = SIZE(var3d, /DIMEN)
  IF N_ELEMENTS(s) NE 3 THEN BEGIN
     PRINT, "ERROR: var3D must be a 3D variable"
     RETURN
  ENDIF
  nx = s[0]
  ny = s[1]
  nz = s[2]
  
  IF NOT is_pow2(nz) THEN BEGIN
     IF is_pow2(nz-1) THEN BEGIN
        ; Remove the final z point
        nz = nz - 1
        data = var3d[*,*,0:(nz-1)]
     ENDIF ELSE BEGIN
        PRINT, "ERROR: var3D z dimension should be a power of 2"
        RETURN
     ENDELSE
  ENDIF ELSE data = var3d

  ; Take FFT of data in Z (toroidal angle), selecting
  ; the desired mode
  varfft = COMPLEXARR(nx, ny)
  
  FOR i=0, nx-1 DO BEGIN
     FOR j=0, ny-1 DO BEGIN
        varfft[i,j] = (FFT(data[i,j,*]))[mode]
     ENDFOR
  ENDFOR

  dz = 2.0*!PI / FLOAT(period*nz)
  
  ; GET THE TOROIDAL SHIFT
  tn = TAG_NAMES(grid)
  tn = STRUPCASE(tn)
  w = WHERE(tn EQ "QINTY", count)
  IF count GT 0 THEN BEGIN
     PRINT, "Using qinty as toroidal shift angle"
     zShift = grid.qinty
  ENDIF ELSE BEGIN
     w = WHERE(tn EQ "ZSHIFT", count)
     IF count GT 0 THEN BEGIN
        PRINT, "Using zShift as toroidal shift angle"
        zShift = grid.zShift
     ENDIF ELSE BEGIN
        PRINT, "ERROR: Can't find qinty or zShift variable"
        RETURN
      ENDELSE
  ENDELSE
  
  
END
