; Produce an ERGOS-style plot of mode structure
; 

PRO ergos_plot, var3d, grid, period=period, mode=mode, noshift=noshift, title=title, _extra=_extra
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
  
  tn = TAG_NAMES(grid)
  tn = STRUPCASE(tn)
  
  ncore = MIN([grid.ixseps1, grid.ixseps2])
  IF nx GT ncore THEN BEGIN
    PRINT, "Only computing modes in the core: "+STR(ncore)+" of "+STR(nx)
    nx = ncore
  ENDIF

  ; Get number of poloidal points in core
  ny = 0
  status = gen_surface(mesh=grid) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=periodic)
    
    IF periodic THEN BEGIN
      ; Found a periodic region. Assume in the core
      ny = N_ELEMENTS(yi)
      BREAK
    ENDIF
  ENDREP UNTIL last
  IF ny EQ 0 THEN BEGIN
    PRINT, "ERROR: No periodic (core) region found!"
    RETURN
  ENDIF
  
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
  

  IF NOT KEYWORD_SET(noshift) THEN BEGIN
    ; Shifting - need the angle
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
  ENDIF

  status = gen_surface(mesh=grid) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=periodic)
    
    IF periodic THEN BEGIN
      IF N_ELEMENTS(yi) NE ny THEN BEGIN
        PRINT, "ERROR: Number of poloidal points in core varies!"
        RETURN
      ENDIF
      
      FOR j=0, ny-1 DO BEGIN
        varfft[xi,j] = (FFT(REFORM(data[xi,yi[j],*])))[mode]
      ENDFOR
      
      IF NOT KEYWORD_SET(noshift) THEN BEGIN
        ; Apply toroidal shift
        kwave = mode * period
        FOR j=0, ny-1 DO BEGIN
          varfft[xi,j] = varfft[xi,j] * COMPLEX(COS(kwave*zshift[xi,yi[j]]), $
                                                -SIN(kwave*zshift[xi,yi[j]]))
        ENDFOR

      ENDIF
    ENDIF
  ENDREP UNTIL last

  zlength = 2.0*!PI / FLOAT(period)
  dz = zlength / FLOAT(nz)
  
  ; Try to get the poloidal angle from the grid file
  w = WHERE(tn EQ "POL_ANGLE", count)
  IF count GT 0 THEN BEGIN
    PRINT, "Using pol_angle from grid file"
    theta2D = grid.pol_angle
  ENDIF ELSE BEGIN
    ; calculate poloidal angle
    
    PRINT, "Sorry: Poloidal angle calculation not implemented yet"
    RETURN
  ENDELSE
  
  theta = FLTARR(nx, ny)
  status = gen_surface(mesh=grid) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=periodic)
    
    IF periodic THEN BEGIN
      theta[xi,*] = theta2D[xi,yi]
    ENDIF
  ENDREP UNTIL last
  
  ; Now take FFT in poloidal angle
  
  nm = FIX(ny/2)-1 ; Number of m frequencies
  
  yfft = COMPLEXARR(nx, 2*nm)
  
  FOR i=0, nx-1 DO BEGIN
    yfft[i,*] = fft_irreg(theta[i,*], varfft[i,*], nf=nm)
  ENDFOR
  
  ; Now re-arrange and take absolute value
  
  IF ny MOD 2 EQ 0 THEN BEGIN
    ; Even number of points
    nm = FIX(ny/2)-1
    
    mvals = FINDGEN(2*nm + 1)-nm
    result = FLTARR(2*nm+1, nx)
    
    FOR i=0, nx-1 DO BEGIN
      result[nm, i] = ABS(yfft[i,0]) ; DC
      FOR m=0, nm-1 DO BEGIN
        ; First negative frequencies
        result[m, i] = ABS(yfft[i, nm+m])
        ; Positive frequencies
        result[nm+1+m, i] = ABS(yfft[i, 1+m])
      ENDFOR
    ENDFOR
  ENDIF ELSE BEGIN
    PRINT, "Sorry: NY odd not implemented yet"
  ENDELSE
  
  ; GET THE TOROIDAL SHIFT
  w = WHERE(tn EQ "QSAFE", count)
  IF count GT 0 THEN BEGIN
    PRINT, "Using qsafe for q profile"
    qsafe = -grid.qinty
  ENDIF ELSE BEGIN
    w = WHERE(tn EQ "SHIFTANGLE", count)
    IF count GT 0 THEN BEGIN
      PRINT, "Using shiftangle for q profile"
      qsafe = -grid.shiftangle
    ENDIF ELSE BEGIN
      PRINT, "WARNING: Can't find q profile"
      qsafe = FLTARR(nx)
    ENDELSE
  ENDELSE
  
  psin = (REFORM(grid.psixy[0:(nx-1),grid.jyseps1_1+1]) - grid.psi_axis) / (grid.psi_bndry - grid.psi_axis)

  qsafe = qsafe[0:(nx-1)] / (2.*!PI)

  tit = "Mode spectrum, n="+STRTRIM(STRING(mode),2)
  IF KEYWORD_SET(title) THEN tit = title + " " + tit
  ;LOADCT, 0
  ;DEVICE, decomposed=0
  CONTOUR2, result, mvals, SQRT(psin), nlev=50, /fill, $
           title=tit, $
           xtitle="Poloidal mode number m", $
           ytitle="SQRT(Poloidal flux psi)", $
           _extra=_extra ;, /cent, /redblue
  
  OPLOT, qsafe*mode, SQRT(psin), color=1, thick=2
  
  ; Mark resonant surfaces
  
  maxm = FIX(MAX(qsafe)) * mode
  qreson = (FINDGEN(maxm)+1) / FLOAT(mode)
  
  ; get x location for each of these resonances
  qpsi = INTERPOL(psin, qsafe, qreson)
  
  oplot, findgen(maxm)+1., SQRT(qpsi), psym=4, color=1
END
