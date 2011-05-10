; Non-orthogonal tokamak grid generator
; =====================================
;
; First creates flux-surfaces, then creates a mesh which conforms to 
; boundary shapes.
;
; Useage:
; ------
; 
; FUNCTION create_nonorthogonal, F, R, Z, settings, critical=critical,
;                                boundary=boundary, fpsi=fpsi
;
; F - psi(nr, nz) 2D array
; R - R(nr)  1D array
; Z - Z(nz)  1D array


FUNCTION create_nonorthogonal, F, R, Z, in_settings, critical=critical, $
                               boundary=boundary, fpsi=fpsi, strict=strict, $
                               nrad_flexible=nrad_flexible, rad_peaking=rad_peaking, $
                               single_rad_grid=single_rad_grid
  
  IF SIZE(nrad_flexible, /TYPE) EQ 0 THEN nrad_flexible = 0
  IF NOT KEYWORD_SET(rad_peaking) THEN rad_peaking = 2.0
  rad_peaking = 1.0 / rad_peaking
  
  IF N_PARAMS() LT 3 THEN BEGIN
     PRINT, "ERROR: Need at least a 2D array of psi values, R and Z arrays"
     RETURN, {error:1}
  ENDIF ELSE IF N_PARAMS() LT 4 THEN BEGIN
                                ; Settings omitted. Set defaults
     PRINT, "Settings not given -> using default values"
     settings = {psi_inner:0.9, $
                 psi_outer:1.1, $
                 nrad:36, $
                 npol:64, $
                 rad_peaking:0.0, $
                 pol_peaking:0.0}
  ENDIF ELSE BEGIN
     PRINT, "Checking settings"
     settings = in_settings     ; So the input isn't changed
     str_check_present, settings, 'psi_inner', 0.9
     str_check_present, settings, 'psi_outer', 1.1
     str_check_present, settings, 'nrad', 36
     str_check_present, settings, 'npol', 64
     str_check_present, settings, 'rad_peaking', 0.0
     str_check_present, settings, 'pol_peaking', 0.0
  ENDELSE

  s = SIZE(F, /DIMENSION)
  IF N_ELEMENTS(s) NE 2 THEN BEGIN
     PRINT, "ERROR: First argument must be 2D array of psi values"
     RETURN, {error:1}
  ENDIF
  nx = s[0]
  ny = s[1]

  s = SIZE(R, /DIMENSION)
  IF (N_ELEMENTS(s) NE 1) OR (s[0] NE nx) THEN BEGIN
     PRINT, "ERROR: Second argument must be 1D array of major radii"
     RETURN, {error:1}
  ENDIF
  
  s = SIZE(Z, /DIMENSION)
  IF (N_ELEMENTS(s) NE 1) OR (s[0] NE ny) THEN BEGIN
     PRINT, "ERROR: Second argument must be 1D array of heights"
     RETURN, {error:1}
  ENDIF

  IF KEYWORD_SET(boundary) THEN BEGIN
     s = SIZE(boundary, /DIMENSION)
     IF (N_ELEMENTS(s) NE 2) OR (s[0] NE 2) THEN BEGIN
        PRINT, "WARNING: boundary must be a 2D array: [2, n]. Ignoring"
        boundary = 0
     ENDIF ELSE BEGIN
        ; Calculate indices
        bndryi = boundary
        bndryi[0,*] = INTERPOL(FINDGEN(nx), R, REFORM(boundary[0,*]))
        bndryi[1,*] = INTERPOL(FINDGEN(ny), Z, REFORM(boundary[1,*]))
     ENDELSE
  ENDIF
  
  ;;;;;;;;;;;;;;; Calculate DCT ;;;;;;;;;;;;;;
  
  PRINT, "Calculating DCT..."
  DCT2Dslow, F, dctF
  PRINT, "Finished DCT"

  ;;;;;;;;;;;;;;;; First plot ;;;;;;;;;;;;;;;;

  nlev = 100
  minf = MIN(f)
  maxf = MAX(f)
  levels = findgen(nlev)*(maxf-minf)/FLOAT(nlev-1) + minf

  safe_colors, /first
  CONTOUR, F, R, Z, levels=levels, color=1, /iso, xstyl=1, ysty=1
  
  IF KEYWORD_SET(boundary) THEN BEGIN
    OPLOT, [REFORM(boundary[0,*]), boundary[0,0]], [REFORM(boundary[1,*]), boundary[1,0]], $
           thick=2,color=2
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Analyse the equilibrium to find O- and X-points
  
  IF NOT KEYWORD_SET(critical) THEN critical = analyse_equil(F, R, Z)
  critical = restrict_psi_range(critical, MAX(settings.psi_outer))

  ; Check that the critical points are inside the boundary
  IF KEYWORD_SET(boundary) THEN critical = critical_bndry(critical, bndryi)

  n_opoint = critical.n_opoint
  n_xpoint = critical.n_xpoint
  primary_opt = critical.primary_opt
  inner_sep   = critical.inner_sep
  opt_ri = critical.opt_ri
  opt_zi = critical.opt_zi
  opt_f  = critical.opt_f
  xpt_ri = critical.xpt_ri
  xpt_zi = critical.xpt_zi
  xpt_f  = critical.xpt_f

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  PRINT, "Refining x-point locations..."
  FOR i=0, n_xpoint-1 DO BEGIN
    PRINT, "  "+STR(i)+": "+STR([xpt_ri[i], xpt_zi[i]])+" "+STR(xpt_f[i])
    refine_xpoint, dctF, xpt_ri[i], xpt_zi[i], ri, zi
    a = local_gradient(dctF, ri, zi)
    PRINT, "   -> "+STR([ri, zi])+" "+STR(a.f)
    xpt_ri[i] = ri
    xpt_zi[i] = zi
    xpt_f[i] = a.f
  ENDFOR

  ; Overplot the separatrices, O-points
  oplot_critical, F, R, Z, critical

  ; Psi normalisation factors
  faxis = critical.opt_f[critical.primary_opt]
  fnorm = critical.xpt_f[critical.inner_sep] - critical.opt_f[critical.primary_opt]

  ; From normalised psi, get range of f
  f_inner = faxis + MIN(settings.psi_inner)*fnorm
  f_outer = faxis + MAX(settings.psi_outer)*fnorm

  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Get flux surface locations in psi
  
  ; Check the number of x-points
  IF critical.n_xpoint EQ 0 THEN BEGIN
    

  ENDIF ELSE BEGIN
    ; Normalised psi value of each separatrix
    xpt_psi = (critical.xpt_f - faxis) / fnorm
    
    si = SORT(xpt_psi) ; Sort separatrices from inside out
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Choose primary x-point.
    ; Determines order of settings arrays
    primary_xpt = si[0]
    
    PRINT, "Primary X-point is number "+STR(primary_xpt)
    PRINT, "   at R = "+STR(INTERPOLATE(R, critical.xpt_ri[primary_xpt])) $
      +" Z = "+STR(INTERPOLATE(Z, critical.xpt_zi[primary_xpt]))
    
    nrad = settings.nrad
    nnrad = N_ELEMENTS(nrad)
    
    IF nnrad NE (critical.n_xpoint + 1) THEN BEGIN
      ; Only given total number of radial points. Decide
      ; distribution automatically
      nrad_flexible = 1 ; Flag used later if nrad not fixed
      
      IF nnrad GT 1 THEN BEGIN
        PRINT, "WARNING: nrad has wrong number of elements ("+STR(nnrad)+")"
        PRINT, "         Should have 1 or "+STR(critical.n_xpoint + 1)+" elements"
      ENDIF
      
      PRINT, "Distributing radial points automatically"
      
      n = TOTAL(nrad,/int)
      fac = 2.*(xpt_f[inner_sep] - f_inner)/(1.+rad_peaking)
      FOR i=1, critical.n_xpoint-1 DO fac = fac + (xpt_f[si[i]] - xpt_f[si[i-1]])/rad_peaking
      fac = fac + 2.*(f_outer - xpt_f[si[critical.n_xpoint-1]])/(1.+rad_peaking)
      dx0 = fac / FLOAT(n)  ; Inner grid spacing
      
      ; Calculate number of grid points
      nrad = LONARR(critical.n_xpoint + 1)
      nrad[0] = FIX( 2.*(xpt_f[inner_sep] - f_inner) / ( (1.+rad_peaking)*dx0 ) + 0.5)
      FOR i=1, critical.n_xpoint-1 DO nrad[i] = FIX((xpt_f[si[i]] - xpt_f[si[i-1]])/(rad_peaking*dx0)-0.5)
      nrad[critical.n_xpoint] = n - TOTAL(nrad,/int)
    ENDIF
  ENDELSE
  
  IF KEYWORD_SET(single_rad_grid) THEN BEGIN
    ; Just produce one grid
    nrad_tot = TOTAL(nrad, /int)
    fvals = radial_grid(nrad_tot, f_inner, f_outer, 1, 1, xpt_f, rad_peaking)
    psi_vals = (fvals - faxis) / fnorm ; Normalised psi
    
    ; Find out how many grid points were used
    nrad = LONARR(critical.n_xpoint + 1)
    tot = 0
    FOR i=0, critical.n_xpoint-1 DO BEGIN
      w = WHERE(psi_vals LT xpt_psi[si[i]], count)
      nrad[i] = count - tot
      tot = tot + nrad[i]
    ENDFOR
    nrad[critical.n_xpoint] = nrad_tot - tot
    
  ENDIF ELSE BEGIN 
    IF critical.n_xpoint GT 1 THEN BEGIN
      ; Between separatrices
      fvals = radial_grid(nrad[1], xpt_f[inner_sep], xpt_f[si[1]], 0, 0, xpt_f, rad_peaking)
      
      FOR i=2, critical.n_xpoint-1 DO fvals = [fvals, radial_grid(nrad[i], xpt_f[si[i-1]], xpt_f[si[i]], 0, 0, xpt_f, rad_peaking)]
      ; Core
      fvals = [radial_grid(nrad[0], f_inner, 2.*xpt_f[inner_sep]-fvals[0], $
                           1, 1, xpt_f, rad_peaking, $
                           out_dp=2.*(fvals[0]-xpt_f[inner_sep]), $
                           in_dp=2.*(fvals[0]-xpt_f[inner_sep])/rad_peaking), fvals]
    ENDIF ELSE BEGIN
      ; Only a single separatrix
      dp0 = (xpt_f[inner_sep] - f_inner)*2./ (FLOAT(nrad[0])*(1. + rad_peaking))
      fvals = radial_grid(nrad[0], f_inner, xpt_f[inner_sep], $
                          1, 0, xpt_f, rad_peaking, $
                          out_dp=rad_peaking*dp0, $
                          in_dp=dp0)
    ENDELSE

    ; SOL
    n = N_ELEMENTS(fvals)
    dpsi = 2.*(xpt_f[si[critical.n_xpoint-1]] - fvals[n-1])
    fvals = [fvals, radial_grid(nrad[critical.n_xpoint], $
                                fvals[n-1]+dpsi, f_outer, 1, 1, xpt_f, rad_peaking, $
                                in_dp=dpsi, out_dp=dpsi/rad_peaking)]  
  ENDELSE
  
  psi_vals = (fvals - faxis) / fnorm ; Normalised psi
  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Get the contours
  ; find where they intersect the midplane and boundary
  
  midplane = [[opt_ri[primary_opt], opt_zi[primary_opt]], $
              [nx-1, opt_zi[primary_opt]]]
  
  
  FOR i=0, N_ELEMENTS(psi_vals)-1 DO BEGIN
    ; Get contour lines 
    contour_lines, F, findgen(nx), findgen(ny), levels=[fvals[i]], $
      path_info=info, path_xy=xy
    
    ; Find a line which crosses the midplane
    gotline = -1 ; Indicate no line found
    FOR j=0, N_ELEMENTS(info)-1 DO BEGIN

      line_ri = REFORM(xy[0,info[j].offset:(info[j].offset+info[j].n-1)])
      line_zi = REFORM(xy[1,info[j].offset:(info[j].offset+info[j].n-1)])

      cpos = line_crossings(line_ri, line_zi, info[j].type, $
                            midplane[0,*], midplane[1,*], 0, ncross=ncross, inds1=start_ind)
      IF ncross GT 0 THEN BEGIN
        gotline = j
        start_ind = start_ind[0]
        BREAK
      ENDIF
    ENDFOR
    
    IF gotline LT 0 THEN BEGIN
      PRINT, "WARNING: No contours at psi = ", psi_vals[i]
      BREAK
    ENDIF
    
    IF KEYWORD_SET(boundary) THEN BEGIN
      ; Find where the line crosses the boundary
      
      bpos = line_crossings(line_ri, line_zi, info[j].type, $
                            bndryi[0,*], bndryi[1,*], 1, ncross=ncross, inds1=b_ind)
      
      IF ncross GT 0 THEN BEGIN
        OPLOT, interpolate(R, bpos[0,*]), interpolate(z, bpos[1,*]), psym=1
        
        si = 0 ; Start index of the line
        w = WHERE(b_ind LT start_ind, count)
        IF count GT 0 THEN si = MAX(b_ind[w]) ; last boundary crossing before the midplane
        
        ei = N_ELEMENTS(line_ri)-1 ; end index
        w = WHERE(b_ind GT start_ind, count)
        IF count GT 0 THEN ei = MIN(b_ind[w]) ; First boundary crossing after midplane
        
        sii = CEIL(si)
        eii = FLOOR(ei)

        line_ri = [INTERPOLATE(line_ri, si), line_ri[sii:eii], INTERPOLATE(line_ri, ei)]
        line_zi = [INTERPOLATE(line_zi, si), line_zi[sii:eii], INTERPOLATE(line_zi, ei)]
        
        oplot, INTERPOLATE(R, line_ri), INTERPOLATE(Z, line_zi), color=4
      ENDIF ELSE BEGIN
        oplot, INTERPOLATE(R, line_ri), INTERPOLATE(Z, line_zi), color=2
      ENDELSE
      
      
    ENDIF
  ENDFOR
  
  FOR i=0, critical.n_xpoint-1 DO BEGIN
    PRINT, "Finding theta location of x-point "+STR(i)
    
    ; Get the separatrices
    legsep = leg_separatrix(dctF, R, Z, xpt_ri[i], xpt_zi[i], $
                            opt_ri[primary_opt], opt_zi[primary_opt], boundary=bndryi)
  ENDFOR

  STOP
  
END

