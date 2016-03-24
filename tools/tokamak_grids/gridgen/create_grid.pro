; Tokamak grid generator
; ======================
; 
; Generates a flux-surface aligned grid from
; an R-Z mesh of psi values. 
;
; Features:
; --------
;  o An arbitrary number of X-points
;  o Automatic default settings when not
;    supplied. 
;
; Author: Ben Dudson, University of York, Nov 2009
; 
;
; NOTE: Throughout, "F" means un-normalised psi,
;       and "psi" means normalised psi
;
; Useage:
; ------
; 
; FUNCTION create_grid, F, R, Z, settings, critical=critical,
;                                boundary=boundary, fpsi=fpsi
;
; F - psi(nr, nz) 2D array
; R - R(nr)  1D array
; Z - Z(nz)  1D array
; 
; Settings is a structure containing:
;   
;   psi_inner, psi_outer  - Range of normalised psi
;   nrad                  - Number of radial grid points
;                           Scalar -> Total number. Distributed
;                                     automatically
;                           Array -> Specified for each section
;   rad_peaking           - Radial separatrix peaking factor
;                           Not supplied -> 0
;                           scalar -> same for all regions
;                           array -> Different for each region
;   npol                  - Number of poloidal points.
;                           Scalar -> Total number. Distributed
;                                     automatically
;                           Array -> Specified for each section
;   pol_peaking           - Poloidal peaking factor
;                           Not supplied -> 0
;                           Scalar -> Same for all poloidal sections
;                           Array -> Specified for each section
; 
; A minimal settings structure must contain 4 scalars:
;    psi_inner, psi_outer, nrad and npol.
;
; Critical is a structure containing critical points:
;  
;   n_opoint, n_xpoint   - Number of O- and X-points
;   primary_opt          - Index of plasma centre O-point
;   inner_sep            - X-point index of inner separatrix
;   opt_ri, opt_zi       - R and Z indices for each O-point
;   opt_f                - Psi value at each O-point
;   xpt_ri, xpt_zi       - R and Z indices for each X-point
;   xpt_f                - Psi value of each X-point
; 
; If critical is not supplied, it is calculated
; 
; Boundary is a 2D array of indices
;   boundary[0, *]       - R values
;   boundary[1, *]       - Z values
;
; fpsi is an (optional) current function as a 2D array
;   fpsi[0,*]            - Psi values
;   fpsi[1,*]            - f values
;
; Return structure contains:
; 
;   error                       - Non-zero if an error occurred
;   psi_inner, psi_outer        - Normalised ranges of psi used
;   nrad[d], npol[d]            - Number of grid points used in each domain
;   yup_xsplit[d]               - Upper Y edge of domain. x < xsplit -> inner
;   yup_xin[d], yup_xout[d]     - Inner and outer domain connections. -1 = none
;   ydown_xsplit[d]             - Lower Y edge of domain. x < xsplit -> inner
;   ydown_xin[d], ydown_xout[d] - Inner and outer domain connections
; 
;   Rixy[x,y], Zixy[x,y]   - 2D arrays of indices into R and Z array
;   Rxy[x,y], Zxy[x,y]     - 2D arrays of (rad,pol) grid-point locations
;   psixy[x,y]             - 2D array of normalised psi at each point
;   faxis, fnorm           - Psi normalisation factors
;   settings               - Structure containing final settings used
;   critical               - Structure describing O- and X-points
;                            (from analyse_equil.pro)
;

PRO swap, a, b
  v = a
  a = b
  b = v
END

FUNCTION remove_ind, arr, ind
  n = N_ELEMENTS(arr)
  IF ind EQ 0 THEN BEGIN
    RETURN, arr[1:*]
  ENDIF ELSE IF ind EQ n-1 THEN BEGIN
    RETURN, arr[0:(ind-1)]
  ENDIF
  
  RETURN, [arr[0:(ind-1)], arr[(ind+1):*]]
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
; Poloidal grid
;
; Divide up a poloidal arc
;

FUNCTION poloidal_grid, interp_data, R, Z, ri, zi, n, fpsi=fpsi, parweight=parweight, $
                        ydown_dist=ydown_dist, yup_dist=yup_dist, $
                        ydown_space=ydown_space, yup_space=yup_space
  
  IF NOT KEYWORD_SET(parweight) THEN parweight = 0.0  ; Default is poloidal distance

  np = N_ELEMENTS(ri)

  IF 0 THEN BEGIN
    ; Calculate poloidal distance along starting line
    drdi = DERIV(INTERPOLATE(R, ri))
    dzdi = DERIV(INTERPOLATE(Z, zi))
    
    dldi = SQRT(drdi^2 + dzdi^2)
    poldist = int_func(findgen(np), dldi) ; Poloidal distance along line
  ENDIF ELSE BEGIN
    rpos = INTERPOLATE(R, ri)
    zpos = INTERPOLATE(Z, zi)
    dd = SQRT((zpos[1:*] - zpos[0:(np-2)])^2 + (rpos[1:*] - rpos[0:(np-2)])^2)
    dd = [dd, SQRT((zpos[0] - zpos[np-1])^2 + (rpos[0] - rpos[np-1])^2)]
    poldist = FLTARR(np)
    FOR i=1,np-1 DO poldist[i] = poldist[i-1] + dd[i-1]
  ENDELSE

  IF SIZE(fpsi, /n_dim) EQ 2 THEN BEGIN
    ; Parallel distance along line
    ; Need poloidal and toroidal field
    ni = N_ELEMENTS(ri)
    bp = FLTARR(ni)
    bt = FLTARR(ni)
    m = interp_data.method
    interp_data.method = 2
    FOR i=0, ni-1 DO BEGIN
      local_gradient, interp_data, ri[i], zi[i], status=status, $
        f=f, dfdr=dfdr, dfdz=dfdz
      ; dfd* are derivatives wrt the indices. Need to multiply by dr/di etc
      dfdr /= INTERPOLATE(DERIV(R),ri[i])
      dfdz /= INTERPOLATE(DERIV(Z),zi[i])
      
      IF i EQ 0 THEN BEGIN
        btr = INTERPOL(REFORM(fpsi[1,*]), REFORM(fpsi[0,*]), f)
      ENDIF
      
      bp[i] = SQRT(dfdr^2 + dfdz^2) / INTERPOLATE(R, ri[i])
      bt[i] = ABS( btr / INTERPOLATE(R, ri[i]))
    ENDFOR
    interp_data.method = m
    b = SQRT(bt^2 + bp^2)
    ddpar = dd * b / bp
    pardist = FLTARR(np)
    FOR i=1,np-1 DO BEGIN
      ip = (i + 1) MOD np
      pardist[i] = pardist[i-1] + ddpar[i] ;0.5*(ddpar[i-1] + ddpar[ip])
    ENDFOR
  ENDIF ELSE pardist = poldist ; Just use the same poloidal distance

  PRINT, "PARWEIGHT: ", parweight

  dist = parweight*pardist + (1. - parweight)*poldist

  ; Divide up distance. No points at the end (could be x-point)
  IF n GE 2 THEN BEGIN
    IF SIZE(ydown_dist, /TYPE) EQ 0 THEN ydown_dist = dist[np-1]* 0.5 / FLOAT(n)
    IF SIZE(yup_dist, /TYPE) EQ 0 THEN yup_dist = dist[np-1] * 0.5 / FLOAT(n)

    IF SIZE(ydown_space, /TYPE) EQ 0 THEN ydown_space = ydown_dist
    IF SIZE(yup_space, /TYPE) EQ 0 THEN yup_space = ydown_dist
    ;dloc = (dist[np-1] - ydown_dist - yup_dist) * FINDGEN(n)/FLOAT(n-1) + ydown_dist  ; Distance locations
    
    
    fn = FLOAT(n-1)
    d = (dist[np-1] - ydown_dist - yup_dist) ; Distance between first and last
    i = FINDGEN(n)
    
    yd = ydown_space < 0.5*d/fn
    yu = yup_space < 0.5*d/fn
    ; Fit to ai + bi^2 + c[i-sin(2pi*i/(n-1))*(n-1)/(2pi)]
    a = yd*2.
    b = (2.*yu - a) / fn
    c = d/fn - a - 0.5*b*fn
    dloc = ydown_dist + a*i + 0.5*b*i^2 + c*[i - SIN(2.*!PI*i / fn)*fn/(2.*!PI)]
    ddloc = a + b*i + c*[1 - COS(2.*!PI*i / fn)]
    
    ; Fit to dist = a*i^3 + b*i^2 + c*i
    ;c = ydown_dist*2.
    ;b = 3.*(d/fn^2 - c/fn) - 2.*yup_dist/fn + c/fn
    ;a = d/fn^3 - c/fn^2 - b/fn
    ;dloc = ydown_dist + c*i + b*i^2 + a*i^3
    
  ENDIF ELSE BEGIN
    PRINT, "SORRY; Need 2 points in each region"
    STOP
  ENDELSE

  ; Get indices in ri, zi
  ind = INTERPOL(FINDGEN(np), dist, dloc)

  RETURN, ind
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Create a grid around a given line
;
; Arguments:
;   interp_data  Data used for interpolation of psi(x,y)
;   R, Z         
;   ri, zi       1D indices into F for starting line
;   f0           Starting f for the line
;   fin, fout    Range of f to grid 
;   npar         Number of perpendicular lines to generate
;   nin, nout    Number of points inside and outside line
FUNCTION grid_region, interp_data, R, Z, $
                      ri, zi, $       ; Starting line to grid.
                      fvals, $        ; Location of the surfaces
                      sind, $         ; Index in fvals of the starting line
                      npar, $         ; Number of points along the line
                      slast=slast, $  ; Index in fvals of last successful point
                      sfirst=sfirst, $
                      oplot=oplot, $
                      boundary=boundary, $
                      ffirst=ffirst, flast=flast, $
                      fpsi=fpsi, $ ; f(psi) = R*Bt optional current function
                      parweight=parweight, $ ; Space equally in parallel (1) or poloidal (0) distance
                      ydown_dist=ydown_dist, yup_dist=yup_dist, $
                      ydown_space=ydown_space, yup_space=yup_space
  
  nsurf = N_ELEMENTS(fvals)
  
  IF sind GE 0 THEN BEGIN
    ; starting position is on one of the output surfaces
    f0 = fvals[sind]
    nin = sind
  ENDIF ELSE BEGIN
    ; Starting position between surfaces
    n = FIX(N_ELEMENTS(ri)/2)
    local_gradient, interp_data, ri[n], zi[n], status=status, f=f0
    
    IF fvals[0] LT fvals[nsurf-1] THEN BEGIN
      w = WHERE(fvals LE f0, nin)
    ENDIF ELSE BEGIN
      w = WHERE(fvals GE f0, nin)
    ENDELSE
  ENDELSE
  nout = nsurf - nin - 1

  sfirst = 0      ; Innermost successful index
  slast = nsurf-1 ; Last successful index 

  ffirst = fvals[sfirst]
  flast = fvals[slast]

  PRINT, "    => Gridding range: ", MIN(fvals), MAX(fvals)
  
  nr = interp_data.nx
  nz = interp_data.ny

  ind = poloidal_grid(interp_data, R, Z, ri, zi, npar, fpsi=fpsi, $
                      ydown_dist=ydown_dist, yup_dist=yup_dist, $
                      ydown_space=ydown_space, yup_space=yup_space, $
                      parweight=parweight)
  
  rii = INTERPOLATE(ri, ind)
  zii = INTERPOLATE(zi, ind)

  ;rii = int_func(SMOOTH(deriv(rii), 3)) + rii[0]
  ;zii = int_func(SMOOTH(deriv(zii), 3)) + zii[0]
  ;STOP
  
  ; Refine the location of the starting point
  FOR i=0, npar-1 DO BEGIN
    follow_gradient, interp_data, R, Z, rii[i], zii[i], f0, ri1, zi1
    rii[i] = ri1
    zii[i] = zi1
  ENDFOR

  ; From each starting point, follow gradient in both directions
  
  rixy = FLTARR(nsurf, npar)
  zixy = FLTARR(nsurf, npar)
  FOR i=0, npar-1 DO BEGIN
    IF sind GE 0 THEN BEGIN
      rixy[nin, i] = rii[i]
      zixy[nin, i] = zii[i]
    ENDIF ELSE BEGIN
      ; fvals[nin] should be just outside the starting position
      ftarg = fvals[nin]
      follow_gradient, interp_data, R, Z, rii[i], zii[i], $
        ftarg, rinext, zinext, status=status
      rixy[nin, i] = rinext
      zixy[nin, i] = zinext
    ENDELSE
    FOR j=0, nout-1 DO BEGIN
      ftarg = fvals[nin+j+1]
      
      follow_gradient, interp_data, R, Z, rixy[nin+j, i], zixy[nin+j, i], $
        ftarg, rinext, zinext, status=status, $
        boundary=boundary, fbndry=fbndry

      IF status EQ 1 THEN BEGIN
        rixy[nin+j+1, i] = -1.0
        IF nin+j LT slast THEN slast = nin+j ; last good surface index
        fbndry = fvals[slast]
        IF (fvals[1] - fvals[0])*(flast - fbndry) GT 0 THEN flast = 0.95*fbndry + 0.05*f0
        BREAK
      ENDIF ELSE IF status EQ 2 THEN BEGIN
        ; Hit a boundary 
        rixy[nin+j+1, i] = rinext
        zixy[nin+j+1, i] = zinext
        IF nin+j LT slast THEN slast = nin+j ; Set the last point
        IF (fvals[1] - fvals[0])*(flast - fbndry) GT 0 THEN flast = 0.95*fbndry + 0.05*f0
        BREAK
      ENDIF ELSE BEGIN
        rixy[nin+j+1, i] = rinext
        zixy[nin+j+1, i] = zinext
      ENDELSE
    ENDFOR
    FOR j=0, nin-1 DO BEGIN
      ftarg = fvals[nin-j-1]
      
      follow_gradient, interp_data, R, Z, rixy[nin-j, i], zixy[nin-j, i], $
        ftarg, rinext, zinext, status=status, $
        boundary=boundary, fbndry=fbndry
      
      IF status EQ 1 THEN BEGIN
        rixy[nin-j-1, i] = -1.0
        IF nin-j GT sfirst THEN sfirst = nin-j
        fbndry = fvals[sfirst]
        IF (fvals[1] - fvals[0])*(ffirst - fbndry) LT 0 THEN ffirst = 0.95*fbndry + 0.05*f0
        BREAK
      ENDIF

      rixy[nin-j-1, i] = rinext
      zixy[nin-j-1, i] = zinext

      IF status EQ 2 THEN BEGIN
        IF nin-j GT sfirst THEN sfirst = nin-j
        IF (fvals[1] - fvals[0])*(ffirst - fbndry) LT 0 THEN ffirst = 0.95*fbndry + 0.05*f0
        BREAK
      ENDIF
    ENDFOR
    IF KEYWORD_SET(oplot) THEN BEGIN
      OPLOT, INTERPOLATE(R, rixy[*, i]), INTERPOLATE(Z, zixy[*, i]), color=4
    ENDIF
  ENDFOR

  RETURN, {rixy:rixy, zixy:zixy, rxy:INTERPOLATE(R, rixy), zxy:INTERPOLATE(Z, zixy)}
END

PRO plot_grid_section, a, _extra=_extra
  s = SIZE(a.rxy, /dimension)
  nrad = s[0]
  npol = s[1]
  
  FOR j=0, nrad-1 DO BEGIN
    OPLOT, a.Rxy[j,*], a.Zxy[j,*], _extra=_extra
  ENDFOR
  
  FOR j=0, npol-1 DO BEGIN
    OPLOT, a.Rxy[*,j], a.Zxy[*,j], _extra=_extra
  ENDFOR
END

; Plot a line from a given starting point to a given f
PRO oplot_line, interp_data, R, Z, ri0, zi0, fto, npt=npt, color=color, _extra=_extra
  
  IF NOT KEYWORD_SET(color) THEN color=1
  
  pos = get_line(interp_data, R, Z, ri0, zi0, fto, npt=npt)
  
  OPLOT, INTERPOLATE(R, pos[*,0]), INTERPOLATE(Z, pos[*,1]), $
    color=color, _extra=_extra
END

FUNCTION line_dist, R, Z, ri, zi
  drdi = DERIV(INTERPOLATE(R, ri))
  dzdi = DERIV(INTERPOLATE(Z, zi))
  dldi = SQRT(drdi^2 + dzdi^2)
  RETURN, int_func(findgen(N_ELEMENTS(dldi)), dldi)
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
; Calculate poloidal distances around x-point

FUNCTION xpt_hthe, dctF, R, Z, sep_info, dist, pf_f, core_f, sol_in_f, sol_out_f, boundary=boundary, psi=psi

  ;;;;;;; Leg 1

  darr = sep_info.leg1_dist
  ind = INTERPOL(FINDGEN(N_ELEMENTS(darr)), darr, dist[0])
  ri = INTERPOLATE(sep_info.leg1_ri, ind)
  zi = INTERPOLATE(sep_info.leg1_zi, ind)

  ; Go inwards to PF
  follow_gradient, dctF, R, Z, ri, zi, $
                   pf_f, pf_ri1, pf_zi1, status=status, fbndry=fbndry, boundary=boundary, psi=psi
  IF status THEN pf_f = fbndry
  
  ; Outwards into in SOL
  follow_gradient, dctF, R, Z, ri, zi, $
                   sol_in_f, sol_in_ri1, sol_in_zi1, status=status, fbndry=fbndry, boundary=boundary, psi=psi
  IF status THEN sol_in_f = fbndry
  
  ;;;;;;; Leg 2
  
  darr = sep_info.leg2_dist
  ind = INTERPOL(FINDGEN(N_ELEMENTS(darr)), darr, dist[3])
  ri = INTERPOLATE(sep_info.leg2_ri, ind)
  zi = INTERPOLATE(sep_info.leg2_zi, ind)
  
  ; Go inwards to PF
  follow_gradient, dctF, R, Z, ri, zi, $
                   pf_f, pf_ri2, pf_zi2, status=status, fbndry=fbndry, boundary=boundary, psi=psi
  IF status THEN pf_f = fbndry
  ; Outwards into in SOL
  follow_gradient, dctF, R, Z, ri, zi, $
                   sol_out_f, sol_out_ri1, sol_out_zi1, status=status, fbndry=fbndry, boundary=boundary, psi=psi
  IF status THEN sol_out_f = fbndry
  ;;;;;;; Core 1
  
  darr = sep_info.core1_dist
  ind = INTERPOL(FINDGEN(N_ELEMENTS(darr)), darr, dist[2])
  ri = INTERPOLATE(sep_info.core1_ri, ind)
  zi = INTERPOLATE(sep_info.core1_zi, ind)
  
  ; Go inwards to core
  follow_gradient, dctF, R, Z, ri, zi, $
                   core_f, core_ri1, core_zi1, status=status, fbndry=fbndry, boundary=boundary, psi=psi
  g = local_gradient(dctF, core_ri1, core_zi1, f=psi)
  
  IF status NE 0 THEN STOP
  ; Outwards into out SOL
  follow_gradient, dctF, R, Z, ri, zi, $
                   sol_out_f, sol_out_ri2, sol_out_zi2, status=status, fbndry=fbndry, boundary=boundary, psi=psi
  IF status THEN sol_out_f = fbndry
  
  ;;;;;;; Core 2
  
  darr = sep_info.core2_dist
  ind = INTERPOL(FINDGEN(N_ELEMENTS(darr)), darr, dist[1])
  ri = INTERPOLATE(sep_info.core2_ri, ind)
  zi = INTERPOLATE(sep_info.core2_zi, ind)
  
  ; Go inwards to core
  follow_gradient, dctF, R, Z, ri, zi, $
                   core_f, core_ri2, core_zi2, status=status, fbndry=fbndry, boundary=boundary, psi=psi
  IF status THEN STOP ;core_f = fbndry
  ; Outwards into in SOL
  follow_gradient, dctF, R, Z, ri, zi, $
                   sol_in_f, sol_in_ri2, sol_in_zi2, status=status, fbndry=fbndry, boundary=boundary, psi=psi
  IF status THEN sol_out_f = fbndry
  
  ; Get distances between end points
  
  OPLOT, INTERPOLATE(R, [pf_ri1, pf_ri2]), $
    INTERPOLATE(Z, [pf_zi1, pf_zi2]), thick=2, color=1

  d1 = (INTERPOLATE(R, pf_ri1) - INTERPOLATE(R, pf_ri2))^2 + $
       (INTERPOLATE(Z, pf_zi1) - INTERPOLATE(Z, pf_zi2))^2
  
  d2 = (INTERPOLATE(R, core_ri1) - INTERPOLATE(R, core_ri2))^2 + $
       (INTERPOLATE(Z, core_zi1) - INTERPOLATE(Z, core_zi2))^2
  
  OPLOT, INTERPOLATE(R, [core_ri1, core_ri2]), $
    INTERPOLATE(Z, [core_zi1, core_zi2]), thick=2, color=2

  d3 = (INTERPOLATE(R, sol_in_ri1) - INTERPOLATE(R, sol_in_ri2))^2 + $
       (INTERPOLATE(Z, sol_in_zi1) - INTERPOLATE(Z, sol_in_zi2))^2

  OPLOT, INTERPOLATE(R, [sol_in_ri1, sol_in_ri2]), $
    INTERPOLATE(Z, [sol_in_zi1, sol_in_zi2]), thick=2, color=3

  d4 = (INTERPOLATE(R, sol_out_ri1) - INTERPOLATE(R, sol_out_ri2))^2 + $
       (INTERPOLATE(Z, sol_out_zi1) - INTERPOLATE(Z, sol_out_zi2))^2

  OPLOT, INTERPOLATE(R, [sol_out_ri1, sol_out_ri2]), $
    INTERPOLATE(Z, [sol_out_zi1, sol_out_zi2]), thick=2, color=4
  
  ;cursor, x, y, /down

  RETURN, SQRT([d1, d3, d2, d4]) ; Clockwise around from PF
END

FUNCTION xpt_hthe_newt, dist
  COMMON xpt_newt_com, F, dctF, R, Z, sep_info, pf_f, core_f, sol_in_f, sol_out_f, xd_target, bndry
  RETURN, xpt_hthe(dctF, R, Z, sep_info, dist, pf_f, core_f, sol_in_f, sol_out_f, boundary=bndry, psi=F) - xd_target
END

FUNCTION solve_xpt_hthe, dctF, R, Z, sep_info, dist0, pf_f, core_f, sol_in_f, sol_out_f, boundary=boundary, psi=psi
  COMMON xpt_newt_com, Fc, dctFc, Rc, Zc, sep_infoc, pf_fc, core_fc, sol_in_fc, sol_out_fc, xd_target, bndry

  dist = dist0

  IF KEYWORD_SET(psi) THEN Fc = psi ELSE Fc = 0
  dctFc = dctF
  Rc = R
  Zc = Z
  sep_infoc = sep_info
  pf_fc = pf_f
  core_fc = core_f
  sol_in_fc = sol_in_f
  sol_out_fc = sol_out_f
  IF KEYWORD_SET(boundary) THEN bndry = boundary
  
  xd = xpt_hthe(dctF, R, Z, sep_info, dist, pf_f, core_f, sol_in_f, sol_out_f, boundary=boundary, psi=psi) 
  xd_target = MAX(xd) ;MEAN(xd)
  
  ; NOTE: IDL throws an internal error 
  ; Internal error: Bad variable type encountered in no_name_var().
  ; when NEWTON is combined with LSODE

  dfdx = FLTARR(4,4)
  delta = 1e-3

  REPEAT BEGIN
  xp0 = xpt_hthe_newt(dist)
  response = FLTARR(4)
  FOR i=0, 3 DO BEGIN
    ; Calculate partial derivatives using finite-differences
    d = FLTARR(4)
    d[i] = delta
    dfdx[*,i] = (xpt_hthe_newt(dist+d) - xp0) / delta
    response[i] = MIN([dfdx[i,i], dfdx[(i+1) MOD 4,i]])
  ENDFOR

  IF MIN(dfdx) LT 0.0 THEN BEGIN
    PRINT, "WARNING: ILL-BEHAVED FITTING"
    RETURN, dist0 ; Don't modify
  ENDIF
  
;  ; Invert using SVD
;  SVDC, dfdx, W, U, V
;  WP = FLTARR(4, 4)
;  for i=0,2 do wp[i,i] = 1.0/w[i]
;  ddist = V ## WP ## TRANSPOSE(U) # xp0
;  
;  ;ddist = INVERT(dfdx) # xp0
;  w = WHERE(ABS(ddist) GT 0.5*dist, count)
;  IF count GT 0 THEN ddist[w] = ddist[w] * 0.5*dist[w] / ABS(ddist[w])
;  dist = dist - ddist
;  
  PRINT, "DIST =", REFORM(dist)
  PRINT, "RESP = ", response
;  PRINT, "CHANGE = ", ddist
  PRINT, "XP0 = ", xp0

  ; Check for responses which are too high or too low
  ; If too low, move out, too high move in

  med = MEDIAN(response)
  w = WHERE(response LT 0.1*med, count1)
  IF count1 GT 0 THEN dist[w] = dist[w] * 2.
  w = WHERE(response GT 10.*med, count2)
  IF count2 GT 0 THEN dist[w] = dist[w] * 0.75
  
  ENDREP UNTIL count1+count2 EQ 0
  
  RETURN, dist
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Try to equalise x-point separations

FUNCTION increase_xpt_hthe, dctF, R, Z, sep_info, dist0, pf_f, core_f, sol_in_f, sol_out_f, boundary=boundary, psi=psi

  dist = dist0
  
  ; Get the initial distances
  xd = xpt_hthe(dctF, R, Z, sep_info, dist, pf_f, core_f, sol_in_f, sol_out_f, boundary=boundary, psi=psi) 
  xd_target = MEAN(xd)
  ; Numbering of xd: 0 - PF, going clockwise

  PRINT, "START:"+STR(dist)+" -> "+STR(xd)
  
  ; Increase spacing until all above a minimum
  REPEAT BEGIN
    m = MIN(xd, ind)
    IF m LT xd_target THEN BEGIN
      ; Try to increase this distance
      ; Can adjust dist[ind] or dist[ind-1]
      
      IF xd[(ind+1) MOD 4] LT xd[(ind+3) MOD 4] THEN BEGIN
        ; Increase dist[ind]
        dist[ind] = dist[ind] * 1.1
      ENDIF ELSE BEGIN
        ; Increase dist[ind-1]
        dist[(ind+3) MOD 4] = dist[(ind+3) MOD 4] * 1.1
      ENDELSE
    ENDIF
    xd = xpt_hthe(dctF, R, Z, sep_info, dist, pf_f, core_f, sol_in_f, sol_out_f, boundary=boundary, psi=psi)
  ENDREP UNTIL (MIN(xd) GE xd_target) OR ( ABS(MIN(xd) - m) LT 1.e-5 )
  
  ; Reduce spacing until one goes below target
  REPEAT BEGIN
    dist_last = dist
    m = MAX(xd, ind)
    IF m GT xd_target THEN BEGIN
      ; Decrease distance 
      IF xd[(ind+1) MOD 4] GT xd[(ind+3) MOD 4] THEN BEGIN
        dist[ind] = dist[ind] / 1.1
      ENDIF ELSE BEGIN
        dist[(ind+3) MOD 4] = dist[(ind+3) MOD 4] / 1.1
      ENDELSE
    ENDIF
    xd = xpt_hthe(dctF, R, Z, sep_info, dist, pf_f, core_f, sol_in_f, sol_out_f, boundary=boundary, psi=psi)
  ENDREP UNTIL MIN(xd) LT xd_target
  dist = dist_last
  
  PRINT, "ADJUST:"+STR(dist)+" -> "+STR(xd)
  
  RETURN, dist
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Main gridding function
;

FUNCTION create_grid, F, R, Z, in_settings, critical=critical, $
                      boundary=boundary, debug=debug, strictbndry=strictbndry, iter=iter, $
                      fpsi = fpsi, $ ; f(psi) = R*Bt current function
                      nrad_flexible=nrad_flexible, $
                      single_rad_grid=single_rad_grid, fast=fast, $
                      xpt_mindist=xpt_mindist, xpt_mul=xpt_mul

  IF SIZE(nrad_flexible, /TYPE) EQ 0 THEN nrad_flexible = 0

  ; Create error handler
  err=0;CATCH, err
  IF err NE 0 THEN BEGIN
    PRINT, "CREATE_GRID failed"
    PRINT, "   Error message: "+!ERROR_STATE.MSG
    CATCH, /cancel
    RETURN, {error:1}
  ENDIF

  IF NOT KEYWORD_SET(iter) THEN iter = 0
  IF iter GT 3 THEN BEGIN
    PRINT, "ERROR: Too many iterations"
    RETURN, {error:1}
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Check the settings
  ; If a setting is missing, set a default value
  
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
                pol_peaking:0.0, $
                parweight:0.0}
  ENDIF ELSE BEGIN
    PRINT, "Checking settings"
    settings = in_settings ; So the input isn't changed
    str_check_present, settings, 'psi_inner', 0.9
    str_check_present, settings, 'psi_outer', 1.1
    str_check_present, settings, 'nrad', 36
    str_check_present, settings, 'npol', 64
    str_check_present, settings, 'rad_peaking', 0.0
    str_check_present, settings, 'pol_peaking', 0.0
    str_check_present, settings, 'parweight', 0.0
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

  ; Get an even number of points for efficient FFTs
  IF nx MOD 2 EQ 1 THEN BEGIN
    ; odd number of points in R. Cut out last point
    R = R[0:(nx-2)]
    F = F[0:(nx-2), *]
    nx = nx - 1
  ENDIF
  IF ny MOD 2 EQ 1 THEN BEGIN
    ; Same for Z
    Z = Z[0:(ny-2)]
    F = F[*,0:(ny-2)]
    ny = ny - 1
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
  
  IF NOT KEYWORD_SET(bndryi) THEN BEGIN
    bndryi = FLTARR(2,4)
    bndryi[0,*] = [1, nx-2, nx-2, 1]
    bndryi[1,*] = [1, 1, ny-2, ny-2]
  ENDIF
  
  ;;;;;;;;;;;;;;; Psi interpolation data ;;;;;;;;;;;;;;
  
  interp_data = {nx:nx, ny:ny, $
                 method:0, $
                 f: F}       ; Always include function
               
  IF KEYWORD_SET(fast) THEN BEGIN
    PRINT, "Using Fast settings"
    interp_data.method = 2
  ENDIF
  
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
    refine_xpoint, interp_data, xpt_ri[i], xpt_zi[i], ri, zi
    local_gradient, interp_data, ri, zi, f=val
    PRINT, "   -> "+STR([ri, zi])+" "+STR(val)
    xpt_ri[i] = ri
    xpt_zi[i] = zi
    xpt_f[i] = val
  ENDFOR

  ; Overplot the separatrices, O-points
  oplot_critical, F, R, Z, critical

  ; Psi normalisation factors
  faxis = critical.opt_f[critical.primary_opt]
  fnorm = critical.xpt_f[critical.inner_sep] - critical.opt_f[critical.primary_opt]

  ; From normalised psi, get range of f
  f_inner = faxis + MIN(settings.psi_inner)*fnorm
  f_outer = faxis + MAX(settings.psi_outer)*fnorm
  
  ; Check the number of x-points
  IF critical.n_xpoint EQ 0 THEN BEGIN
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Grid entirely in the core

    PRINT, "Generating grid entirely in the core"

    nrad = TOTAL(settings.nrad,/int) ; Add up all points
    npol = TOTAL(settings.npol,/int)
    rad_peaking = settings.rad_peaking[0] ; Just the first region
    pol_peaking = settings.pol_peaking[0]

    ; work out where to put the surfaces in the core
    fvals = radial_grid(nrad, f_inner, f_outer, 1, 1, [xpt_f[inner_sep]], rad_peaking)

    ; Create a starting surface
    sind = FIX(nrad / 2)
    start_f = fvals[sind]
    
    contour_lines, F, findgen(nx), findgen(ny), levels=[start_f], $
      path_info=info, path_xy=xy
  
    IF N_ELEMENTS(info) GT 1 THEN BEGIN
      ; Find the surface closest to the o-point
      
      ind = closest_line(info, xy, opt_ri[primary_opt], opt_zi[primary_opt])
      info = info[ind]
    ENDIF ELSE info = info[0]
    
    start_ri = REFORM(xy[0,info.offset:(info.offset+info.n-1)])
    start_zi = REFORM(xy[1,info.offset:(info.offset+info.n-1)])

    ; Make sure that the line goes clockwise
    
    m = MAX(INTERPOLATE(Z, start_zi), ind)
    IF (DERIV(INTERPOLATE(R, start_ri)))[ind] LT 0.0 THEN BEGIN
      ; R should be increasing at the top. Need to reverse
      start_ri = REVERSE(start_ri)
      start_zi = REVERSE(start_zi)
    ENDIF
    
    ; Last point should be the same as the first
    start_ri = [start_ri, start_ri[0]]
    start_zi = [start_zi, start_zi[0]]
    
    ; Smooth and refine the starting location
    np = N_ELEMENTS(start_ri)
    s = 3
    start_ri = (SMOOTH([start_ri[(np-s):(np-1)], start_ri, start_ri[0:(s-1)]], s))[s:(np-1+s)]
    start_zi = (SMOOTH([start_zi[(np-s):(np-1)], start_zi, start_zi[0:(s-1)]], s))[s:(np-1+s)]
    
    FOR i=0, np-1 DO BEGIN
      follow_gradient, interp_data, R, Z, start_ri[i], start_zi[i], start_f, ri1, zi1
      start_ri[i] = ri1
      start_zi[i] = zi1
    ENDFOR
    
    oplot_contour, info, xy, R, Z, /periodic, color=2, thick=1.5
    
    a = grid_region(interp_data, R, Z, $
                    start_ri, start_zi, $
                    fvals, $
                    sind, $
                    npol, fpsi=fpsi, $
                    parweight=settings.parweight, $
                    /oplot)
    
    OPLOT, [REFORM(a.rxy[0,*]), a.rxy[0,0]], [REFORM(a.zxy[0,*]), a.zxy[0,0]], color=4
      
    FOR i=1, nrad-1 DO BEGIN
      OPLOT, [REFORM(a.rxy[i,*]), a.rxy[i,0]], [REFORM(a.zxy[i,*]), a.zxy[i,0]], color=4
    ENDFOR

    FOR i=0, npol-1 DO BEGIN
      OPLOT, a.rxy[*,i], a.zxy[*,i], color=4
    ENDFOR
    
    ; Get other useful variables
    Psixy = FLTARR(nrad, npol)
    FOR i=0, npol-1 DO psixy[*,i] = (fvals - faxis)/fnorm ; to get normalised psi
    
    ; Calculate magnetic field components
    dpsidR = FLTARR(nrad, npol)
    dpsidZ = dpsidR
    
    interp_data.method = 2

    FOR i=0,nrad-1 DO BEGIN
      FOR j=0,npol-1 DO BEGIN
        local_gradient, interp_data, a.Rixy[i,j], a.Zixy[i,j], status=status, $
          dfdr=dfdr, dfdz=dfdz
        ; dfd* are derivatives wrt the indices. Need to multiply by dr/di etc
        dpsidR[i,j] = dfdr/INTERPOLATE(DERIV(R),a.Rixy[i,j]) 
        dpsidZ[i,j] = dfdz/INTERPOLATE(DERIV(Z),a.Zixy[i,j]) 
      ENDFOR
    ENDFOR

    ; Set topology to connect in the core
    yup_xsplit = [nrad]
    ydown_xsplit = [nrad]
    yup_xin = [0]
    yup_xout = [-1]
    ydown_xin = [0]
    ydown_xout = [-1]

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Create result structure

    result = {error:0, $ ; Signal success
              psi_inner:settings.psi_inner, psi_outer:settings.psi_outer, $ ; Unchanged psi range
              nrad:nrad, npol:npol, $  ; Number of points in radial and poloidal direction
              Rixy:a.Rixy, Zixy:a.Zixy, $  ; Indices into R and Z of each point
              Rxy:a.Rxy, Zxy:a.Zxy, $  ; Location of each grid point
              psixy:psixy, $ ; Normalised psi for each point
              dpsidR:dpsidR, dpsidZ:dpsidZ, $ ; Psi derivatives (for Bpol)
              faxis:faxis, fnorm:fnorm, $ ; Psi normalisation factors
              settings:settings, $ ; Settings used to create grid
              critical:critical, $ ; Critical points
              yup_xsplit:yup_xsplit, $ ; X index where domain splits (number of points in xin)
              ydown_xsplit:ydown_xsplit, $
              yup_xin:yup_xin, yup_xout:yup_xout, $ ; Domain index to connect to
              ydown_xin:ydown_xin, ydown_xout:ydown_xout}

    RETURN, result
    
  ENDIF ELSE BEGIN

    rad_peaking = settings.rad_peaking
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Grid contains at least one x-point

    ; Normalised psi value of each separatrix
    xpt_psi = (critical.xpt_f - faxis) / fnorm
    
    si = SORT(xpt_psi) ; Sort separatrices from inside out

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Choose primary x-point.
    ; Determines order of settings arrays
    primary_xpt = si[0]

    PRINT, "Primary X-point is number "+STR(primary_xpt)
    PRINT, "   at R = "+STR(INTERPOLATE(R, critical.xpt_ri[primary_xpt])) $
      +" Z = "+STR(INTERPOLATE(Z, critical.xpt_zi[primary_xpt]))
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; work out where to put the surfaces

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
      
      ;fvals = radial_grid(TOTAL(nrad), f_inner, f_outer, 1, 1, xpt_f, rad_peaking)
      
      
    ENDIF
    
    ; Now got number of points in each region
    
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
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Create arrays of psi values for each region
    
    sol_psi_vals = FLTARR(critical.n_xpoint, TOTAL(nrad,/int))
    pf_psi_vals  = FLTARR(critical.n_xpoint, 2, TOTAL(nrad,/int))
    FOR i=0, critical.n_xpoint-1 DO BEGIN
      sol_psi_vals[i,*]  = psi_vals
      pf_psi_vals[i,0,*] = psi_vals
      pf_psi_vals[i,1,*] = psi_vals
    ENDFOR
    
    IF N_ELEMENTS(settings.psi_inner) EQ (critical.n_xpoint + 1) THEN BEGIN
      psi_inner = settings.psi_inner
    ENDIF ELSE BEGIN
      IF N_ELEMENTS(settings.psi_inner) GT 1 THEN BEGIN
        PRINT, "WARNING: psi_inner has wrong number of elements ("+STR(N_ELEMENTS(settings.psi_inner))+")"
        PRINT, "         Should have 1 or "+STR(critical.n_xpoint+1)+" elements"
      ENDIF
      PRINT, "Keeping same inner psi for all regions"
      
      psi_inner = FLTARR(critical.n_xpoint+1) + MIN(settings.psi_inner)
    ENDELSE
    
    IF N_ELEMENTS(settings.psi_outer) EQ critical.n_xpoint THEN BEGIN
      psi_outer = settings.psi_outer
    ENDIF ELSE BEGIN
      IF N_ELEMENTS(settings.psi_outer) GT 1 THEN BEGIN
        PRINT, "WARNING: psi_outer has wrong number of elements ("+STR(N_ELEMENTS(settings.psi_outer))+")"
        PRINT, "         Should have 1 or "+STR(critical.n_xpoint)+" elements"
      ENDIF 
      PRINT, "Keeping same outer psi for all regions"
      
      psi_outer = FLTARR(critical.n_xpoint) + MAX(settings.psi_outer)
    ENDELSE
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Need to work out where the x-points are
    ; and create lines just inside the separatrices
    
    ; take a contour surface just inside the innermost separatrix
    ;sind = FIX(nrad[0]/2) ;nrad[0]-1 ; the last point inside core
    ;sind = nrad[0]-1
    ;f_cont = fvals[sind]
    f_cont = faxis + fnorm*(0.1*psi_inner[0] + 0.9)

    contour_lines, F, findgen(nx), findgen(ny), levels=[f_cont], $
      path_info=info, path_xy=xy
    
    IF N_ELEMENTS(info) GT 0 THEN BEGIN
      ; Find the surface closest to the o-point

      info = info[closest_line(info, xy, $
                               critical.opt_ri[critical.primary_opt], critical.opt_zi[critical.primary_opt])]
    ENDIF ELSE info = info[0]
    
    oplot_contour, info, xy, R, Z, /periodic, color=3, thick=1.5

    start_ri = REFORM(xy[0,info.offset:(info.offset+info.n-1)])
    start_zi = REFORM(xy[1,info.offset:(info.offset+info.n-1)])

    ; Make sure that the line goes clockwise
    
    m = MAX(INTERPOLATE(Z, start_zi), ind)
    IF (DERIV(INTERPOLATE(R, start_ri)))[ind] LT 0.0 THEN BEGIN
      ; R should be increasing at the top. Need to reverse
      start_ri = REVERSE(start_ri)
      start_zi = REVERSE(start_zi)
    ENDIF

    ; now have (start_ri, start_zi). For each x-point, find the radial
    ; line going through the x-point
    
    fri = FFT(start_ri) ; for interpolating periodic functions
    fzi = FFT(start_zi)

    xpt_ind = FLTARR(critical.n_xpoint)  ; index into start_*i
    
    pf_info = PTRARR(critical.n_xpoint)
    
    
    FOR i=0, critical.n_xpoint-1 DO BEGIN
      PRINT, "Finding theta location of x-point "+STR(i)
      
      ; Get the separatrices
      legsep = leg_separatrix2(interp_data, R, Z, xpt_ri[i], xpt_zi[i], $
                               opt_ri[primary_opt], opt_zi[primary_opt], boundary=bndryi)
      
      ; Go a little way along each core separatrix and follow
      follow_gradient, interp_data, R, Z, $
                       legsep.core1[2,0], legsep.core1[2,1], $
                       0.95 * f_cont + 0.05*opt_f[primary_opt], $
                       rhit, zhit, $
                       boundary=TRANSPOSE([[start_ri], [start_zi]]), ibndry=hit_ind1
      follow_gradient, interp_data, R, Z, $
                       legsep.core2[2,0], legsep.core2[2,1], $
                       0.95 * f_cont + 0.05*opt_f[primary_opt], $
                       rhit, zhit, $
                       boundary=TRANSPOSE([[start_ri], [start_zi]]), ibndry=hit_ind2
      
      ni = N_ELEMENTS(start_ri)
      
      ; Refine the theta index of the X-point using divide and conquer
      REPEAT BEGIN
        IF MIN([ni - hit_ind2 + hit_ind1, ni - hit_ind1 + hit_ind2]) LT ABS(hit_ind2 - hit_ind1) THEN BEGIN
          ; One at the beginning and one at the end (across the join)
          mini = (hit_ind2 + hit_ind1 - ni) / 2.
          IF mini LT 0. THEN mini = mini + ni
        ENDIF ELSE mini = (hit_ind1 + hit_ind2) / 2.
        
        ;OPLOT, [INTERPOLATE(R[start_ri], hit_ind1)], [INTERPOLATE(Z[start_zi], hit_ind1)], psym=2, color=2
        ;OPLOT, [INTERPOLATE(R[start_ri], hit_ind2)], [INTERPOLATE(Z[start_zi], hit_ind2)], psym=2, color=2
        ;OPLOT, [INTERPOLATE(R[start_ri], mini)], [INTERPOLATE(Z[start_zi], mini)], psym=2, color=4
        
        PRINT, "Theta location: " + STR(hit_ind1) + "," + STR(hit_ind2) + " -> " + STR(mini)

        ;  Get line a little bit beyond the X-point
        pos = get_line(interp_data, R, Z, $
                       INTERPOLATE(start_ri, mini), INTERPOLATE(start_zi, mini), $
                       critical.xpt_f[i] + (critical.xpt_f[i] - opt_f[primary_opt]) * 0.05)
        
        ;OPLOT, INTERPOLATE(R, pos[*,0]), INTERPOLATE(Z, pos[*,1]), color=4, thick=2
        
        ; Find which separatrix line this intersected with
        cpos = line_crossings([xpt_ri[i], legsep.core1[*,0]], $
                              [xpt_zi[i], legsep.core1[*,1]], 0, $
                              pos[*,0], pos[*,1], 0, $
                              ncross=ncross, inds1=inds)
        IF ncross GT 0 THEN BEGIN
          hit_ind1 = mini
        ENDIF ELSE BEGIN
          hit_ind2 = mini
        ENDELSE
        dist = MIN([ni - hit_ind2 + hit_ind1, ni - hit_ind1 + hit_ind2, ABS([hit_ind2 - hit_ind1])])
      ENDREP UNTIL dist LT 0.1
      IF MIN([ni - hit_ind2 + hit_ind1, ni - hit_ind1 + hit_ind2]) LT ABS(hit_ind2 - hit_ind1) THEN BEGIN
        ; One at the beginning and one at the end (across the join)
        mini = (hit_ind2 + hit_ind1 - ni) / 2.
        IF mini LT 0. THEN mini = mini + ni
      ENDIF ELSE mini = (hit_ind1 + hit_ind2) / 2.
        
      xpt_ind[i] = mini  ; Record the index
        
      ; Plot the line to the x-point
      oplot_line, interp_data, R, Z, $
        fft_interp(fri, mini), fft_interp(fzi, mini), critical.xpt_f[i]
      oplot_line, interp_data, R, Z, $
        fft_interp(fri, mini), fft_interp(fzi, mini), f_inner

      ; Get tangent vector
      drdi = INTERPOLATE((DERIV(INTERPOLATE(R, start_ri))), mini)
      dzdi = INTERPOLATE((DERIV(INTERPOLATE(Z, start_zi))), mini)
      
      tmp = {core_ind:mini, drdi:drdi, dzdi:dzdi, $ ; Core index and tangent vector
             sol:LONARR(2)} ; Array to store SOL indices
      
      pf_info[i] = PTR_NEW(tmp)
    ENDFOR

    ; Sort x-points by index along this core surface
    ci = SORT(xpt_ind)
    
    ; Extract starting lines for each section
    sol_info = PTRARR(critical.n_xpoint)
    FOR i=0, critical.n_xpoint-1 DO BEGIN
      IF i NE (critical.n_xpoint-1) THEN BEGIN
        ri = [ fft_interp(fri,xpt_ind[ci[i]]), $
               start_ri[FIX(xpt_ind[ci[i]]+1.0):FIX(xpt_ind[ci[i+1]])], $
               fft_interp(fri,xpt_ind[ci[i+1]]) ]
        
        zi = [ fft_interp(fzi,xpt_ind[ci[i]]), $
               start_zi[FIX(xpt_ind[ci[i]]+1.0):FIX(xpt_ind[ci[i+1]])], $
               fft_interp(fzi,xpt_ind[ci[i+1]]) ]
      ENDIF ELSE BEGIN
        ; Index wraps around
        IF xpt_ind[ci[i]] GT N_ELEMENTS(start_ri)-2 THEN BEGIN
          ri = [ fft_interp(fri,xpt_ind[ci[i]]), $
                 start_ri[0:FIX(xpt_ind[ci[0]])], $
                 fft_interp(fri,xpt_ind[ci[0]]) ]
          
          zi = [ fft_interp(fzi,xpt_ind[ci[i]]), $
                 start_zi[0:FIX(xpt_ind[ci[0]])], $
                 fft_interp(fzi,xpt_ind[ci[0]]) ]
          
        ENDIF ELSE BEGIN
          ri = [ fft_interp(fri,xpt_ind[ci[i]]), $
                 start_ri[FIX(xpt_ind[ci[i]]+1.0):*], $
                 start_ri[0:FIX(xpt_ind[ci[0]])], $
                 fft_interp(fri,xpt_ind[ci[0]]) ]
          
          zi = [ fft_interp(fzi,xpt_ind[ci[i]]), $
                 start_zi[FIX(xpt_ind[ci[i]]+1.0):*], $
                 start_zi[0:FIX(xpt_ind[ci[0]])], $
                 fft_interp(fzi,xpt_ind[ci[0]]) ]
        ENDELSE
      ENDELSE
      
      ; Calculate length of the line
      drdi = DERIV(INTERPOLATE(R, ri))
      dzdi = DERIV(INTERPOLATE(Z, zi))
      dldi = SQRT(drdi^2 + dzdi^2)
      length = INT_TABULATED(findgen(N_ELEMENTS(dldi)), dldi)

      ; Change the grid in the far SOL
      w = WHERE(ci EQ primary_xpt)
      solid = (i - w[0] + critical.n_xpoint) MOD critical.n_xpoint
      xpt_psi_max = MAX([xpt_psi[ci[i]], xpt_psi[ci[(i+1) MOD critical.n_xpoint]]])
      w = WHERE(sol_psi_vals[i,*] GT xpt_psi_max, nsol)
      PRINT, "Generating sol psi: ", xpt_psi_max, psi_outer[solid]
      IF xpt_psi_max GT psi_outer[solid] THEN BEGIN
        PRINT, "**Far SOL cannot include both x-points"
        PRINT, " Probably caused by intersection with boundaries."
        
        IF KEYWORD_SET(strictbndry) THEN BEGIN
          PRINT, "**Re-running, removing the strict keyword"
        ENDIF ELSE BEGIN
          
          IF N_ELEMENTS(boundary[0,*]) GT 4 THEN BEGIN
            PRINT, "**Re-running, simplifying boundary"
            boundary = TRANSPOSE([ [MIN(boundary[0,*]), MAX(boundary[0,*]), MAX(boundary[0,*]), MIN(boundary[0,*])], $
                                   [MIN(boundary[1,*]), MIN(boundary[1,*]), MAX(boundary[1,*]), MAX(boundary[1,*])] ])
          ENDIF ELSE BEGIN
            PRINT, "**Re-running, removing boundary"
            boundary = 0
          ENDELSE
        ENDELSE
          
        IF KEYWORD_SET(nrad_flexible) THEN nrad = TOTAL(nrad,/int) ; Allow nrad to change again

        new_settings = {psi_inner:psi_inner, psi_outer:(max(xpt_psi)+0.02), $
                        nrad:nrad, npol:settings.npol, $
                        rad_peaking:settings.rad_peaking, pol_peaking:settings.pol_peaking}
        RETURN, create_grid(F, R, Z, new_settings, critical=critical, $
                            boundary=boundary, iter=iter+1, nrad_flexible=nrad_flexible, $
                            single_rad_grid=single_rad_grid, fast=fast)
      ENDIF
      dpsi = sol_psi_vals[i,TOTAL(nrad,/int)-nsol-1] - sol_psi_vals[i,TOTAL(nrad,/int)-nsol-2]
      sol_psi_vals[i,(TOTAL(nrad,/int)-nsol):*] = radial_grid(nsol, $
                                                         sol_psi_vals[i,(TOTAL(nrad,/int)-nsol-1)]+dpsi, $
                                                         psi_outer[solid], $
                                                         1, 1, [xpt_psi_max], $
                                                         settings.rad_peaking, $
                                                         in_dp=dpsi)

      tmp = {xpt1:ci[i], xpt2:ci[(i+1) MOD critical.n_xpoint], $ ; X-point indices
             ri:ri, zi:zi, length:length, $  ; R and Z points and the line length
             nsol:nsol}
      
      sol_info[i] = PTR_NEW(tmp)

      ; Add this SOL to the PF structure
      (*(pf_info[ci[i]])).sol[i MOD 2] = i
      (*(pf_info[ci[(i+1) MOD critical.n_xpoint]])).sol[i MOD 2] = i
    ENDFOR
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Do the same for each PF region

    sep_info = PTRARR(critical.n_xpoint) ; Separatrix info

    FOR i=0, critical.n_xpoint-1 DO BEGIN
      ; Change grid in PF regions (inner psi can vary)
      ; psi_inner sorted in separatrix psi from inside out
      xind = si[i]
      
      ; Get number of points in this PF region (npf)
      w = WHERE(psi_vals LT xpt_psi[xind], npf)
      w = WHERE(ci EQ xind)
      id = w[0]
      
      IF KEYWORD_SET(single_rad_grid) THEN BEGIN
        ; Gridding as one region
        IF (npf+1) LT TOTAL(nrad,/int) THEN BEGIN
          dpsi = pf_psi_vals[xind,0,npf+1] - pf_psi_vals[xind,0,npf]
          pf_psi_out = (pf_psi_vals[xind,0,npf] - 0.5*dpsi) < xpt_psi[xind]
          pf_psi_vals[xind,0,0:(npf-1)] = radial_grid(npf, psi_inner[id+1], $
                                                      pf_psi_out, $
                                                      1, 0, $
                                                      [xpt_psi[xind]], settings.rad_peaking, $
                                                      out_dp=dpsi)
        ENDIF ELSE BEGIN
          pf_psi_out = (pf_psi_vals[xind,0,npf] - 0.5*dpsi) < xpt_psi[xind]
          pf_psi_vals[xind,0,0:(npf-1)] = radial_grid(npf, psi_inner[id+1], $
                                                      pf_psi_out, $
                                                      1, 0, $
                                                      [xpt_psi[xind]], settings.rad_peaking)
        ENDELSE
      ENDIF ELSE BEGIN
        ; Gridding in multiple regions. Ensure equal spacing around separatrix
        dpsi = 2.*(pf_psi_vals[xind,0,npf] - xpt_psi[xind])
        pf_psi_out = xpt_psi[xind] - 0.5*dpsi
        
        pf_psi_vals[xind,0,0:(npf-1)] = radial_grid(npf, psi_inner[id+1], $
                                                    pf_psi_out, $
                                                    1, 1, $
                                                    [xpt_psi[xind]], settings.rad_peaking, $
                                                    out_dp=dpsi)
      ENDELSE
      
      pf_psi_vals[xind,1,0:(npf-1)] = pf_psi_vals[xind,0,0:(npf-1)]
      
      
      ; Get the lines from the x-point to the target plates
      legsep = leg_separatrix2(interp_data, R, Z, xpt_ri[i], xpt_zi[i], $
                               opt_ri[primary_opt], opt_zi[primary_opt], boundary=bndryi)
      
      pf_ri = [REVERSE(legsep.leg1[*,0]), xpt_ri[i], legsep.leg2[*,0]]
      pf_zi = [REVERSE(legsep.leg1[*,1]), xpt_zi[i], legsep.leg2[*,1]]
      mini = N_ELEMENTS(legsep.leg1[*,0])
      
      ; Use the tangent vector to determine direction
      ; relative to core and so get direction of positive theta
      
      drdi = INTERPOLATE((DERIV(INTERPOLATE(R, pf_ri))), mini)
      dzdi = INTERPOLATE((DERIV(INTERPOLATE(Z, pf_zi))), mini)
      
      IF drdi * (*pf_info[xind]).drdi + dzdi * (*pf_info[xind]).dzdi GT 0.0 THEN BEGIN
        ; Line is parallel to the core. Need to reverse
        pf_ri = REVERSE(pf_ri)
        pf_zi = REVERSE(pf_zi)
        mini = N_ELEMENTS(pf_ri) - 1. - mini

        ; Structure for x-point grid spacing info
        xpt_sep = {leg1_ri:[xpt_ri[i], legsep.leg2[*,0]], leg1_zi:[xpt_zi[i], legsep.leg2[*,1]], $
                   leg2_ri:[xpt_ri[i], legsep.leg1[*,0]], leg2_zi:[xpt_zi[i], legsep.leg1[*,1]], $
                   core1_ri:[xpt_ri[i], legsep.core2[*,0]], core1_zi:[xpt_zi[i], legsep.core2[*,1]], $
                   core2_ri:[xpt_ri[i], legsep.core1[*,0]], core2_zi:[xpt_zi[i], legsep.core1[*,1]]}
      ENDIF ELSE BEGIN
        xpt_sep = {leg1_ri:[xpt_ri[i], legsep.leg1[*,0]], leg1_zi:[xpt_zi[i], legsep.leg1[*,1]], $
                   leg2_ri:[xpt_ri[i], legsep.leg2[*,0]], leg2_zi:[xpt_zi[i], legsep.leg2[*,1]], $
                   core1_ri:[xpt_ri[i], legsep.core1[*,0]], core1_zi:[xpt_zi[i], legsep.core1[*,1]], $
                   core2_ri:[xpt_ri[i], legsep.core2[*,0]], core2_zi:[xpt_zi[i], legsep.core2[*,1]]}
      ENDELSE

      leg1_dist = line_dist(R, Z, xpt_sep.leg1_ri, xpt_sep.leg1_zi)
      leg2_dist = line_dist(R, Z, xpt_sep.leg2_ri, xpt_sep.leg2_zi)
      core1_dist = line_dist(R, Z, xpt_sep.core1_ri, xpt_sep.core1_zi)
      core2_dist = line_dist(R, Z, xpt_sep.core2_ri, xpt_sep.core2_zi)
      xpt_sep = CREATE_STRUCT(xpt_sep, "leg1_dist", leg1_dist, "leg2_dist", leg2_dist, "core1_dist", core1_dist, "core2_dist", core2_dist)
      sep_info[i] = PTR_NEW(xpt_sep)

      ; Make sure that the PF sections point to the correct SOLs
      w = WHERE(ci EQ xind)
      i1 = w[0]
      IF (*pf_info[xind]).sol[0] NE w[0] THEN BEGIN
        (*pf_info[xind]).sol = REVERSE((*pf_info[xind]).sol)
      ENDIF
      
      ; Copy radial grid in the SOL
      solind = (*pf_info[xind]).sol[0]
      nsol = (*sol_info[solind]).nsol
      pf_psi_vals[xind, 0, (TOTAL(nrad,/int)-nsol):*] = sol_psi_vals[solind, (TOTAL(nrad,/int)-nsol):*]
      
      solind = (*pf_info[xind]).sol[1]
      nsol = (*sol_info[solind]).nsol
      pf_psi_vals[xind, 1, (TOTAL(nrad,/int)-nsol):*] = sol_psi_vals[solind, (TOTAL(nrad,/int)-nsol):*]

      ; Put the starting line into the pf_info structure
      tmp = CREATE_STRUCT(*(pf_info[xind]), $
                          'npf', npf, $ ; Number of radial points in this PF region
                          'ri0', [pf_ri[0:mini], INTERPOLATE(pf_ri, mini)], $
                          'zi0', [pf_zi[0:mini], INTERPOLATE(pf_zi, mini)], $
                          'ri1', [INTERPOLATE(pf_ri, mini), pf_ri[(mini+1):*]], $
                          'zi1', [INTERPOLATE(pf_zi, mini), pf_zi[(mini+1):*]])

      ; Calculate length of each section
      dldi = SQRT(DERIV(INTERPOLATE(R, tmp.ri0))^2 + DERIV(INTERPOLATE(Z, tmp.zi0))^2)
      len0 = INT_TABULATED(FINDGEN(N_ELEMENTS(dldi)), dldi)
      dldi = SQRT(DERIV(INTERPOLATE(R, tmp.ri1))^2 + DERIV(INTERPOLATE(Z, tmp.zi1))^2)
      len1 = INT_TABULATED(FINDGEN(N_ELEMENTS(dldi)), dldi)
      tmp = CREATE_STRUCT(tmp, 'len0', len0, 'len1', len1)
      
      PTR_FREE, pf_info[xind]
      pf_info[xind] = PTR_NEW(tmp)
    ENDFOR
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Decide poloidal resolutions
    ; If automatic, divide up based on length, with a
    ; minimum of 2 points per region
    
    npol = settings.npol
    nnpol = N_ELEMENTS(npol)
    IF nnpol EQ 1 THEN npol = npol[0]
    
    IF nnpol NE 3*critical.n_xpoint THEN BEGIN
      IF nnpol GT 1 THEN BEGIN
        PRINT, "WARNING: npol has wrong number of elements ("+STR(nnpol)+")"
        PRINT, "         Should have 1 or "+STR(3*critical.n_xpoint)+" elements"
        npol = TOTAL(npol,/int)
      ENDIF
      
      IF npol LT 6*critical.n_xpoint THEN BEGIN
        PRINT, "WARNING: supplied npol ("+STR(npol)+") too small"
        npol = 6*critical.n_xpoint
        PRINT, "   => Increasing npol to "+ STR(npol)
      ENDIF
    
      nnpol = npol
      npol = LONARR(3*critical.n_xpoint) + 2
      nnpol = nnpol - 6*critical.n_xpoint ; Extra points to divide up
      
      ; Get lengths
      length = FLTARR(3*critical.n_xpoint)
      FOR i=0, critical.n_xpoint-1 DO BEGIN
        ; PF regions
        length[i]                     = (*pf_info[i]).len0
        length[critical.n_xpoint + i] = (*pf_info[i]).len1
        ; SOL
        length[2*critical.n_xpoint + i] = (*sol_info[i]).length
      ENDFOR
      
      FOR i=0, nnpol-1 DO BEGIN
        ; Add an extra point to the longest length
        
        dl = length / FLOAT(npol)
        dl[0:(2*critical.n_xpoint-1)] = length[0:(2*critical.n_xpoint-1)] $
          / (FLOAT(npol[0:(2*critical.n_xpoint-1)]) - 0.5)
        
        m = MAX(dl, ind)
        npol[ind] = npol[ind] + 1
      ENDFOR
      
      ; Now sort out the order of npol. Starts from innermost xpoint
      ; and goes clockwise around: PF, SOL, PF, PF, SOL, PF
      
      npol2 = npol
      xpt = si[0] ; X-point index to start with
      FOR i=0, critical.n_xpoint-1 DO BEGIN
        npol2[3*i] = npol[xpt]  ; PF part 0
        
        ; Get the SOL ID
        solid = (*pf_info[xpt]).sol[0]
        IF (*sol_info[solid]).xpt1 NE xpt THEN BEGIN
          PRINT, "ERROR: Indexing inconsistent => Bug!"
          STOP
        ENDIF
        npol2[3*i+1] = npol[2*critical.n_xpoint + solid]
        
        ; Get the next x-point
        xpt = (*sol_info[solid]).xpt2
        npol2[3*i+2] = npol[critical.n_xpoint + xpt] ; PF part 1
      ENDFOR
      
      npol = npol2
    ENDIF
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Poloidal spacing. Need to ensure regular spacing
    ; around x-points.

    ; Calculate distance for equal spacing in each region
    xpt = si[0] ; Start with the innermost x-point
    xpt_dist = FLTARR(critical.n_xpoint, 4) ; Distance between x-point and first grid point
    
    ; NOTE: xpt_dist indices go clockwise around the X-point, starting
    ; from the lower left leg when the core is at the top.
    ; For the lower x-point,
    ;     0 = inner leg, 1 = inner SOL, 2 = outer sol, 3 = outer leg
    ; For the upper x-point (if any)
    ;     0 = outer leg, 1 = outer sol, 2 = iner sol, 3 = inner leg
    
    FOR i=0, critical.n_xpoint-1 DO BEGIN
      ; Grid the lower PF region
      
      ; Calculate poloidal distance along starting line
      ; NOTE: (ri0, zi0) is the line going from x-point down the leg
      ; to the left, if the core is at the top.
      ; - If xpt is the lower x-point, then this corresponds to the
      ;   lower inner leg.
      ; - If xpt is the upper x-point then this is the upper outer leg
      
      poldist = line_dist(R, Z, (*pf_info[xpt]).ri0, (*pf_info[xpt]).zi0) ; Poloidal distance along line
      xdist = MAX(poldist) * 0.5 / FLOAT(npol[3*i]) ; Equal spacing
      
      xpt_dist[xpt, 0] = xdist
      
      ; SOL
      solid = (*pf_info[xpt]).sol[0]
      
      poldist = line_dist(R, Z, (*sol_info[solid]).ri, (*sol_info[solid]).zi)
      xdist = MAX(poldist) * 0.5 / FLOAT(npol[3*i+1])

      PRINT, "S :", solid, max(poldist), npol[3*i+1]
      
      xpt2 = (*sol_info[solid]).xpt2
      
      xpt_dist[xpt, 1] = xdist
      xpt_dist[xpt2, 2] = xdist

      ; Second PF region
      xpt = xpt2
      
      poldist = line_dist(R, Z, (*pf_info[xpt]).ri1, (*pf_info[xpt]).zi1)
      xdist = MAX(poldist) * 0.5 / FLOAT(npol[3*i+2])
      
      xpt_dist[xpt, 3] = xdist
    ENDFOR
    ;FOR i=0, critical.n_xpoint-1 DO xpt_dist[i,*] = MEAN(xpt_dist[i,*])
    
    FOR i=0, critical.n_xpoint-1 DO BEGIN
      md = MAX(xpt_dist[i,*])
      xpt_dist[i,*] = 0.5*xpt_dist[i,*] + 0.5*md
    ENDFOR
    
    ; Try to equalise
    FOR i=0, critical.n_xpoint-1 DO BEGIN
      pf_f = faxis + fnorm*pf_psi_vals[i,0,0]
      sol_in_f = faxis + fnorm*pf_psi_vals[i,0,TOTAL(nrad,/int)-1]
      sol_out_f = faxis + fnorm*pf_psi_vals[i,1,TOTAL(nrad,/int)-1]
      core_f = faxis + fnorm*psi_inner[0]
      PRINT, "PSI", psi_inner[0], core_f
      
      IF KEYWORD_SET(xpt_mindist) THEN xpt_dist[i,*] = xpt_dist[i,*] > xpt_mindist
    ENDFOR

    ; Each x-point now has a distance. Could multiply by a scaling
    ; factor to adjust x-point resolution
    IF KEYWORD_SET(xpt_mul) THEN BEGIN
      PRINT, "Multiplying x-point distance by ", xpt_mul
      xpt_dist = xpt_dist * xpt_mul
    ENDIF

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Calculate distances along starting line
    xpt = si[0] ; Start with the innermost x-point
    FOR i=0, critical.n_xpoint-1 DO BEGIN
      
      solid = (*pf_info[xpt]).sol[0]
      xpt2 = (*sol_info[solid]).xpt2 ; the next x-point
      
      ydown_dist = xpt_dist[xpt, 1]
      ; Locate this point on the separatrix (core2)
      dist = line_dist(R, Z, (*sep_info[xpt]).core2_ri, (*sep_info[xpt]).core2_zi); Distance along the separatrix
      sepi = INTERPOL(findgen(N_ELEMENTS(dist)), dist, ydown_dist) ; Index into separatrix
      ; Follow from sep_info[i]->core2 to just inside starting f
      line = get_line(interp_data, R, Z, $
                      INTERPOLATE((*sep_info[xpt]).core2_ri, sepi), $
                      INTERPOLATE((*sep_info[xpt]).core2_zi, sepi), $
                      0.95*f_cont + 0.05*faxis, npt=30)
      ; Find intersection of this line with starting line
      cpos = line_crossings((*sol_info[solid]).ri, (*sol_info[solid]).zi, 0, $
                            line[*,0], line[*,1], 0, ncross=ncross, inds1=start_ind)

      IF ncross NE 1 THEN BEGIN
        PRINT, "WARNING: PROBLEM MAPPING STARTING LOCATION"
      ENDIF ELSE BEGIN
        start_ind = start_ind[0] ; Got index into the starting line
        ; Find out distance along starting line
        dist = line_dist(R, Z, (*sol_info[solid]).ri, (*sol_info[solid]).zi)
        d = INTERPOLATE(dist, start_ind)
        ydown_dist = MIN([d, dist[N_ELEMENTS(dist)-1] - d])
      ENDELSE

      xpt_dist[xpt, 1] = ydown_dist
      
      ; Repeat for yup
      yup_dist = xpt_dist[xpt2, 2]
      ; Locate this point on the separatrix (core1)
      dist = line_dist(R, Z, (*sep_info[xpt2]).core1_ri, (*sep_info[xpt2]).core1_zi); Distance along the separatrix
      sepi = INTERPOL(findgen(N_ELEMENTS(dist)), dist, yup_dist) ; Index into separatrix
      ; Follow from sep_info[i]->core1 to just inside starting f
      line = get_line(interp_data, R, Z, $
                      INTERPOLATE((*sep_info[xpt2]).core1_ri, sepi), $
                      INTERPOLATE((*sep_info[xpt2]).core1_zi, sepi), $
                      0.95*f_cont + 0.05*faxis, npt=20)
      ; Find intersection of this line with starting line
      cpos = line_crossings((*sol_info[solid]).ri, (*sol_info[solid]).zi, 0, $
                            line[*,0], line[*,1], 0, ncross=ncross, inds1=start_ind)
      IF ncross NE 1 THEN BEGIN
        PRINT, "WARNING: PROBLEM MAPPING STARTING LOCATION"
      ENDIF ELSE BEGIN
        start_ind = start_ind[0] ; Got index into the starting line
        ; Find out distance along starting line
        dist = line_dist(R, Z, (*sol_info[solid]).ri, (*sol_info[solid]).zi)
        yup_dist = MAX(dist) - INTERPOLATE(dist, start_ind)
      ENDELSE
      
      xpt_dist[xpt2, 2] = yup_dist

      xpt = xpt2 ; Move to the next x-point
    ENDFOR
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Now have all the radial locations for each region
    ; and number of poloidal points
    ; 
    ; => Grid each section
    ; 
    
    IF KEYWORD_SET(bndryi) THEN BEGIN
      IF KEYWORD_SET(strictbndry) THEN BEGIN
        ; No part of the grid is allowed to be outside boundary
        PRINT, "Keeping the grid stictly within the boundary"
        gridbndry = bndryi
      ENDIF ELSE BEGIN
        ; Grid can leave boundary
        gridbndry = FLTARR(2,4)
        gridbndry[0,*] = [0, 0, nx-1, nx-1]
        gridbndry[1,*] = [0, ny-1, ny-1, 0]
      ENDELSE
    ENDIF
    
    ; Create 2D arrays for the grid
    Rxy = FLTARR(TOTAL(nrad,/int), TOTAL(npol,/int))
    Zxy = Rxy
    Rixy = Rxy
    Zixy = Rxy
    Psixy = Rxy
    
    ; Create arrays for the topology connections
    
    yup_xsplit   = LONARR(nnpol) ; X index where domain is split on the upper Y side
    ydown_xsplit = LONARR(nnpol) ; X index " " " lower Y side
    yup_xin      = LONARR(nnpol) ; The region number connected to on inner X, upper Y
    yup_xout     = LONARR(nnpol)
    ydown_xin    = LONARR(nnpol)
    ydown_xout   = LONARR(nnpol)
    
    xpt = si[0] ; Start with the innermost x-point
    ypos = 0
    rerun = 0   ; Flag. If 1 then have to re-run the grid generator
    FOR i=0, critical.n_xpoint-1 DO BEGIN
      ; Calculate maximum psi of x-point
      xpt_psi_max = MAX([xpt_psi[ci[i]], xpt_psi[ci[(i+1) MOD critical.n_xpoint]]])
      PRINT, "Gridding regions "+STR(3*i)+" to " +STR(3*i+2)
      ; Grid the lower PF region
      PRINT, "   x-point index ", xpt
      a = grid_region(interp_data, R, Z, $
                      (*pf_info[xpt]).ri0, (*pf_info[xpt]).zi0, $
                      faxis + fnorm*pf_psi_vals[xpt,0,*], $
                      (*pf_info[xpt]).npf-1, $
                      npol[3*i], $
                      sfirst=sfirst1, $
                      slast=slast1, $
                      boundary=gridbndry, $
                      ffirst=ffirst, flast=flast1, fpsi=fpsi, yup_dist=xpt_dist[xpt, 0], /oplot)
      Rxy[*, ypos:(ypos+npol[3*i]-1)] = a.Rxy
      Zxy[*, ypos:(ypos+npol[3*i]-1)] = a.Zxy
      Rixy[*, ypos:(ypos+npol[3*i]-1)] = a.Rixy
      Zixy[*, ypos:(ypos+npol[3*i]-1)] = a.Zixy
      FOR j=ypos, ypos+npol[3*i]-1 DO Psixy[*, j] = pf_psi_vals[xpt,0,*]
      ypos = ypos + npol[3*i]
      
      ; Set topology
      ydown_xsplit[3*i] = (*pf_info[xpt]).npf
      yup_xsplit[3*i]   = (*pf_info[xpt]).npf
      ydown_xin[3*i]    = -1 ; Target plates
      ydown_xout[3*i]   = -1 
      yup_xin[3*i]  = 3*((i-1+critical.n_xpoint) MOD critical.n_xpoint) + 2 ; Previous PF
      yup_xout[3*i] = 3*i+1 ; The following SOL region
      
      IF (flast1 - faxis)/fnorm LT xpt_psi_max THEN BEGIN
        PRINT, "WARNING: Due to intersections with the boundary,"
        PRINT, "         the SOL can't cover both x-points"
        IF KEYWORD_SET(strictbndry) THEN BEGIN
          PRINT, "** Switching off strict boundary"
          strictbndry = 0
        ENDIF ELSE IF rerun EQ 0 THEN BEGIN
          ; Only if haven't already marked as restarting
          PRINT, "** Removing the boundary"
          boundary = 0
        ENDIF
        
        flast1 = faxis + fnorm*MAX(pf_psi_vals[xpt,0,*]) ; don't restrict range
        rerun = 1
      ENDIF

      plot_grid_section, a, color=4

      IF sfirst1 GT 0 THEN BEGIN
        PRINT, "  PF region "+STR(3*i)+" incomplete"
        PRINT, "    => Inner only good to f = "+STR(ffirst)
        pfirst = (ffirst - faxis)/fnorm
        PRINT, "    => Inner normalised psi = "+STR(pfirst)
        
        ; Change the inner psi
        w = WHERE(ci EQ xpt)
        IF pfirst GT psi_inner[w[0]+1] THEN BEGIN
          psi_inner[w[0]+1] = pfirst
        ENDIF
        rerun = 1  ; Signal that the grid needs to be rerun
      ENDIF

      ; SOL region
      solid = (*pf_info[xpt]).sol[0]
      
      PRINT, "   SOL index ", solid

      xpt2 = (*sol_info[solid]).xpt2 ; the second x-point

      ; Start and end positions.
      ydown_dist = xpt_dist[xpt, 1]
      yup_dist = xpt_dist[xpt2, 2]
      ; Grid spacing
      ydown_space = MAX([xpt_dist[xpt, 1], xpt_dist[xpt, 2]]) ;0.5*(xpt_dist[xpt, 1] + xpt_dist[xpt, 2])
      yup_space   = MAX([xpt_dist[xpt2, 1], xpt_dist[xpt2, 2]]) ;0.5*(xpt_dist[xpt2, 1] + xpt_dist[xpt2, 2])
      
      a = grid_region(interp_data, R, Z, $
                      (*sol_info[solid]).ri, (*sol_info[solid]).zi, $
                      faxis + fnorm*sol_psi_vals[solid,*], $
                      nrad[0]-1, $
                      npol[3*i+1], $
                      sfirst=sfirst2, $
                      slast=slast2, $
                      boundary=gridbndry, $
                      ffirst=ffirst, flast=flast2, fpsi=fpsi, $
                      ydown_dist=ydown_dist, yup_dist=yup_dist, $
                      ydown_space=ydown_space, yup_space=yup_space, /oplot)
      
      Rxy[*, ypos:(ypos+npol[3*i+1]-1)] = a.Rxy
      Zxy[*, ypos:(ypos+npol[3*i+1]-1)] = a.Zxy
      Rixy[*, ypos:(ypos+npol[3*i+1]-1)] = a.Rixy
      Zixy[*, ypos:(ypos+npol[3*i+1]-1)] = a.Zixy
      FOR j=ypos, ypos+npol[3*i+1]-1 DO Psixy[*, j] = sol_psi_vals[solid,*]
      ypos = ypos + npol[3*i+1]

      ydown_xsplit[3*i+1] = (*pf_info[xpt]).npf
      yup_xsplit[3*i+1]   = (*pf_info[(*sol_info[solid]).xpt2]).npf
      ydown_xin[3*i+1]    = 3*( (i-1+critical.n_xpoint) MOD critical.n_xpoint )+1 ; Previous SOL region
      ydown_xout[3*i+1]   = 3*i ; Previous PF region
      yup_xin[3*i+1]  = 3*( (i+1) MOD critical.n_xpoint )+1 ; Next SOL region
      yup_xout[3*i+1] = 3*i + 2 ; Next PF region

      IF (flast2 - faxis)/fnorm LT xpt_psi_max THEN BEGIN
        PRINT, "WARNING: Due to intersections with the boundary,"
        PRINT, "         the SOL can't cover both x-points"
        IF KEYWORD_SET(strictbndry) THEN BEGIN
          PRINT, "** Switching off strict boundary"
          strictbndry = 0
        ENDIF ELSE IF rerun EQ 0 THEN BEGIN
          PRINT, "** Removing the boundary"
          boundary = 0
        ENDIF
        
        flast2 = faxis + fnorm*MAX(sol_psi_vals[solid,*]) ; don't restrict range
        rerun = 1
      ENDIF

      plot_grid_section, a, color=4
      
      ; Second PF region
      xpt = (*sol_info[solid]).xpt2
      PRINT, "   x-point index ", xpt
      a = grid_region(interp_data, R, Z, $
                      (*pf_info[xpt]).ri1, (*pf_info[xpt]).zi1, $
                      faxis + fnorm*pf_psi_vals[xpt,1,*], $
                      (*pf_info[xpt]).npf-1, $
                      npol[3*i+2], $
                      sfirst=sfirst3, $
                      slast=slast3, $
                      boundary=gridbndry, $
                      ffirst=ffirst, flast=flast3, fpsi=fpsi, $
                      ydown_dist=xpt_dist[xpt, 3], /oplot)
      Rxy[*, ypos:(ypos+npol[3*i+2]-1)] = a.Rxy
      Zxy[*, ypos:(ypos+npol[3*i+2]-1)] = a.Zxy
      Rixy[*, ypos:(ypos+npol[3*i+2]-1)] = a.Rixy
      Zixy[*, ypos:(ypos+npol[3*i+2]-1)] = a.Zixy
      FOR j=ypos, ypos+npol[3*i+2]-1 DO Psixy[*, j] = pf_psi_vals[xpt,1,*]
      ypos = ypos + npol[3*i+2]

      ; Set topology
      ydown_xsplit[3*i+2] = (*pf_info[xpt]).npf
      yup_xsplit[3*i+2]   = (*pf_info[xpt]).npf
      ydown_xin[3*i+2]    = 3*( (i+1) MOD critical.n_xpoint ) ; Next PF region
      ydown_xout[3*i+2]   = 3*i+1 ; Previous SOL region
      yup_xin[3*i+2]  = -1 ; Target plates
      yup_xout[3*i+2] = -1

      IF (flast3 - faxis)/fnorm LT xpt_psi_max THEN BEGIN
        PRINT, "WARNING: Due to intersections with the boundary,"
        PRINT, "         the SOL can't cover both x-points"
        IF KEYWORD_SET(strictbndry) THEN BEGIN
          PRINT, "** Switching off strict boundary"
          strictbndry = 0
        ENDIF ELSE IF rerun EQ 0 THEN BEGIN
          PRINT, "** Removing the boundary"
          boundary = 0
        ENDIF
        
        flast2 = faxis + fnorm*MAX(pf_psi_vals[xpt,1,*]) ; don't restrict range
        rerun = 1
      ENDIF

      plot_grid_section, a, color=4
      
      IF sfirst3 GT 0 THEN BEGIN
        PRINT, "  PF region "+STR(3*i+2)+" incomplete"
        PRINT, "    => Inner only good to f = "+STR(ffirst)
        pfirst = (ffirst - faxis)/fnorm
        PRINT, "    => Inner normalised psi = "+STR(pfirst)
        
        ; Change the inner psi. Make sure hasn't already been reduced
        w = WHERE(ci EQ xpt)
        IF pfirst GT psi_inner[w[0]+1] THEN BEGIN
          psi_inner[w[0]+1] = pfirst
        ENDIF
        rerun = 1  ; Signal that the grid needs to be rerun
      ENDIF
      
      ; Check if any of the outer edges failed
      slast = MIN([slast1, slast2, slast3])
      IF slast LT (TOTAL(nrad,/int)-1) THEN BEGIN
        plast = MIN(([flast1, flast2, flast3] - faxis) / fnorm)
        flast = faxis + fnorm*plast
        PRINT, "   => These regions only good until f = "+STR(flast)
        PRINT, "   => i.e. Normalised psi of "+STR( plast )
        
        ; Set the new outer psi for these regions
        w = WHERE(ci EQ primary_xpt)
        id = (solid - w[0] + critical.n_xpoint) MOD critical.n_xpoint
        psi_outer[id] = plast
        PRINT, "SETTING PSI_OUT", id, solid
        rerun = 1
      ENDIF
    ENDFOR

    IF rerun THEN BEGIN
      ; Create new settings
      
      IF nrad_flexible THEN nrad = TOTAL(nrad, /integer)  ; Allow nrad to change again

      new_settings = {psi_inner:psi_inner, psi_outer:psi_outer, $ ; New ranges
                      nrad:nrad, npol:npol, $
                      rad_peaking:settings.rad_peaking, pol_peaking:settings.pol_peaking}
      
      IF iter GT 0 THEN strictbndry = 0

      PRINT, "** Re-running grid generator with changed settings"
      PRINT, "psi outer = ", psi_outer
      RETURN, create_grid(F, R, Z, new_settings, critical=critical, $
                          boundary=boundary, strictbndry=strictbndry, $
                          iter=iter+1, nrad_flexible=nrad_flexible, $
                          single_rad_grid=single_rad_grid, fast=fast)
      
    ENDIF
    
    ; Successfully created grid
    
    new_settings = {psi_inner:psi_inner, psi_outer:psi_outer, $ ; New ranges
                    nrad:nrad, npol:npol, $
                    rad_peaking:settings.rad_peaking, pol_peaking:settings.pol_peaking}
    
    ; Calculate magnetic field components
    dpsidR = FLTARR(TOTAL(nrad, /int), TOTAL(npol, /int))
    dpsidZ = dpsidR

    interp_data.method = 2

    FOR i=0,TOTAL(nrad,/int)-1 DO BEGIN
      FOR j=0,TOTAL(npol,/int)-1 DO BEGIN
        local_gradient, interp_data, Rixy[i,j], Zixy[i,j], status=status, $
          dfdr=dfdr, dfdz=dfdz
        ; dfd* are derivatives wrt the indices. Need to multiply by dr/di etc
        dpsidR[i,j] = dfdr/INTERPOLATE(DERIV(R),Rixy[i,j]) 
        dpsidZ[i,j] = dfdz/INTERPOLATE(DERIV(Z),Zixy[i,j]) 
      ENDFOR
    ENDFOR

    result = {error:0, $ ; Signals success
              psi_inner:psi_inner, psi_outer:psi_outer, $ ; Range of psi
              nrad:nrad, npol:npol, $  ; Number of points in each domain
              Rixy:Rixy, Zixy:Zixy, $  ; Indices into R and Z of each point
              Rxy:Rxy, Zxy:Zxy, $ ; Location of each grid point
              psixy:psixy, $ ; Normalised psi for each point
              dpsidR:dpsidR, dpsidZ:dpsidZ, $ ; Psi derivatives (for Bpol)
              faxis:faxis, fnorm:fnorm, $ ; Psi normalisation factors
              settings:new_settings, $ ; Settings used to create grid
              critical:critical, $ ; Critical points
              yup_xsplit:yup_xsplit, $ ; X index where domain splits (number of points in xin)
              ydown_xsplit:ydown_xsplit, $
              yup_xin:yup_xin, yup_xout:yup_xout, $ ; Domain index to connect to
              ydown_xin:ydown_xin, ydown_xout:ydown_xout}
              
              
  ENDELSE

  CATCH, /cancel

  RETURN, result
END
