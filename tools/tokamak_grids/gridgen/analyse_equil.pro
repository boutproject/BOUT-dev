;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Equilibrium analysis routine
; 
; Takes a RZ psi grid, and finds x-points and o-points
; 
; F - F(nr, nz) 2D array of psi values
; R - R(nr) 1D array of major radii
; Z - Z(nz) 1D array of heights
;
; Returns a structure of critical points containing:
;
;   n_opoint, n_xpoint   - Number of O- and X-points
;   primary_opt          - Index of plasma centre O-point
;   inner_sep            - X-point index of inner separatrix
;   opt_ri, opt_zi       - R and Z indices for each O-point
;   opt_f                - Psi value at each O-point
;   xpt_ri, xpt_zi       - R and Z indices for each X-point
;   xpt_f                - Psi value of each X-point
; 

FUNCTION remove_ind, var, ind
  n = N_ELEMENTS(var)
  IF ind EQ 0 THEN BEGIN
    RETURN, var[1:*]
  ENDIF ELSE IF ind EQ (n-1) THEN BEGIN
    RETURN, var[0:(n-2)]
  ENDIF
  RETURN, [var[0:(ind-1)], var[(ind+1):*]]
END

FUNCTION analyse_equil, F, R, Z
  s = SIZE(F, /DIMENSION)
  nx = s[0]
  ny = s[1]
  
  ;;;;;;;;;;;;;;;; Find critical points ;;;;;;;;;;;;;
  ;
  ; Need to find starting locations for O-points (minima/maxima)
  ; and X-points (saddle points)
  ;
  
  dfdr = DBLARR(nx, ny)
  FOR j=0, ny-1 DO BEGIN
    dfdr[*,j] = diff(f[*,j])
  ENDFOR
  
  dfdz = DBLARR(nx, ny)
  FOR i=0, nx-1 DO BEGIN
    dfdz[i,*] = diff(f[i,*])
  ENDFOR

  ; Use contour to get crossing-points where dfdr = dfdz = 0
  
  contour_lines, dfdr, findgen(nx), findgen(ny), levels=[0.0D], $
    path_info=rinfo, path_xy=rxy
  
  contour_lines, dfdz, findgen(nx), findgen(ny), levels=[0.0D], $
    path_info=zinfo, path_xy=zxy
  
  ; Find where these two cross
  nextrema = 0
  
  FOR i=0, N_ELEMENTS(rinfo)-1 DO BEGIN
    FOR j=0, N_ELEMENTS(zinfo)-1 DO BEGIN
      cross = line_crossings(rxy[0,rinfo[i].offset:(rinfo[i].offset + rinfo[i].n - 1)], $
                             rxy[1,rinfo[i].offset:(rinfo[i].offset + rinfo[i].n - 1)], $
                             rinfo[i].type, $
                             zxy[0,zinfo[j].offset:(zinfo[j].offset + zinfo[j].n - 1)], $
                             zxy[1,zinfo[j].offset:(zinfo[j].offset + zinfo[j].n - 1)], $
                             zinfo[j].type, $
                             ncross=ncross)
      IF ncross GT 0 THEN BEGIN
        IF nextrema EQ 0 THEN BEGIN
          rex = REFORM(cross[0,*])
          zex = REFORM(cross[1,*])
        ENDIF ELSE BEGIN
          rex = [rex, REFORM(cross[0,*])]
          zex = [zex, REFORM(cross[1,*])]
        ENDELSE
        nextrema = nextrema + ncross
      ENDIF
    ENDFOR
  ENDFOR
  
  ; Check for points too close to the edges
  w = WHERE((rex GT 2) AND (rex LT nx-3) AND $
            (zex GT 2) AND (zex LT ny-3), nextrema)
  
  rex = rex[w]
  zex = zex[w]
  
  rex = ROUND(rex)
  zex = ROUND(zex)
  
  PRINT, ""
  PRINT, "Number of critical points:", nextrema
  
  ;;;;;;;;;;;;;;; Characterise extrema ;;;;;;;;;;;;;;;;;
  ; Fit a surface through local points using 6x6 matrix
  ; This is to determine the type of extrema, and to
  ; refine the location
  ; 
  
  n_opoint = 0
  n_xpoint = 0

  ; Calculate inverse matrix
  rio = [-1, 0, 0, 0, 1, 1] ; R index offsets
  zio = [ 0,-1, 0, 1, 0, 1] ; Z index offsets
  
  ; Fitting a + br + cz + drz + er^2 + fz^2
  A = TRANSPOSE([[DBLARR(6)+1], $
                 [rio], $
                 [zio], $
                 [rio*zio], $
                 [rio^2], $
                 [zio^2]])
  SVDC, A,W,U,V

  FOR e=0, nextrema-1 DO BEGIN
    ; Fit in index space so result is index number
    
    PRINT, "Critical point "+STRTRIM(STRING(e),2)
    
    valid = 1

    localf = DBLARR(6)
    FOR i=0, 5 DO BEGIN
      ; Get the f value in a stencil around this point
      xi = ((rex[e]+rio[i]) > 0) < (nx-1) ; Zero-gradient at edges
      yi = ((zex[e]+zio[i]) > 0) < (ny-1)
      localf[i] = F[xi, yi]
    ENDFOR
    res=SVSOL(U,W,V,localf)
    
    ; Res now contains [a,b,c,d,e,f]
    ;                  [0,1,2,3,4,5]

    ; This determines whether saddle or extremum
    det = 4.D*res[4]*res[5] - res[3]^2
    
    IF det LT 0.0D THEN BEGIN
      PRINT, "   X-point"
    ENDIF ELSE BEGIN
      PRINT, "   O-point"
    ENDELSE
    
    ; Get location (2x2 matrix of coefficients)
    
    rinew = (res[3]*res[2] - 2.D*res[1]*res[5]) / det
    zinew = (res[3]*res[1] - 2.D*res[4]*res[2]) / det

    IF (ABS(rinew) GT 1.D) OR (ABS(zinew) GT 1.0D) THEN BEGIN
      ; Method has gone slightly wrong. Try a different method.
      ; Get a contour line starting at this point. Should
      ; produce a circle around the real o-point. 
      PRINT, "   Fitted location deviates too much"
      IF det LT 0.0D THEN BEGIN
        PRINT, "   => X-point probably not valid"
        PRINT, "      deviation = "+STR(rinew)+","+STR(zinew)
        ;valid = 0
      ENDIF ELSE BEGIN
        
        contour_lines, F, findgen(nx), findgen(ny), levels=[F[rex[e], zex[e]]], $
          path_info=info, path_xy=xy
        
        IF N_ELEMENTS(info) GT 1 THEN BEGIN
          ; More than one contour. Select the one closest
          ind = closest_line(info, xy, rex[e], zex[e])
          info = info[ind]
        ENDIF ELSE info = info[0]
        
        rinew = 0.5D*(MAX(xy[0, info.offset:(info.offset + info.n - 1)]) + $
                     MIN(xy[0, info.offset:(info.offset + info.n - 1)])) - rex[e]
        zinew = 0.5D*(MAX(xy[1, info.offset:(info.offset + info.n - 1)]) + $
                     MIN(xy[1, info.offset:(info.offset + info.n - 1)])) - zex[e]
        
        IF (ABS(rinew) GT 2.D) OR (ABS(zinew) GT 2.0D) THEN BEGIN
          PRINT, "   Backup method also failed. Keeping initial guess"
          rinew = 0.D
          zinew = 0.D
        ENDIF
      ENDELSE
    ENDIF

    IF valid THEN BEGIN
      fnew = res[0] + res[1]*rinew + res[2]*zinew $
        + res[3]*rinew*zinew + res[4]*rinew^2 + res[5]*zinew^2
      
      rinew = rinew + rex[e]
      zinew = zinew + zex[e]
      
      PRINT, "   Starting index: " + STR(rex[e])+", "+STR(zex[e])
      PRINT, "   Refined  index: " + STR(rinew)+", "+STR(zinew)
      
      rnew = INTERPOLATE(R, rinew, /DOUBLE)
      znew = INTERPOLATE(Z, zinew, /DOUBLE)
      
      PRINT, "   Position: " + STR(rnew)+", "+STR(znew)
      PRINT, "   F = "+STR(fnew)
      
      IF det LT 0.0D THEN BEGIN
        
        IF n_xpoint EQ 0 THEN BEGIN
          xpt_ri = [rinew]
          xpt_zi = [zinew]
          xpt_f = [fnew]
          n_xpoint = n_xpoint + 1
        ENDIF ELSE BEGIN
          ; Check if this duplicates an existing point
          
          m = MIN((xpt_ri - rinew)^2 + (xpt_zi - zinew)^2, ind)
          IF m LT 2.D THEN BEGIN
            PRINT, "   Duplicates existing X-point."
          ENDIF ELSE BEGIN
            xpt_ri = [xpt_ri, rinew]
            xpt_zi = [xpt_zi, zinew]
            xpt_f = [xpt_f, fnew]
            n_xpoint = n_xpoint + 1
          ENDELSE
        ENDELSE
        
      ENDIF ELSE BEGIN
        
        IF n_opoint EQ 0 THEN BEGIN
          opt_ri = [rinew]
          opt_zi = [zinew]
          opt_f = [fnew]
          n_opoint = n_opoint + 1
        ENDIF ELSE BEGIN
          ; Check if this duplicates an existing point
          
          m = MIN((opt_ri - rinew)^2 + (opt_zi - zinew)^2, ind)
          IF m LT 2.D THEN BEGIN
            PRINT, "   Duplicates existing O-point"
          ENDIF ELSE BEGIN
            opt_ri = [opt_ri, rinew]
            opt_zi = [opt_zi, zinew]
            opt_f = [opt_f, fnew]
            n_opoint = n_opoint + 1
          ENDELSE
        ENDELSE
      ENDELSE
    ENDIF
  ENDFOR

  PRINT, "Number of O-points: "+STR(n_opoint)
  PRINT, "Number of X-points: "+STR(n_xpoint)

  IF n_opoint EQ 0 THEN BEGIN
    PRINT, "No O-points! Giving up on this equilibrium"
    RETURN, {n_opoint:0, n_xpoint:0, primary_opt:-1}
  ENDIF

  ;;;;;;;;;;;;;;; Find plasma centre ;;;;;;;;;;;;;;;;;;;
  ; Find the O-point closest to the middle of the grid
  dR = R[1] - R[0]
  dZ = Z[1] - Z[0]
  mind = dR^2 * (opt_ri[0] - (DOUBLE(nx)/2.D))^2 + dZ^2*(opt_zi[0] - (DOUBLE(ny)/2.D))^2
  ind = 0
  FOR i=1, n_opoint-1 DO BEGIN
    d = dR^2*(opt_ri[i] - (DOUBLE(nx)/2.D))^2 + dZ^2*(opt_zi[i] - (DOUBLE(ny)/2.D))^2
    IF d LT mind THEN BEGIN
      ind = i
      mind = d
    ENDIF
  ENDFOR
  
  primary_opt = ind
  PRINT, "Primary O-point is at "+STR(INTERPOLATE(R, opt_ri[ind], /DOUBLE)) + $
    ", " + STR(INTERPOLATE(Z, opt_zi[ind], /DOUBLE))
  PRINT, ""
  
  IF n_xpoint GT 0 THEN BEGIN
    
    ; Find the primary separatrix
    
    ; First remove non-monotonic separatrices
    nkeep = 0
    FOR i=0, n_xpoint-1 DO BEGIN
      ; Draw a line between the O-point and X-point
      
      n = 100 ; Number of points
      farr = DBLARR(n)
      dr = (xpt_ri[i] - opt_ri[primary_opt]) / DOUBLE(n)
      dz = (xpt_zi[i] - opt_zi[primary_opt]) / DOUBLE(n)
      FOR j=0, n-1 DO BEGIN
        ; interpolate f at this location
        farr[j] = INTERPOLATE(F, opt_ri[primary_opt] + dr*DOUBLE(j), opt_zi[primary_opt] + dz*DOUBLE(j), /DOUBLE)
      ENDFOR
      
      IF farr[n-1] LT farr[0] THEN BEGIN
        farr *= -1.0D ; Reverse, so maximum is always at the X-point
      ENDIF
      ; farr should be monotonic, and shouldn't cross any other separatrices
      
      ma = MAX(farr, maxind)
      mi = MIN(farr, minind)
      
      ; Discard if there is more than a 5% discrepancy in normalised
      ; psi between the maximum and the X-point, or the minimum and
      ; the O-point.
      IF (ma - farr[n-1])/(ma - farr[0]) GT 0.05D THEN continue
      IF (farr[0] - mi)/(farr[n-1] - mi) GT 0.05D THEN continue
      
      ; Monotonic, so add this to a list of x-points to keep
      IF nkeep EQ 0 THEN keep = [i] ELSE keep = [keep, i]
      nkeep = nkeep + 1
    ENDFOR
    
    IF nkeep GT 0 THEN BEGIN
      PRINT, "Keeping x-points ", keep
      xpt_ri = xpt_ri[keep]
      xpt_zi = xpt_zi[keep]
      xpt_f = xpt_f[keep]
    ENDIF ELSE PRINT, "No x-points kept"
    n_xpoint = nkeep

    ; Check for duplicates
    
    IF n_xpoint GT 1 THEN BEGIN
      i = 1
      REPEAT BEGIN
        ; note here MIN() sets the value of 'ind' to the index where the minimum was found
        m = MIN((xpt_ri[0:(i-1)] - xpt_ri[i])^2 + (xpt_zi[0:(i-1)] - xpt_zi[i])^2, ind)
        IF m LT 4.D THEN BEGIN
          PRINT, "Duplicates: ", i, ind
          
          IF ABS(opt_f[primary_opt] - xpt_f[i]) LT ABS(opt_f[primary_opt] - xpt_f[ind]) THEN BEGIN
            ; i is closer to O-point than ind.
            
            xpt_ri[ind] = xpt_ri[n_xpoint-1]
            xpt_zi[ind] = xpt_zi[n_xpoint-1]
            xpt_f[ind] = xpt_f[n_xpoint-1]
          ENDIF ELSE BEGIN
            xpt_ri[i] = xpt_ri[n_xpoint-1]
            xpt_zi[i] = xpt_zi[n_xpoint-1]
            xpt_f[i] = xpt_f[n_xpoint-1]
          ENDELSE
          xpt_ri = xpt_ri[0:(n_xpoint - 2)]
          xpt_zi = xpt_zi[0:(n_xpoint - 2)]
          xpt_f = xpt_f[0:(n_xpoint - 2)]
          n_xpoint = n_xpoint - 1
        ENDIF
        i = i + 1
      ENDREP UNTIL i GE n_xpoint
    ENDIF
    
    ; Now find x-point closest to primary O-point
    s = SORT(ABS(opt_f[primary_opt] - xpt_f))
    xpt_ri = xpt_ri[s]
    xpt_zi = xpt_zi[s]
    xpt_f = xpt_f[s]
    inner_sep = 0
    
  ENDIF ELSE BEGIN
    ; No x-points. Pick mid-point in f
   
    xpt_f = 0.5D*(MAX(F) + MIN(F))
    
    PRINT, "WARNING: No X-points. Setting separatrix to F = "+STR(xpt_f)

    xpt_ri = 0
    xpt_zi = 0
    inner_sep = 0
  ENDELSE

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Put results into a structure
  
  result = {n_opoint:n_opoint, n_xpoint:n_xpoint, $ ; Number of O- and X-points
            primary_opt:primary_opt, $ ; Which O-point is the plasma centre
            inner_sep:inner_sep, $ ; Innermost X-point separatrix
            opt_ri:opt_ri, opt_zi:opt_zi, opt_f:opt_f, $ ; O-point location (indices) and psi values
            xpt_ri:xpt_ri, xpt_zi:xpt_zi, xpt_f:xpt_f}  ; X-point locations and psi values
  
  RETURN, result
END
