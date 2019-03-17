;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Separatrix line routine
; 
; Finds lines going through a given x-point. Uses this to
; return lines between the x-point and boundary
;
; INPUTS
;
; interp_data   Input to local_gradient, describes interpolation
; R, Z          1D arrays of major radius, height [m]
; xpt_ri, xpt_zi   Indices of the X-point into R and Z arrays
; opt_ri, opt_zi   Indices of the primary O-point into R and Z arrays
; 
; OPTIONAL INPUTS
; 
; boundary (in)  Optional 2D array [2,n] of points on the boundary
;
; OUTPUTS
; 

FUNCTION leg_separatrix2, interp_data, R, Z, xpt_ri, xpt_zi, $
                          opt_ri, opt_zi, $ ; Location of primary O-point
                          status=status, $
                          boundary=boundary, $
                          debug=debug, xpt_psi=xpt_psi
  
  psi = interp_data.f
  nr = interp_data.nx
  nz = interp_data.ny

  IF NOT KEYWORD_SET(boundary) THEN BEGIN
    bndry = DBLARR(2,4)
    bndry[0,*] = [1, nr-2, nr-2, 1]
    bndry[1,*] = [1, 1, nz-2, nz-2]
  ENDIF ELSE bndry = boundary
  
  drpdi = R[1] - R[0]
  dzpdi = Z[1] - Z[0]

  ; Get value of psi at the X-point
  
  IF KEYWORD_SET(xpt_psi) THEN f0 = xpt_psi ELSE BEGIN
    local_gradient, interp_data, xpt_ri, xpt_zi, f=f0
  ENDELSE
  
  ; Get contour lines at this level, i.e. the separatrix

  contour_lines, psi, levels=[f0], path_info=info, path_xy=xy
  nsep = N_ELEMENTS(info) ; Will be split into two or more lines
  
  ncore = 0
  npf = 0

  ; Find where the separatrix intersects a circle around the X-point
  FOR i=0, nsep-1 DO BEGIN
    sep_ri = REFORM(xy[0,info[i].offset:(info[i].offset+info[i].n-1)])
    sep_zi = REFORM(xy[1,info[i].offset:(info[i].offset+info[i].n-1)])
    
    ; Find smallest distance to x-point
    md = SQRT( MIN((sep_ri - xpt_ri)^2 + (sep_zi - xpt_zi)^2) )
    
    ; Create a circle around the x-point
    di = 2 ; Radius of 2 grid points
    IF (md GT di) AND (md LT 6) THEN di = md * 1.25D
    dthe = 2.D*!DPI*FINDGEN(6)/6.D
    ri = xpt_ri + di*COS(dthe)
    zi = xpt_zi + di*SIN(dthe)
    
    IF KEYWORD_SET(debug) THEN BEGIN
      plot, sep_ri, sep_zi
      oplot, ri, zi, color=2
      
      oplot, [xpt_ri], [xpt_zi], psym=1, color=4
      
      oplot, bndry[0,*], bndry[1,*], color=2, thick=2
    ENDIF

    cpos = line_crossings(sep_ri, sep_zi, info[i].type, $
                          ri, zi, 1, $
                          ncross=ncross, inds1=inds)

    PRINT, "Intersections: ", ncross
    FOR j = 0, ncross-1 DO BEGIN
      ; Intersection. Get location
      cri = INTERPOLATE(sep_ri, inds[j], /DOUBLE)
      czi = INTERPOLATE(sep_zi, inds[j], /DOUBLE)
      
      IF KEYWORD_SET(debug) THEN oplot, [cri], [czi], psym=2

      ;Get direction of line
      ; DERIV uses one-sided differences at the ends of the array. These may be
      ; very inaccurate because the separatrix points are not evenly spaced.
      ; Therefore extend the ends of the array by cycling around a couple of
      ; points from the other end. This is incorrect for the leg contours,
      ; which are open, but since they are open they should not have ends near
      ; the X-point, so this should not affect them.
      ; Need to add +2 to inds[j] because we added two points onto the
      ; beginning of sep_ri/sep_zi.
      drdi = INTERPOLATE(DERIV([sep_ri[-2:*],sep_ri,sep_ri[0:2]]), inds[j]+2, /DOUBLE)
      dzdi = INTERPOLATE(DERIV([sep_zi[-2:*],sep_zi,sep_zi[0:2]]), inds[j]+2, /DOUBLE)
      
      ; First check if this is towards or away from the X-point
      dir = 1 ; direction to go away from x-point
      
      d = drdi*(cri - xpt_ri) + dzdi*(czi - xpt_zi) ; Dot-product
      IF d LT 0.D THEN dir = -1
      
      ; Get the indices for a line radiating from the x-point
      si = [inds[j]]
      IF dir GT 0 THEN BEGIN
        in = CEIL(si[0])
        IF in LT N_ELEMENTS(sep_ri) THEN BEGIN
          ; normal case
          si = [si, in + indgen(N_ELEMENTS(sep_ri) - in), indgen(in - 1)]
        ENDIF ELSE BEGIN
          ; can't call indgen(0) so handle this specially
          si = [si, indgen(in - 1)]
        ENDELSE
      ENDIF ELSE BEGIN
        in = FLOOR(si[0])
        ; contour is closed, so we can loop around: start at si, add elements
        ; until the beginning of the contour, then add elements starting from the
        ; end of the contour
        IF in LT N_ELEMENTS(sep_ri)-1 THEN BEGIN
          ; normal case
          si = [si, reverse(indgen(in+1)), reverse(indgen(N_ELEMENTS(sep_ri)-in-1)) + in + 1]
        ENDIF ELSE BEGIN
          ; can't call indgen(0) so handle this specially
          si = [si, reverse(indgen(in+1))]
        ENDELSE
      ENDELSE
      sepri = INTERPOLATE(sep_ri, si, /DOUBLE)
      sepzi = INTERPOLATE(sep_zi, si, /DOUBLE)
      
      ; Then check if this is towards the O-point (core) or away (PF)
      
      d = drdi*dir * (opt_ri - xpt_ri)*drpdi^2 + dzdi*dir * (opt_zi - xpt_zi) * dzpdi^2
      ;OPLOT, INTERPOLATE(R, xpt_ri, /DOUBLE) + [0.D, drdi*dir*drpdi * 100], INTERPOLATE(Z, xpt_zi, /DOUBLE) + [0.D, dzdi*dir*dzpdi * 100], color=2
      
      IF d GT 0 THEN BEGIN
        ; Core
        
        ; Find where both ri and zi have had an extrema.
        ; This will take the contour at most half way around the closed part of
        ; the separatrix.
        ; Taking the first extremum in ri or zi (old behaviour) risks having a
        ; very short contour on one leg if the X-point is not quite at the
        ; bottom or top of the core region
        n = N_ELEMENTS(sepri)
        dr = DERIV(sepri)
        dr = dr[1:*] * dr[0:(n-2)]
        dz = DERIV(sepzi)
        dz = dz[1:*] * dz[0:(n-2)]

        inr = MIN(WHERE(dr[1:*] LE 0.0D)) + 1
        inz = MIN(WHERE(dz[1:*] LE 0.0D)) + 1
        in = MAX([inr, inz])
        
        if in GT 0 THEN BEGIN
          ; if in<=0 then there is no extremum, so don't truncate sepri/sepzi
          sepri = sepri[0:in]
          sepzi = sepzi[0:in]
        ENDIF
        
        IF KEYWORD_SET(debug) THEN OPLOT, sepri, sepzi, color=4
        
        IF ncore EQ 0 THEN BEGIN
          core1 = [[sepri], [sepzi]]
          ncore = 1
        ENDIF ELSE core2 = [[sepri], [sepzi]]
        
      ENDIF ELSE BEGIN
        ; PF
        
        sepri = INTERPOLATE(sep_ri, si, /DOUBLE)
        sepzi = INTERPOLATE(sep_zi, si, /DOUBLE)
        
        ; Find where it crosses a boundary
        
        cpos = line_crossings(sepri, sepzi, 0, $
                              bndry[0,*], bndry[1,*], 1, $
                              ncross=ncross, inds1=in)

        ; Find where it crosses the edge of the grid
        cpos = line_crossings(sepri, sepzi, 0, $
                              [0, 0, nr, nr], [0, nz, nz, 0], 1, $
                              ncross=ncrossgrid, inds1=ingrid)
        
        IF ncross GT 0 THEN BEGIN
          ; Keep points past the boundary in case we want to include y-boundary
          ; guard cells
          IF ncrossgrid GT 0 THEN BEGIN
            ; never go further than the edge of the grid
            lastgridind = FLOOR(ingrid[0])
          ENDIF ELSE BEGIN
            lastgridind = N_ELEMENTS(sepri)-1
          ENDELSE

          lastind = FLOOR(in[0])
          sepri = [sepri[0:lastind], INTERPOLATE(sepri, in[0], /DOUBLE), sepri[lastind+1:lastgridind]]
          sepzi = [sepzi[0:lastind], INTERPOLATE(sepzi, in[0], /DOUBLE), sepzi[lastind+1:lastgridind]]
        ENDIF ELSE BEGIN
          lastind = N_ELEMENTS(sepri)-1
        ENDELSE
        
        IF KEYWORD_SET(debug) THEN OPLOT, sepri, sepzi, color=3
        
        IF npf EQ 0 THEN BEGIN
          pf1 = [[sepri], [sepzi]]
          pf1_lastind = lastind+1
          npf = 1
        ENDIF ELSE BEGIN
          pf2 = [[sepri], [sepzi]]
          pf2_lastind = lastind+1
        ENDELSE
        
      ENDELSE
    ENDFOR
  ENDFOR

  ; Check core1 and pf1 are aligned
  IF ABS( (core1[0,0] - xpt_ri)*(pf2[0,0] - xpt_ri) + (core1[0,1] - xpt_zi)*(pf2[0,1] - xpt_zi) ) GT ABS( (core2[0,0] - xpt_ri)*(pf2[0,0] - xpt_ri) + (core2[0,1] - xpt_zi)*(pf2[0,1] - xpt_zi) ) THEN BEGIN
    ; swap core1 and core2
    tmp = core1
    core1 = core2
    core2 = tmp
  ENDIF
  
  IF KEYWORD_SET(debug) THEN BEGIN
    STOP
  ENDIF ELSE BEGIN
    OPLOT, INTERPOLATE(R, pf1[0:pf1_lastind,0], /DOUBLE), INTERPOLATE(Z, pf1[0:pf1_lastind,1], /DOUBLE), color=3, thick=2
    OPLOT, INTERPOLATE(R, pf2[0:pf2_lastind,0], /DOUBLE), INTERPOLATE(Z, pf2[0:pf2_lastind,1], /DOUBLE), color=4, thick=2
    
    OPLOT, INTERPOLATE(R, core1[*,0], /DOUBLE), INTERPOLATE(Z, core1[*,1], /DOUBLE), color=3, thick=2
    OPLOT, INTERPOLATE(R, core2[*,0], /DOUBLE), INTERPOLATE(Z, core2[*,1], /DOUBLE), color=4, thick=2
  ENDELSE

  RETURN, {leg1:pf1, leg1_lastind:pf1_lastind, leg2:pf2, leg2_lastind:pf2_lastind, core1:core1, core2:core2, ri:xpt_ri, zi:xpt_zi}
END
