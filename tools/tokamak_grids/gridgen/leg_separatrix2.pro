;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Separatrix line routine
; 
; Finds lines going through a given x-point. Uses this to
; return lines between the x-point and boundary
;
; INPUTS
; R, Z  
;
; OPTIONAL INPUTS
; 
; boundary (in)  Optional 2D array [2,n] of points on the boundary
;
; OUTPUTS
; 

FUNCTION leg_separatrix2, dctpsi, R, Z, xpt_ri, xpt_zi, $
                          opt_ri, opt_zi, $ ; Location of primary O-point
                          status=status, $
                          boundary=boundary, $
                          f=f, $
                          debug=debug
  
  IF NOT KEYWORD_SET(f) THEN BEGIN
    ; Need psi
    
    psi = DCT2D(dctpsi, /inverse)
  ENDIF ELSE psi = f

  s = SIZE(psi, /DIMENSION)
  nr = s[0]
  nz = s[1]
  
  ; Get value of psi at the X-point
  
  g = local_gradient(dctpsi, xpt_ri, xpt_zi, f=psi)
  f0 = g.f 
  
  ; Get contour lines at this level, i.e. the separatrix

  contour_lines, psi, levels=[f0], path_info=info, path_xy=xy
  nsep = N_ELEMENTS(info) ; Will be split into two or more lines
  
  ; Create a circle around the x-point
  
  di = 2 ; Radius of 2 grid points
  dthe = 2.*!PI*FINDGEN(6)/6.
  ri = xpt_ri + di*COS(dthe)
  zi = xpt_zi + di*SIN(dthe)

  ncore = 0
  npf = 0

  ; Find where the separatrix intersects the circle
  FOR i=0, nsep-1 DO BEGIN
    sep_ri = REFORM(xy[0,info[i].offset:(info[i].offset+info[i].n-1)])
    sep_zi = REFORM(xy[1,info[i].offset:(info[i].offset+info[i].n-1)])
    
    IF KEYWORD_SET(debug) THEN BEGIN
      plot, sep_ri, sep_zi
      oplot, ri, zi, color=2
      IF KEYWORD_SET(boundary) THEN BEGIN
        oplot, boundary[0,*], boundary[1,*], color=2, thick=2
      ENDIF
    ENDIF

    cpos = line_crossings(sep_ri, sep_zi, 0, $
                          ri, zi, 1, $
                          ncross=ncross, inds1=inds)
    PRINT, "Intersections: ", ncross
    FOR j = 0, ncross-1 DO BEGIN
      ; Intersection. Get location
      cri = INTERPOLATE(sep_ri, inds[j])
      czi = INTERPOLATE(sep_zi, inds[j])
      
      IF KEYWORD_SET(debug) THEN oplot, [cri], [czi], psym=2

      ;Get direction of line
      drdi = INTERPOLATE(DERIV(sep_ri), inds[j])
      dzdi = INTERPOLATE(DERIV(sep_zi), inds[j])
      
      ; First check if this is towards or away from the X-point
      dir = 1 ; direction to go away from x-point
      
      d = drdi*(cri - xpt_ri) + dzdi*(czi - xpt_zi) ; Dot-product
      IF d LT 0. THEN dir = -1
      
      ; Get the indices for a line radiating from the x-point
      si = [inds[j]]
      IF dir GT 0 THEN BEGIN
        in = CEIL(si[0])
        si = [si, in + indgen(N_ELEMENTS(sep_ri) - in)]
      ENDIF ELSE BEGIN
        in = FLOOR(si[0])
        si = [si, reverse(indgen(in+1))]
      ENDELSE
      sepri = INTERPOLATE(sep_ri, si)
      sepzi = INTERPOLATE(sep_zi, si)
      
      ; Then check if this is towards the O-point (core) or away (PF)
      
      d = drdi*dir * (opt_ri - xpt_ri) + dzdi*dir * (opt_zi - xpt_zi)
      
      IF d GT 0 THEN BEGIN
        ; Core
        
        ; find where either ri or zi has an extrema
        n = N_ELEMENTS(sepri)
        dr = DERIV(sepri)
        dr = dr[1:*] * dr[0:(n-2)]
        dz = DERIV(sepzi)
        dz = dz[1:*] * dz[0:(n-2)]
        in = MIN(WHERE((dr[1:*] LE 0.0) OR (dz[1:*] LE 0.0))) + 1
        
        sepri = sepri[0:in]
        sepzi = sepzi[0:in]
        
        IF KEYWORD_SET(debug) THEN OPLOT, sepri, sepzi, color=4
        
        IF ncore EQ 0 THEN BEGIN
          core1 = [[sepri], [sepzi]]
          ncore = 1
        ENDIF ELSE core2 = [[sepri], [sepzi]]
        
      ENDIF ELSE BEGIN
        ; PF
        
        sepri = INTERPOLATE(sep_ri, si)
        sepzi = INTERPOLATE(sep_zi, si)
        
        ; Find where it crosses a boundary
        
        cpos = line_crossings(sepri, sepzi, 0, $
                              boundary[0,*], boundary[1,*], 1, $
                              ncross=ncross, inds1=in)
        
        IF ncross GT 0 THEN BEGIN
          sepri = [sepri[0:FLOOR(in[0])], INTERPOLATE(sepri, in[0])]
          sepzi = [sepzi[0:FLOOR(in[0])], INTERPOLATE(sepzi, in[0])]
        ENDIF
        
        IF KEYWORD_SET(debug) THEN OPLOT, sepri, sepzi, color=3
        
        IF npf EQ 0 THEN BEGIN
          pf1 = [[sepri], [sepzi]]
          npf = 1
        ENDIF ELSE pf2 = [[sepri], [sepzi]]
        
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
    OPLOT, INTERPOLATE(R, pf1[*,0]), INTERPOLATE(Z, pf1[*,1]), color=3, thick=2
    OPLOT, INTERPOLATE(R, pf2[*,0]), INTERPOLATE(Z, pf2[*,1]), color=4, thick=2
    
    OPLOT, INTERPOLATE(R, core1[*,0]), INTERPOLATE(Z, core1[*,1]), color=3, thick=2
    OPLOT, INTERPOLATE(R, core2[*,0]), INTERPOLATE(Z, core2[*,1]), color=4, thick=2
  ENDELSE

  RETURN, {leg1:pf1, leg2:pf2, core1:core1, core2:core2, ri:xpt_ri, zi:xpt_zi}
END
