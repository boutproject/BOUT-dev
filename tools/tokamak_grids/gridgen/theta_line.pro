FUNCTION theta_line, dctF, ri0, zi0, di0, nstep, boundary=boundary, dir=dir
  COMMON td_com, fdata, lastgoodpos
  
  ri = [ri0]
  zi = [zi0]
  fdata = dctF
  
  pos = [ri0,zi0]
  
  di = di0
  IF KEYWORD_SET(dir) THEN BEGIN
    ; Set direction to go in
    
    dt = theta_differential(0., pos)
    IF TOTAL(dt*dir) LT 0. THEN di = -di
  ENDIF
  
  FOR i=1, nstep DO BEGIN
    CATCH, theError
    IF theError EQ 0 THEN BEGIN
      pos = LSODE(pos, 0, di, 'theta_differential')
      ri = [ri, pos[0]]
      zi = [zi, pos[1]]
    ENDIF ELSE BEGIN
      ; Error occurred. Should have been caused by hitting a boundary
      ; and/or going off the end of the domain
      CATCH, /cancel
      ri = [ri, lastgoodpos[0]]
      zi = [zi, lastgoodpos[1]]
      RETURN, [[ri], [zi]]
    ENDELSE

    IF KEYWORD_SET(boundary) THEN BEGIN
      ; Check if crossed a boundary
      n = N_ELEMENTS(ri)
      cpos = line_crossings(ri[(n-2):*], zi[(n-2):*], 0, $
                            boundary[0,*], boundary[1,*], 1, ncross=ncross)
      IF ncross GT 0 THEN BEGIN
        ; Change the last point to the intersection location
        
        ri[n-1] = cpos[0,0]
        zi[n-1] = cpos[1,0]
        
        RETURN, [[ri], [zi]]
      ENDIF
    ENDIF
  ENDFOR
  RETURN, [[ri], [zi]]
END
