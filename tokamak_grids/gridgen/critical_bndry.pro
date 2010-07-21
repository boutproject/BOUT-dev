;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make sure all critical points are inside the boundary
;
; bndryi[2,n]   r and z indices of the boundary
FUNCTION critical_bndry, critical, bndryi
  ; First the o-points
  
  primary_opt = critical.primary_opt
  opt_ri = critical.opt_ri
  opt_zi = critical.opt_zi
  opt_f  = critical.opt_f
  
  IF critical.n_opoint GT 1 THEN BEGIN
    w = [critical.primary_opt]
    primary_opt = 0
    
    FOR i=0,critical.n_opoint-1 DO BEGIN
      IF i NE critical.primary_opt THEN BEGIN
        ; Test if outside boundary by drawing a line to the primary o-point
        cp = line_crossings([opt_ri[i], opt_ri[primary_opt]], $
                            [opt_zi[i], opt_zi[primary_opt]], $
                            0, $
                            bndryi[0,*], bndryi[1,*], 1, $
                            ncross=ncross)
        IF ncross EQ 0 THEN w = [w,i]
      ENDIF
    ENDFOR
    
    n_opoint = N_ELEMENTS(w)
    opt_ri = opt_ri[w]
    opt_zi = opt_zi[w]
    opt_f  = opt_f[w]
  ENDIF ELSE BEGIN
    n_opoint = 1
    opt_ri = critical.opt_ri
    opt_zi = critical.opt_zi
    opt_f  = critical.opt_f
  ENDELSE

  ; Check the x-points
  n_xpoint = 0
  FOR i=0, critical.n_xpoint-1 DO BEGIN
    ; Test if outside boundary by drawing a line to the primary o-point
    cp = line_crossings([critical.xpt_ri[i], opt_ri[primary_opt]], $
                        [critical.xpt_zi[i], opt_zi[primary_opt]], $
                        0, $
                        bndryi[0,*], bndryi[1,*], 1, $
                        ncross=ncross)
    IF ncross EQ 0 THEN BEGIN
      ; hasn't crossed -> inside boundary. Add index to w
      IF n_xpoint EQ 0 THEN w = [i] ELSE w = [w,i]
      n_xpoint = n_xpoint + 1
    ENDIF
  ENDFOR
  
  IF n_xpoint EQ 0 THEN BEGIN
    ; Keep the inner sep (used for normalisation)
    w = [critical.inner_sep]
    inner_sep = 0
  ENDIF ELSE BEGIN
    ; Need to update inner_sep
    w2 = WHERE(w EQ critical.inner_sep)
    inner_sep = w2[0]
  ENDELSE

  ; Select x-points
  xpt_ri = critical.xpt_ri[w]
  xpt_zi = critical.xpt_zi[w]
  xpt_f  = critical.xpt_f[w]

  result = {n_opoint:n_opoint, n_xpoint:n_xpoint, $
            primary_opt:primary_opt, $
            inner_sep:inner_sep, $
            opt_ri:opt_ri, opt_zi:opt_zi, opt_f:opt_f, $
            xpt_ri:xpt_ri, xpt_zi:xpt_zi, xpt_f:xpt_f}
  RETURN, result
END
