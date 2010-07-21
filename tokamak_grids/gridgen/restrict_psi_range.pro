;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Select only those X- and O-points within the given range of
; normalised psi
FUNCTION restrict_psi_range, critical, psi_outer
  
  IF critical.n_xpoint EQ 0 THEN RETURN, critical

  xpt_psi = (critical.xpt_f - critical.opt_f[critical.primary_opt]) $
    / (critical.xpt_f[critical.inner_sep] - critical.opt_f[critical.primary_opt])
  
  PRINT, "X-point locations: "
  FOR i=0, critical.n_xpoint-1 DO BEGIN
    PRINT, "  "+STR(i)+": "+STR(xpt_psi[i])
  ENDFOR

  ; Select the x-points within this range
  w = WHERE(xpt_psi LT psi_outer, count)
  PRINT, "Number of x-points in range: "+STR(count)
  PRINT, ""
  
  IF count EQ 0 THEN BEGIN
    ; Keep the inner sep
    w = [critical.inner_sep]
    inner_sep = 0
  ENDIF ELSE BEGIN
    ; Need to update inner_sep
    w2 = WHERE(w EQ critical.inner_sep)
    inner_sep = w2[0]
  ENDELSE
  
  result = {n_opoint:critical.n_opoint, n_xpoint:count, $
            primary_opt:critical.primary_opt, $
            inner_sep:inner_sep, $
            opt_ri:critical.opt_ri, opt_zi:critical.opt_zi, opt_f:critical.opt_f, $
            xpt_ri:critical.xpt_ri[w], xpt_zi:critical.xpt_zi[w], xpt_f:critical.xpt_f[w]}
  RETURN, result
END
