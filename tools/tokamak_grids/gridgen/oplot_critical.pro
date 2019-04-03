; Overplots the output of analyse_equil
PRO oplot_critical, F, R, Z, a
  
  ; Plot X-points and separatrices
  FOR i=0, a.n_xpoint-1 DO BEGIN
    ; plot the separatrix contour
    CONTOUR, F, R, Z, levels=[a.xpt_f[i]], c_colors=2, /overplot
    oplot, [INTERPOLATE(R, a.xpt_ri[i], /DOUBLE)], [INTERPOLATE(Z, a.xpt_zi[i], /DOUBLE)], psym=7, color=2
  ENDFOR

  ; Plot O-points
  FOR i=0, a.n_opoint-1 DO BEGIN
    oplot, [INTERPOLATE(R, a.opt_ri[i], /DOUBLE)], [INTERPOLATE(Z, a.opt_zi[i], /DOUBLE)], psym=7, color=3
  ENDFOR
END
