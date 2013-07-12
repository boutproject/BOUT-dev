; Plots an R-Z equilibrium structure

PRO plot_rz_equil, data
  nlev = 100
  minf = MIN(data.psi)
  maxf = MAX(data.psi)
  levels = findgen(nlev)*(maxf-minf)/FLOAT(nlev-1) + minf
  
  safe_colors, /first
  CONTOUR, data.psi, data.r, data.z, levels=levels, /iso, color=1
  IF data.nlim GT 2 THEN OPLOT, [data.rlim, data.rlim[0]], $
                                [data.zlim, data.zlim[0]], $
                                color = 2, thick=2
  critical = data.critical
  IF data.nlim GT 2 THEN BEGIN
    OPLOT, [REFORM(data.rlim), data.rlim[0]], [REFORM(data.zlim), data.zlim[0]], $
           thick=2,color=2
    
    ; Check that the critical points are inside the boundary
    bndryi = FLTARR(2, data.nlim)
    bndryi[0,*] = INTERPOL(FINDGEN(data.nr), data.R, data.rlim)
    bndryi[1,*] = INTERPOL(FINDGEN(data.nz), data.Z, data.zlim)
    
    critical = critical_bndry(critical, bndryi)
  ENDIF
  ; Overplot the separatrices, O-points
  oplot_critical, data.psi, data.r, data.z, critical
END
