; Plots an R-Z equilibrium structure

PRO plot_rz_equil, data, _extra=_extra
  nlev = 100
  minf = MIN(data.psi)
  maxf = MAX(data.psi)
  levels = findgen(nlev)*(maxf-minf)/DOUBLE(nlev-1) + minf
  
  safe_colors, /first
  CONTOUR, data.psi, data.r, data.z, levels=levels, /iso, color=1, $
    /xsty,/ysty, xtit="Major radius [m]", ytit="Height [m]", _extra=_extra
  
  IF in_struct(data, "critical") THEN BEGIN
    critical = data.critical
    IF in_struct(data, "nlim") THEN BEGIN
      IF data.nlim GT 2 THEN BEGIN
        
        OPLOT, [REFORM(data.rlim), data.rlim[0]], [REFORM(data.zlim), data.zlim[0]], $
          thick=2,color=2
        
        ; Check that the critical points are inside the boundary
        bndryi = DBLARR(2, data.nlim)
        bndryi[0,*] = INTERPOL(FINDGEN(data.nr), data.R, data.rlim)
        bndryi[1,*] = INTERPOL(FINDGEN(data.nz), data.Z, data.zlim)
        
        critical = critical_bndry(critical, bndryi)
      ENDIF
      ; Overplot the separatrices, O-points
      oplot_critical, data.psi, data.r, data.z, critical
    ENDIF
  ENDIF
END
