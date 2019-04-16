; Plot a mesh as produced by create_grid.pro

PRO plot_mesh, mesh, overplot=overplot, _extra=_extra
  
  over = 0
  IF KEYWORD_SET(overplot) THEN over = 1
  
  ; Plot flux surfaces
  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface_hypnotoad(last=last, xi=xi, period=period)
    IF period THEN yi = [yi, yi[0]]
    IF over EQ 0 THEN BEGIN
      PLOT, mesh.Rxy[xi, yi], mesh.Zxy[xi, yi], $
        xr=[MIN(mesh.Rxy), MAX(mesh.Rxy)], yr=[MIN(mesh.Zxy), MAX(mesh.Zxy)], $
        /iso, thick=0.5D, _extra=_extra
      over = 1
    ENDIF ELSE BEGIN
      OPLOT, mesh.Rxy[xi, yi], mesh.Zxy[xi, yi], thick=0.5D
    ENDELSE  
  ENDREP UNTIL last
  
  ; Plot radial lines
  
  FOR i=0,TOTAL(mesh.npol)-1 DO BEGIN
    OPLOT, mesh.Rxy[*,i], mesh.Zxy[*,i], thick=1.5D
  ENDFOR
END
