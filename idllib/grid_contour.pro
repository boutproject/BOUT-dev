; A contour plot using BOUT-06 style input grid (for tokamaks)

FUNCTION indgen2, min, max
  n = max - min + 1
  RETURN, indgen(n) + min
END

PRO grid_contour, g, data2d, nlev=nlev, _extra=_extra
  
  IF NOT KEYWORD_SET(nlev) THEN nlev=50

  nx = g.nx
  ny = g.ny

  j1_1 = g.jyseps1_1
  j1_2 = g.jyseps1_2
  
  j2_1 = g.jyseps2_1
  j2_2 = g.jyseps2_2
  
  ix1 = g.ixseps1
  ix2 = g.ixseps2
  
  ; get data range
  mind = MIN(data2d)
  maxd = MAX(data2d)
  
  IF KEYWORD_SET(centre) THEN BEGIN
    ; make zero white
    IF mind + maxd GT 0.0 THEN mind = -maxd ELSE maxd = -mind
  ENDIF

  lev=mind + (maxd-mind)*indgen(nLev)/(nLev-1)
  col=2+253*indgen(nLev)/(nLev-1)

  ; plot to get ranges
  plot, g.rxy, g.zxy, /iso, /nodata, _extra=_extra
  
  IF ix2 EQ ix1 THEN BEGIN
    ; Single null or connected double 
    
    ; core
    yinds = [indgen2(j1_1+1, j1_2), indgen2(j2_1+1, j2_2), j1_1+1]
    
    contour, data2d[0:ix1, yinds], $
      g.rxy[0:ix1, yinds], g.zxy[0:ix1, yinds], $
      /over, /fill, lev=lev, c_col=col
    
    IF j1_2 EQ j2_1 THEN BEGIN
      ; No upper PF region. Only one SOL
      yinds = indgen2(0, ny-1)
      
      contour, data2d[ix1:*, yinds], $
        g.rxy[ix1:*, yinds], g.zxy[ix1:*, yinds], $
        /over, /fill, lev=lev, c_col=col
      
    ENDIF ELSE BEGIN
      ; Inner SOL
      
      
    ENDELSE
    
    ; Lower PF region
    yinds = [indgen2(0, j1_1), indgen2(j2_2+1,ny-1)]
    contour, data2d[0:ix1, yinds], $
      g.rxy[0:ix1, yinds], g.zxy[0:ix1, yinds], $
      /over, /fill, lev=lev, c_col=col
    
    ; Lower x-point
    d = [[data2d[ix1,j1_1], data2d[ix1,j2_2+1]], [data2d[ix1,j1_1+1], data2d[ix1,j2_2]]]
    r = [[g.rxy[ix1,j1_1], g.rxy[ix1,j2_2+1]], [g.rxy[ix1,j1_1+1], g.rxy[ix1,j2_2]]]
    z = [[g.zxy[ix1,j1_1], g.zxy[ix1,j2_2+1]], [g.zxy[ix1,j1_1+1], g.zxy[ix1,j2_2]]]

    contour, d, r, z, $
      /over, /fill, lev=lev, c_col=col
  ENDIF ELSE PRINT, "Sorry, can't handle unbalanced nulls yet"
  
END
