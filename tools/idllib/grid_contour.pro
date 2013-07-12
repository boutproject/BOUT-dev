; A contour plot using BOUT-06 style input grid (for tokamaks)
;
; Inputs
; ------
; 
;  g          A structure containing a grid input file
;  data2d     2D data [nx, ny] to display
;
;  nlev       Sets number of colour levels
;  scale      If set, displays a side-bar with the scale
;
;  Extra keywords passed through to the initial call to plot
; 

PRO grid_contour, g, data2d, nlev=nlev, scale=scale, maxd=maxd, mind=mind, _extra=_extra
  
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
  IF NOT KEYWORD_SET(mind) THEN mind = MIN(data2d) ELSE mind = FLOAT(mind)
  IF NOT KEYWORD_SET(maxd) THEN maxd = MAX(data2d) ELSE maxd = FLOAT(maxd)
  
  

  IF KEYWORD_SET(centre) THEN BEGIN
    ; make zero white
    IF mind + maxd GT 0.0 THEN mind = -maxd ELSE maxd = -mind
  ENDIF

  lev=mind + (maxd-mind)*indgen(nLev)/(nLev-1)
  col=2+253*indgen(nLev)/(nLev-1)

  ; plot to get ranges
  plot, g.rxy, g.zxy, /iso, /nodata, _extra=_extra
  
  IF ix2 EQ nx THEN ix2 = ix1

  IF ix2 EQ ix1 THEN BEGIN
    ; Single null or connected double 
    
    IF ix1 GT nx-1 THEN BEGIN
      ix1 = nx-1
      ix2 = ix1
    ENDIF

    ; core
    yinds = [range(j1_1+1, j1_2), range(j2_1+1, j2_2), j1_1+1]
    
    contour, data2d[0:ix1, yinds], $
      g.rxy[0:ix1, yinds], g.zxy[0:ix1, yinds], $
      /over, /fill, lev=lev, c_col=col
    
    IF ix1 LT nx-1 THEN BEGIN
      ; X-point, otherwise just core
      
      IF j1_2 EQ j2_1 THEN BEGIN
        ; No upper PF region. Only one SOL
        yinds = range(0, ny-1)
        
        contour, data2d[ix1:*, yinds], $
          g.rxy[ix1:*, yinds], g.zxy[ix1:*, yinds], $
          /over, /fill, lev=lev, c_col=col
        
      ENDIF ELSE BEGIN
        ; Inner SOL
        
        
      ENDELSE
      
      ; Lower PF region
      yinds = [range(0, j1_1), range(j2_2+1,ny-1)]
      contour, data2d[0:ix1, yinds], $
        g.rxy[0:ix1, yinds], g.zxy[0:ix1, yinds], $
        /over, /fill, lev=lev, c_col=col
      
      ; Lower x-point
      d = [[data2d[ix1,j1_1], data2d[ix1,j2_2+1]], [data2d[ix1,j1_1+1], data2d[ix1,j2_2]]]
      r = [[g.rxy[ix1,j1_1], g.rxy[ix1,j2_2+1]], [g.rxy[ix1,j1_1+1], g.rxy[ix1,j2_2]]]
      z = [[g.zxy[ix1,j1_1], g.zxy[ix1,j2_2+1]], [g.zxy[ix1,j1_1+1], g.zxy[ix1,j2_2]]]
      
      contour, d, r, z, $
        /over, /fill, lev=lev, c_col=col
    ENDIF
  ENDIF ELSE BEGIN
    ; Unbalanced double null
    
    ny_inner = g.ny_inner
    
    ; core
    yinds = [range(j1_1+1, j2_1), range(j1_2+1, j2_2), j1_1+1]
    
    contour, data2d[0:ix1, yinds], $
      g.rxy[0:ix1, yinds], g.zxy[0:ix1, yinds], $
      /over, /fill, lev=lev, c_col=col
    
    ; Lower PF region
    
    yinds = [range(0, j1_1), range(j2_2+1, ny-1)]
    
    contour, data2d[0:ix1, yinds], $
      g.rxy[0:ix1, yinds], g.zxy[0:ix1, yinds], $
      /over, /fill, lev=lev, c_col=col

    ; Upper PF region
    
    yinds = [range(ny_inner, j1_2), range(j2_1+1, ny_inner-1)]
    contour, data2d[0:ix2, yinds], $
      g.rxy[0:ix2, yinds], g.zxy[0:ix2, yinds], $
      /over, /fill, lev=lev, c_col=col

    ; intra-separatrix region
    IF ix1 LT ix2 THEN BEGIN
      ; Lower 
      
      yinds = [range(0, j2_1), range(j1_2+1, ny-1)]
      contour, data2d[ix1:ix2, yinds], $
        g.rxy[ix1:ix2, yinds], g.zxy[ix1:ix2, yinds], $
        /over, /fill, lev=lev, c_col=col
      
    ENDIF ELSE BEGIN
      ; Upper
      
      yinds = [range(ny_inner, j2_2), range(j1_1+1, ny_inner-1)]
      
      contour, data2d[ix2:ix1, yinds], $
        g.rxy[ix2:ix1, yinds], g.zxy[ix2:ix1, yinds], $
        /over, /fill, lev=lev, c_col=col
      
    ENDELSE
    
    ; Inner SOL
    ix = MAX([ix1, ix2])
    
    yinds = range(0, ny_inner-1)
    contour, data2d[ix:*, yinds], $
      g.rxy[ix:*, yinds], g.zxy[ix:*, yinds], $
      /over, /fill, lev=lev, c_col=col
    
    ; Outer SOL
    
    yinds = range(ny_inner, ny-1)
    contour, data2d[ix:*, yinds], $
      g.rxy[ix:*, yinds], g.zxy[ix:*, yinds], $
      /over, /fill, lev=lev, c_col=col
    
    ; Lower x-point
    d = [[data2d[ix1,j1_1], data2d[ix1,j2_2+1]], [data2d[ix1,j1_1+1], data2d[ix1,j2_2]]]
    r = [[g.rxy[ix1,j1_1], g.rxy[ix1,j2_2+1]], [g.rxy[ix1,j1_1+1], g.rxy[ix1,j2_2]]]
    z = [[g.zxy[ix1,j1_1], g.zxy[ix1,j2_2+1]], [g.zxy[ix1,j1_1+1], g.zxy[ix1,j2_2]]]

    contour, d, r, z, $
      /over, /fill, lev=lev, c_col=col

    ; Upper x-point
    d = [[data2d[ix2,j2_1], data2d[ix2,j1_2+1]], [data2d[ix2,j2_1+1], data2d[ix2,j1_2]]]
    r = [[g.rxy[ix2,j2_1], g.rxy[ix2,j1_2+1]], [g.rxy[ix2,j2_1+1], g.rxy[ix2,j1_2]]]
    z = [[g.zxy[ix2,j2_1], g.zxy[ix2,j1_2+1]], [g.zxy[ix2,j2_1+1], g.zxy[ix2,j1_2]]]

    contour, d, r, z, $
      /over, /fill, lev=lev, c_col=col

  ENDELSE

  IF KEYWORD_SET(scale) THEN BEGIN
    ; Save the axes
    xaxis = !x
    yaxis = !y
    p = !p
    
    ; Location for color bar
    loc = [0.90, 0.15, 0.95, 0.95]
    
    bar = REPLICATE(1B, 10) # col
    
    xsize = (loc(2) - loc(0)) * !D.X_VSIZE
    ysize = (loc(3) - loc(1)) * !D.Y_VSIZE 
    xstart = loc(0) * !D.X_VSIZE
    ystart = loc(1) * !D.Y_VSIZE 
    
    ;bar = BYTSCL(bar, TOP=nLev-1)
    
    
    IF !D.NAME EQ 'PS' THEN $
      TV, bar, xstart, ystart, XSIZE=xsize, YSIZE=ysize ELSE $
      TV, CONGRID(bar, xsize, ysize), xstart, ystart
    
    PLOT, lev, POSITION=loc, /nodata, /noerase, xticks=1, xminor=1, xstyle=4, ystyle=1, color=color
    
    ; Restore the contour plot axes
    !x = xaxis
    !y = yaxis
    !p = p
  ENDIF
  
END
