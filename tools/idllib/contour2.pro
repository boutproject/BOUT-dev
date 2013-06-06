; A replacement for contour for plotting BOUT results
; Uses a blue-white-red (blue -ve, red +ve) color scheme
; 
; Data can be either 2 or 3D
; If 3D, scales color scheme based on entire range
; 

PRO contour2, data, x, y, t=t, nlev=nlev, centre=centre, redblue=redblue, color=color, $
              revcolor=revcolor, _extra=_extra

  IF NOT KEYWORD_SET(t) THEN t = 0
  IF NOT KEYWORD_SET(nlev) THEN nlev=100

  np = N_PARAMS() ; number of parameters used
  dim = SIZE(data, /dim)
  
  IF np LT 3 THEN y = FINDGEN(dim[1])
  IF np LT 2 THEN x = FINDGEN(dim[0])

  ; get data range
  mind = MIN(data)
  maxd = MAX(data)
  
  IF KEYWORD_SET(centre) THEN BEGIN
    ; make zero white
    IF mind + maxd GT 0.0 THEN mind = -maxd ELSE maxd = -mind
  ENDIF
  
  lev=mind + (maxd-mind)*indgen(nLev)/(nLev-1)
  col=255*indgen(nLev)/(nLev-1)

  IF !D.NAME EQ 'X' THEN DEVICE,decomposed=0
  
  IF KEYWORD_SET(redblue) THEN BEGIN
    ; Define red-blue color table
    common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
    
    red = BYTARR(256)
    green = red
    blue = red
    
    ; need to keep color[0] = white, color[1] = black
    red[0] = 255
    green[0] = 255
    blue[0] = 255
    
    ; now create scale
    
    FOR i=2, 255 DO BEGIN
      green[i] = 256 - 2*ABS(i - 128.5)
      
      IF i GT 129 THEN blue[i] = 256 - 2*ABS(i - 128.5) ELSE blue[i] = 255
      IF i LE 129 THEN red[i] = 256 - 2*ABS(i - 128.5) ELSE red[i] = 255
      
    ENDFOR
    
    tvlct,red,green,blue
    
    col=2+253*indgen(nLev)/(nLev-1)
    
    color = 1
  ENDIF

  IF KEYWORD_SET(revcolor) THEN col = REVERSE(col)

  nd = SIZE(data, /n_dim)
  IF nd EQ 2 THEN BEGIN
    d = data
  ENDIF ELSE IF nd EQ 3 THEN BEGIN
    d = data[*,*,t]
  ENDIF ELSE BEGIN
    PRINT, "ERROR; incorrect number of dimensions"
  ENDELSE
  
  ; Plot data
  
  CONTOUR, d, x, y, $
    /fil, lev=lev, c_col=col, _extra=_extra, $
    POSITION=[0.15, 0.15, 0.95, 0.8], color=color

  ; Save the axes
  xaxis = !x
  yaxis = !y
  p = !p

  ; Location for color bar
  loc = [0.15, 0.90, 0.95, 0.95]
  
  bar = col # REPLICATE(1B, 10)
  
  xsize = (loc(2) - loc(0)) * !D.X_VSIZE
  ysize = (loc(3) - loc(1)) * !D.Y_VSIZE 
  xstart = loc(0) * !D.X_VSIZE
  ystart = loc(1) * !D.Y_VSIZE 
  
  ;bar = BYTSCL(bar, TOP=nLev-1)
  
  IF !D.NAME EQ 'PS' THEN $
    TV, bar, xstart, ystart, XSIZE=xsize, YSIZE=ysize ELSE $
    TV, CONGRID(bar, xsize, ysize), xstart, ystart
  
  ;PLOTS, [loc(0), loc(0), loc(2), loc(2), loc(0)], $
  ;  [loc(1), loc(3), loc(3), loc(1), loc(1)], /NORMAL, color=color
  
  PLOT, lev, [0,1], POSITION=loc, /nodata, /noerase, yticks=1, yminor=1, ystyle=4, xstyle=1, color=color
  
  ; Restore the contour plot axes
  !x = xaxis
  !y = yaxis
  !p = p
END
