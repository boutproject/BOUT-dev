; A replacement for contour for plotting BOUT results
; Uses a blue-white-red (blue -ve, red +ve) color scheme
; 
; Data can be either 2 or 3D
; If 3D, scales color scheme based on entire range
; 

PRO contour2, data, x, y, t=t, nlev=nlev, centre=centre, $
              xrange=xrange, yrange=yrange, $
              xtitle=xtitle, ytitle=ytitle, $
              charsize=charsize, $
              xstyle=xstyle, ystyle=ystyle, _extra=_extra

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
  col=2+253*indgen(nLev)/(nLev-1)
  
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

  PLOT, x, y, /nodata, color=1, $
    xtitle=xtitle, ytitle=ytitle, $
    xrange=xrange, yrange=yrange, $
    charsize=charsize, $
    xstyle=xstyle, ystyle=ystyle

  nd = SIZE(data, /n_dim)
  IF nd EQ 2 THEN BEGIN
    CONTOUR, data, x, y, /over, /fil, lev=lev, c_col=col
  ENDIF ELSE IF nd EQ 3 THEN BEGIN
    CONTOUR, data[*,*,t], x, y, /over, /fil, lev=lev, c_col=col, _extra=_extra
  ENDIF ELSE BEGIN
    PRINT, "ERROR; incorrect number of dimensions"
  ENDELSE

END
