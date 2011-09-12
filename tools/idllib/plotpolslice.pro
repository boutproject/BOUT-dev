; PLOTS 2D poloidal cross-section
; specify toroidal angle zangle in degrees


; interpolates a 1D periodic function
FUNCTION zinterp, v, zind
  v = REFORM(v)
  nz = N_ELEMENTS(v)
  z0 = ROUND(zind)

  p = zind - FLOAT(z0)          ; between -0.5 and 0.5
  
  IF p LT 0.0 THEN BEGIN
      z0 = z0 - 1
      p = p + 1.0
  ENDIF
  
  z0 = ((z0 MOD (nz-1)) + (nz-1)) MOD (nz-1)
  
  ; for now 3-point interpolation
          
  zp = (z0 + 1) MOD (nz - 1)
  zm = (z0 - 1 + (nz-1)) MOD (nz - 1)
  
  result = 0.5*p*(p-1.0)*v[zm] $
    + (1.0 - p*p)*v[z0] $
    + 0.5*p*(p+1.0)*v[zp]

  RETURN, result
END

PRO plotpolslice, var3d_in, uedge, period=period, zangle=zangle, $
                  nlev=nlev, yr=yr, _extra=_extra, profile=profile, $
                  output=output, lines=lines, linecol=linecol, filter=filter

  ON_ERROR, 2

  IF KEYWORD_SET(filter) THEN BEGIN
    var3d = zfilter(var3d_in, filter)
  ENDIF ELSE var3d = var3d_in

  IF max(var3d) - min(var3d) LT 1.e-30 THEN BEGIN
    PRINT, "Data has no range top plot!"
    RETURN
  ENDIF

  s = SIZE(var3d, /dimensions)
  IF N_ELEMENTS(s) NE 3 THEN BEGIN
      PRINT, "Error: Variable must be 3 dimensional"
      RETURN
  ENDIF
  IF NOT KEYWORD_SET(period) THEN period = 1 ; default = full torus
  IF NOT KEYWORD_SET(zangle) THEN zangle = 0.0
  IF NOT KEYWORD_SET(nlev) THEN nlev = 100
  IF NOT KEYWORD_SET(linecol) THEN linecol = 1 ; color for extra lines

  IF KEYWORD_SET(output) THEN BEGIN
      PRINT, "Outputting to "+output
      SET_PLOT, 'PS'
      DEVICE, file=output, /color
  ENDIF

  period = FIX(period) ; make sure it's an integer
  IF period LT 1 THEN period = 1
  zangle = FLOAT(zangle) * !PI / 180. ; convert from degrees to radians
  
  nx = s[0]
  ny = s[1]
  nz = s[2]

  dz = 2.0*!PI / FLOAT(period*(nz-1))

  ; GET THE TOROIDAL SHIFT
  tn = TAG_NAMES(uedge)
  tn = STRUPCASE(tn)
  w = WHERE(tn EQ "QINTY", count)
  IF count GT 0 THEN BEGIN
      PRINT, "Using qinty as toroidal shift angle"
      zShift = uedge.qinty
  ENDIF ELSE BEGIN
      w = WHERE(tn EQ "ZSHIFT", count)
      IF count GT 0 THEN BEGIN
          PRINT, "Using zShift as toroidal shift angle"
          zShift = uedge.zShift
      ENDIF ELSE BEGIN
          PRINT, "ERROR: Can't find qinty or zShift variable"
          RETURN
      ENDELSE
  ENDELSE

  ; calculate number of poloidal points

  ny2 = ny  ; start with just the original points
  nskip = FLTARR(ny-1)
  FOR i=0, ny-2 DO BEGIN
      ; get maximum number of extra points in-between i and i+1
      ip = (i + 1) MOD ny
      
      nskip[i] = 0
      FOR x=0, nx-1 DO BEGIN
          ns = MAX(ABS(zShift[x,ip] - zShift[x,i]))/dz - 1
          IF ns GT nskip[i] THEN nskip[i] = ns
      ENDFOR
  ENDFOR
 
  nskip = LONG(ROUND(nskip))
  ny2 = LONG(ny2 + TOTAL(nskip))

  PRINT, "Number of poloidal points in output:", ny2

  var2d = FLTARR(nx, ny2)
  rxy = FLTARR(nx, ny2)
  zxy = FLTARR(nx, ny2)

  ypos = 0l
  FOR y=0, ny-2 DO BEGIN
      ; put in the original points
      FOR x=0, nx-1 DO BEGIN
          zind = (zangle - zShift[x,y])/dz
          var2d[x,ypos] = zinterp(var3d[x,y,*], zind)
          IF KEYWORD_SET(profile) THEN var2d[x,ypos] = var2d[x,ypos] + profile[x,y]
          rxy[x,ypos] = uedge.rxy[x,y]
          zxy[x,ypos] = uedge.zxy[x,y]
      ENDFOR
      ypos = ypos + 1
      
      PRINT, y, ypos

      ; and the extra points

      FOR x=0, nx-1 DO BEGIN
          zi0 = (zangle - zShift[x,y])/dz
          zip = (zangle - zShift[x,y+1])/dz

          dzi = (zip - zi0) / (nskip[y] + 1)
          
          FOR i=0l, nskip[y]-1 DO BEGIN
              zi = zi0 + FLOAT(i+1)*dzi ; zindex 
              w = FLOAT(i+1)/FLOAT(nskip[y]+1) ; weighting
              
              var2d[x,ypos+i] = w*zinterp(var3d[x,y+1,*], zi) + (1.0-w)*zinterp(var3d[x,y,*], zi)
              IF KEYWORD_SET(profile) THEN var2d[x,ypos+i] = var2d[x,ypos+i] + w*profile[x,y+1] + (1.0-w)*profile[x,y]
              rxy[x,ypos+i] = w*uedge.rxy[x,y+1] + (1.0-w)*uedge.rxy[x,y]
              zxy[x,ypos+i] = w*uedge.zxy[x,y+1] + (1.0-w)*uedge.zxy[x,y]
          ENDFOR
      ENDFOR
      
      ypos = ypos + nskip[y]
  ENDFOR

  ; FINAL POINT
  
  FOR x=0, nx-1 DO BEGIN
      zind = (zangle - zShift[x,ny-1])/dz
      var2d[x,ypos] = zinterp(var3d[x,ny-1,*], zind)
      IF KEYWORD_SET(profile) THEN var2d[x,ypos] = var2d[x,ypos] + profile[x,ny-1]
      rxy[x,ypos] = uedge.rxy[x,ny-1]
      zxy[x,ypos] = uedge.zxy[x,ny-1]
  ENDFOR

  ; Get data range
  IF KEYWORD_SET(yr) THEN BEGIN
      mind = yr[0]
      maxd = yr[1]
  ENDIF ELSE BEGIN
      mind = MIN(var2d)
      maxd = MAX(var2d)

      ; balance the range
      IF maxd GT -1.0*mind THEN mind = -1.*maxd ELSE maxd = -1.*mind
  ENDELSE

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

  ;safe_colors, /first
  
  plot, rxy, zxy, psym=3,/iso, _extra=_extra, /nodata, color=1

  ;STOP 
  CONTOUR, var2d, rxy, zxy, /over, /fil, lev=lev, c_col=col
  
  y = ny2-1
  yp = 0
  d = [[reform(var2d[*,y])], [reform(var2d[*,yp])]]
  r = [[reform(rxy[*,y])], [reform(rxy[*,yp])]]
  z = [[reform(zxy[*,y])], [reform(zxy[*,yp])]]
; create a solid/"filled" plot of the simulation domain
  ;CONTOUR, d, r, z, /over, /fil, lev=lev, c_col=col
  
  ; draw the inner boundary of the simulation domain
  ;oplot, [reform(rxy[0,*]), rxy[0,0]], [reform(zxy[0,*]),zxy[0,0]], color=1
  
  ; draw the outer boundary of the simulation domain
  oplot, [reform(rxy[nx-1,*]),rxy[nx-1,0]],[reform(zxy[nx-1,*]),zxy[nx-1,0]], color=1

  ;FOR y=0, ny2-1 DO BEGIN
  ;    IF y EQ ny2-1 THEN yp = 0 ELSE yp = y + 1
  ;    d = [[reform(var2d[*,y])], [reform(var2d[*,yp])]]
  ;    r = [[reform(rxy[*,y])], [reform(rxy[*,yp])]]
  ;    z = [[reform(zxy[*,y])], [reform(zxy[*,yp])]]
  ;    CONTOUR, d, r, z, /over, /irr, /fil, lev=lev, c_col=co
  ;ENDFOR
  
  ;STOP

  IF KEYWORD_SET(lines) THEN BEGIN
      FOR i=0, N_ELEMENTS(lines) - 1 DO BEGIN
          OPLOT, rxy[lines[i], *], zxy[lines[i], *], lines=2, color=linecol
      ENDFOR
  ENDIF

  IF KEYWORD_SET(output) THEN BEGIN
      DEVICE, /close
      SET_PLOT, 'X'
  ENDIF
END
