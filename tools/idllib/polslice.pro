; code to produce poloidal cross-section of BOUT data
; handles single or double null

FUNCTION indgen2, a, b
  IF b GT a THEN BEGIN
    n = b - a + 1
    i = indgen(n)+a
  ENDIF ELSE BEGIN
    n = a-b+1
    i = REVERSE(indgen(n))+b
  ENDELSE
  RETURN, i
END

PRO polslice, data, gridfile, xstart=xstart, ystart=ystart, $
              nlevel=nlevel, level=level, color=color, $
              title=title, output=output
  
  ;barpos = [1.7,0.5] ; position of bottom-left corner

  s = SIZE(gridfile, /TYPE)
  IF s EQ 8 THEN BEGIN
      ; gridfile contains the grid data
      g = gridfile
  ENDIF ELSE IF s EQ 7 THEN BEGIN
      ; read the data
      g = file_import(gridfile)
  ENDIF ELSE BEGIN
      PRINT, "ERROR: gridfile must be either a structure or a filename"
      RETURN
  ENDELSE

  IF NOT KEYWORD_SET(xstart) THEN xstart = 0
  IF NOT KEYWORD_SET(ystart) THEN ystart = 0

  if not keyword_set(NLEVEL) then NLEVEL=100
  if not keyword_set(LEVEL) then level=min(data)+(max(data)-min(data))*findgen(NLEVEL+1)/NLEVEL
  if not keyword_set(COLOR) then color=20+2.4*100*findgen(NLEVEL+1)/NLEVEL

  sg = SIZE(g.rxy, /dimensions)
  sd = SIZE(data, /dimensions)
  
  rxy = g.rxy
  zxy = g.zxy

  ;;;;;;;;;;;;;; SET LAYOUT ;;;;;;;;;;;;;;;;;

  barwidth = (MAX(rxy) - MIN(rxy))/8.0
  barheight = (MAX(zxy) - MIN(zxy))*3./4.0
  barpos = [MAX(rxy)+barwidth,MIN(zxy) + (MAX(zxy) - MIN(zxy))/8.0]

  xrange=[MAX([0.0,MIN(rxy)-3.0*barwidth]), MAX(rxy)+3.0*barwidth]
  yrange=[MIN(zxy), MAX(zxy)]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  tn = TAG_NAMES(g)
  tn = STRUPCASE(tn)
  
  w = WHERE(tn EQ "IXSEPS1", count)
  IF count GT 0 THEN ixseps = g.ixseps1 - xstart ELSE ixseps = g.nx - ystart
  IF ixseps+1 GE g.nx THEN ixseps = g.nx-2
  w = WHERE(tn EQ "JYSEPS1_1", count)
  IF count GT 0 THEN jyseps1_1 = g.jyseps1_1 - ystart ELSE jyseps1_1 = -1 - ystart
  w = WHERE(tn EQ "JYSEPS1_2", count)
  IF count GT 0 THEN jyseps1_2 = g.jyseps1_2 - ystart ELSE jyseps1_2 = (g.ny/2) - ystart
  w = WHERE(tn EQ "JYSEPS2_1", count)
  IF count GT 0 THEN jyseps2_1 = g.jyseps2_1 - ystart ELSE jyseps2_1 = (g.ny/2) - ystart
  w = WHERE(tn EQ "JYSEPS2_2", count)
  IF count GT 0 THEN jyseps2_2 = g.jyseps2_2 - ystart ELSE jyseps2_2 = (g.ny-1) - ystart
  w = WHERE(tn EQ "IXLB2", count)
  IF count GT 0 THEN ixlb2 = g.ixlb2 - ystart ELSE ixlb2 = g.ny/2 - ystart
 

  IF sd[0] NE sg[0] THEN BEGIN
    IF sd[0] GT sg[0] THEN BEGIN
      PRINT, "Data x dimension larger than grid!"
      RETURN
    ENDIF ELSE BEGIN
      rxy = rxy[xstart:(sg[0]-1-xstart), *]
      zxy = zxy[xstart:(sg[0]-1-xstart), *]
    ENDELSE
  ENDIF

  IF sd[1] NE sg[1] THEN BEGIN
    IF sd[1] GT sg[1] THEN BEGIN
      PRINT, "Data y dimension larger than grid!"
      RETURN
    ENDIF ELSE BEGIN
      rxy = rxy[*,ystart:(sg[1]-1-ystart)]
      zxy = zxy[*,ystart:(sg[1]-1-ystart)]
    ENDELSE
  ENDIF

  nx = sd[0]
  ny = sd[1]

  xs = xrange[1] - xrange[0]
  ys = yrange[1] - yrange[0]

  IF KEYWORD_SET(output) THEN BEGIN
    !P.FONT=0
    SET_PLOT, 'PS'
    DEVICE, filename=output, bits=16, XSIZE=xs*5, YSIZE=ys*5, /color
  ENDIF

  safe_colors, /first

  ; plot core
  
  yinds = [indgen2(jyseps1_1+1, jyseps2_1), indgen2(jyseps1_2+1, jyseps2_2), jyseps1_1+1]
  
  contour, data[0:(ixseps+1), yinds], rxy[0:(ixseps+1), yinds], zxy[0:(ixseps+1), yinds], $
    nlevel=nlevel, /fill, /iso, c_col=color, level=level, $
    xrange=xrange, yrange=yrange, title=title

  ;cursor, x, y, /down
  
  IF jyseps1_1 GE 0 THEN BEGIN

      IF jyseps1_2 EQ jyseps2_1 THEN BEGIN
                                ; Single null
          
          PRINT, "SINGLE NULL"
          
                                ; SOL
          contour, data[ixseps:*, *], rxy[ixseps:*,*], zxy[ixseps:*,*], /overplot, $
            /fill, nlevel=nlevel, level=level, c_col=color
          
                                ;cursor, x, y, /down
          
                                ; PF
          yinds = [indgen2(0, jyseps1_1), indgen2(jyseps2_2+1, ny-1)]
          contour, data[0:ixseps, yinds], rxy[0:ixseps, yinds], zxy[0:ixseps, yinds], $
            /overplot, /fill, nlevel=nlevel, level=level, c_col=color
          
                                ; cursor, x, y, /down
          
                                ; fill in the hole
          
          
          xind = [25,25,25,25]
          yind = [jyseps1_1, jyseps1_1+1, jyseps2_2, jyseps2_2+1]
          
          d = [[data[25,jyseps1_1], data[25, jyseps1_1+1]], [data[25,jyseps2_2+1], data[25,jyseps2_2]]]
          
          r = [[rxy[25,jyseps1_1], rxy[25, jyseps1_1+1]], [rxy[25,jyseps2_2+1], rxy[25,jyseps2_2]]]
          z = [[zxy[25,jyseps1_1], zxy[25, jyseps1_1+1]], [zxy[25,jyseps2_2+1], zxy[25,jyseps2_2]]]
          
          contour, d, r, z, $
      /overplot, /fill, nlevel=nlevel, level=level, c_col=color
          
                                ; plot edges
          
          yinds = [indgen2(jyseps1_1+1, jyseps1_2), indgen2(jyseps2_1+1, jyseps2_2), jyseps1_1+1]
          oplot, rxy[0,yinds], zxy[0, yinds], thick=2, color=1
          
          oplot, rxy[nx-1,*], zxy[nx-1, *], thick=2, color=1
    
          oplot, rxy[*,0], zxy[*,0], thick=2, color=1
          oplot, rxy[*,ny-1], zxy[*,ny-1], thick=2, color=1

          yinds = [indgen2(0, jyseps1_1), indgen2(jyseps2_2+1, ny-1)]
          
          oplot, rxy[0,yinds], zxy[0,yinds], thick=2, color=1
          
                                ; plot separatrix
          
          yinds = [indgen2(0,jyseps1_1), indgen2(jyseps2_2,jyseps1_1+1), indgen2(jyseps2_2+1, ny-1)]
          
          oplot, rxy[25,yinds], zxy[25,yinds], thick=3, color=1, lines=1
          
      ENDIF ELSE BEGIN
                                ; Double NULL
          PRINT, "DOUBLE NULL"
          
                                ;cursor, x, y, /down
          
                                ; Outer SOL
          contour, data[ixseps:*, (ixlb2-2):*], rxy[ixseps:*,(ixlb2-2):*], $
            zxy[ixseps:*,(ixlb2-2):*], /overplot, $
            /fill, nlevel=nlevel, level=level, c_col=color
          
                                ;cursor, x, y, /down
          
                                ; Inner SOL
          
          contour, data[ixseps:*, 0:(ixlb2-3)], rxy[ixseps:*,0:(ixlb2-3)], $
            zxy[ixseps:*,0:(ixlb2-3)], /overplot, $
            /fill, nlevel=nlevel, level=level, c_col=color
          
                                ; cursor, x,y, /down
          
                                ; lower PF
          yinds = [indgen2(0, jyseps1_1), indgen2(jyseps2_2+1, ny-1)]
          contour, data[0:ixseps, yinds], rxy[0:ixseps, yinds], zxy[0:ixseps, yinds], $
            /overplot, /fill, nlevel=nlevel, level=level, c_col=color
          
                                ; upper PF
          
                                ;cursor, x, y, /down
          
          yinds = [indgen2(ixlb2-3, jyseps2_1+1), indgen2(jyseps1_2, ixlb2-2)]
          contour, data[0:ixseps, yinds], rxy[0:ixseps, yinds], zxy[0:ixseps, yinds], $
            /overplot, /fill, nlevel=nlevel, level=level, c_col=color

                                ;cursor, x, y, /down
          
                                ; fill gaps
          
          d = [[data[25,jyseps1_1], data[25, jyseps1_1+1]], [data[25,jyseps2_2+1], data[25,jyseps2_2]]]
          r = [[rxy[25,jyseps1_1], rxy[25, jyseps1_1+1]], [rxy[25,jyseps2_2+1], rxy[25,jyseps2_2]]]
          z = [[zxy[25,jyseps1_1], zxy[25, jyseps1_1+1]], [zxy[25,jyseps2_2+1], zxy[25,jyseps2_2]]]
          
          contour, d, r, z, $
            /overplot, /fill, nlevel=nlevel, level=level, c_col=color
          
          d = [[data[25,jyseps2_1], data[25, jyseps2_1+1]], [data[25,jyseps1_2+1], data[25,jyseps1_2]]]
          r = [[rxy[25,jyseps2_1], rxy[25, jyseps2_1+1]], [rxy[25,jyseps1_2+1], rxy[25,jyseps1_2]]]
          z = [[zxy[25,jyseps2_1], zxy[25, jyseps2_1+1]], [zxy[25,jyseps1_2+1], zxy[25,jyseps1_2]]]
          
          contour, d, r, z, $
            /overplot, /fill, nlevel=nlevel, level=level, c_col=color
          
      ENDELSE
  ENDIF

  ; add color bar

  bx = [[0, barwidth], [0, barwidth]] + barpos[0]
  by = [[0,0],[barheight, barheight]] + barpos[1]

  bdata = [[level[0], level[0]], [level[nlevel-1], level[nlevel-1]]]
  
  contour, bdata, bx, by, $
    /overplot, /fill, nlevel=nlevel, level=level, c_col=color

  ; print labels

  xyouts, barpos[0]+barwidth, barpos[1], STRING(level[0])
  xyouts, barpos[0]+barwidth, barpos[1]+barheight, STRING(level[nlevel-1])

  IF KEYWORD_SET(output) THEN BEGIN
    DEVICE, /close
    SET_PLOT, 'X'
  ENDIF
END

