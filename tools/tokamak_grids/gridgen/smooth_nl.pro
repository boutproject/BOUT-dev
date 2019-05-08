; nonlinear smoothing

FUNCTION smooth_nl, input, mesh, iter=iter
  
  IF NOT KEYWORD_SET(iter) THEN iter=50

  s = SIZE(input, /dim)
  nx = s[0]
  ny = s[1]
  
  output = input

  tmp = output

  markx = DBLARR(nx, ny)
  marky = markx
  
  mxn = markx
  myn = mxn
  
  it = 0
  REPEAT BEGIN  
    status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
    REPEAT BEGIN
      ; Get the next domain
      yi = gen_surface_hypnotoad(period=period, last=last, x=x)
      
      IF x GT 0 AND x LT nx-1 THEN BEGIN
        n = N_ELEMENTS(yi)
        FOR j=0, n-1 DO BEGIN
          IF period THEN BEGIN
            jm = (j-1 + n) MOD n
            jp = (j+1) MOD n
          ENDIF ELSE BEGIN
            jm = (j-1) > 0
            jp = (j+1) < (n-1)
          ENDELSE
          y = yi[j]
          ym = yi[jm]
          yp = yi[jp]
          
          dxm = output[x,y] - output[x-1,y]
          dxp = output[x+1,y] - output[x,y]
          
          dym = output[x,y] - output[x,ym]
          dyp = output[x,yp] - output[x,y]
          
          mxn[x,y] = 0.5D*(ABS(dxm) + ABS(dxp))
          myn[x,y] = 0.5D*(ABS(dym) + ABS(dyp))
          
        ENDFOR
      ENDIF
    ENDREP UNTIL last
    
    markx = (0.5D*mxn / MEAN(mxn)) < 1.0D
    marky = (0.5D*myn / MEAN(myn)) < 1.0D
    
    status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
    REPEAT BEGIN
      ; Get the next domain
      yi = gen_surface_hypnotoad(period=period, last=last, x=x)
      
      IF x GT 0 AND x LT nx-1 THEN BEGIN
        n = N_ELEMENTS(yi)
        FOR j=0, n-1 DO BEGIN
          IF period THEN BEGIN
            jm = (j-1 + n) MOD n
            jp = (j+1) MOD n
          ENDIF ELSE BEGIN
            jm = (j-1) > 0
            jp = (j+1) < (n-1)
          ENDELSE
          y = yi[j]
          ym = yi[jm]
          yp = yi[jp]
          
          ; Smooth the smoothing mask
          mx = 0.1D*(markx[x,y] + $
                    markx[x-1,y] + markx[x+1,y] + $
                    markx[x,ym] + markx[x, yp])
          
          my = 0.1D*(marky[x,y] + $
                     marky[x-1,y] + marky[x+1,y] + $
                     marky[x,ym] + marky[x, yp])
          
          tmp[x,y] = (1.0D - mx-my)*output[x,y] $
            + mx*0.5D*(output[x-1,y] + output[x+1,y])  $
            + my*0.5D*(output[x,ym] + output[x,yp])
        ENDFOR
      ENDIF
    ENDREP UNTIL last
    
    FOR y=0,ny-1 DO BEGIN
      tmp[0,y] = tmp[1,y]
      tmp[nx-1,y] = tmp[nx-2,y]
    ENDFOR
    
    change = MAX(ABS(tmp - output))
    output = tmp
    WRITEU, -1, 13, "Smoothing iteration"+STR(it)+" max change: "+STR(change)
    
    it = it + 1
  ENDREP UNTIL (change LT 1e-3) OR (it GE iter)
  
  PRINT, ""

  RETURN, output
END
