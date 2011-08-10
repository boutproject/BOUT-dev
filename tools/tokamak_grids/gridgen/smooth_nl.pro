; nonlinear smoothing

FUNCTION smooth_nl, input, mesh
  
  s = SIZE(input, /dim)
  nx = s[0]
  ny = s[1]
  
  output = input

  tmp = output

  REPEAT BEGIN
    
    mark = FLTARR(nx, ny)
    
    status = gen_surface(mesh=mesh) ; Start generator
    REPEAT BEGIN
      ; Get the next domain
      yi = gen_surface(period=period, last=last, x=x)
      
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
          
          mark[x,y] = ABS(dxm - dxp); + ABS(dym - dyp)
          
        ENDFOR
      ENDIF
    ENDREP UNTIL last
    
    mark = mark / MAX(mark)
    
    status = gen_surface(mesh=mesh) ; Start generator
    REPEAT BEGIN
      ; Get the next domain
      yi = gen_surface(period=period, last=last, x=x)
      
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
          m = 0.2*(mark[x,y] + $
                   mark[x-1,y] + mark[x+1,y] + $
                   mark[x,ym] + mark[x, yp])
          
          tmp[x,y] = (1.0-m)*output[x,y] + $
            m*0.25*(output[x-1,y] + output[x+1,y] + output[x,ym] + output[x,yp])
        ENDFOR
      ENDIF
    ENDREP UNTIL last
        
    change = MAX(ABS(tmp - output))
    output = tmp
    ;PRINT, change
  ENDREP UNTIL change LT 1e-2
  
  RETURN, output
END
