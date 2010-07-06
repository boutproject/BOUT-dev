; Monotone interpolation using Hermite splines
; 
; x, y - input points
; u - points where result is needed

FUNCTION h00, t
  RETURN, (1. + 2.*t)*(1. - t)^2
END

FUNCTION h10, t
  RETURN, t*(1.-t)^2
END

FUNCTION h01, t
  RETURN, t^2*(3. - 2.*t)
END

FUNCTION h11, t
  RETURN, t^2*(t - 1.)
END

FUNCTION spline_mono, x, y, u, yp0=yp0, ypn_1=ypn_1
  n = N_ELEMENTS(x)

  IF n LT 3 THEN BEGIN
    ; Just linear interpolate
    RETURN, interpol(y,x,u)
  ENDIF

  ; Calculate delta
  D = (y[1:*] - y[0:(n-1)]) / (x[1:*] - x[0:(n-1)])
  
  IF N_ELEMENTS(yp0) EQ 0 THEN yp0 = D[0]
  IF N_ELEMENTS(ypn_1) EQ 0 THEN ypn_1 = D[n-2]
  
  m = FLTARR(n)
  m[1:(n-2)] = 0.5*(D[0:(n-3)] + D[1:(n-2)])
  m[0] = yp0
  
  FOR i=0, n-2 DO BEGIN
    IF ABS(D[i]) LT 1.e-6 THEN BEGIN
      m[i] = 0.
      m[i+1 < (n-1)] = 0.
    ENDIF ELSE BEGIN
      a = m[i] / D[i]
      b = m[i+1] / D[i]
      
      c = SQRT(a^2 + b^2)
      IF c GT 3. THEN BEGIN
        t = 3. / c
        m[i] = t * a * D[i]
        m[i+1] = t* b * D[i]
      ENDIF
    ENDELSE
  ENDFOR

  nout = N_ELEMENTS(u)
  
  result = FLTARR(nout)
  
  FOR i=0, nout-1 DO BEGIN
    xup = MIN(where(x GT u[i],count))
    IF count LE 0 THEN BEGIN
      result[i] = y[n-1]
      CONTINUE
    ENDIF
    IF xup EQ 0 THEN BEGIN
      result[i] = y[0]
      CONTINUE
    ENDIF
    xlow = xup - 1
    
    h = FLOAT(x[xup] - x[xlow])
    t = (FLOAT(u[i]) - FLOAT(x[xlow])) / h
    
    result[i] = y[xlow] * h00(t) + $
      h*m[xlow]*h10(t) + $
      y[xup]*h01(t) + $
      h*m[xup]*h11(t)
  ENDFOR
  
  RETURN, result
END
