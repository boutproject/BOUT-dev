; Discrete Cosine Transform
; 
; Based on algorithm here:
;
; http://fourier.eng.hmc.edu/e161/lectures/dct/node2.html
;
; Ben Dudson, Feb 2010
;
;

; DCT for any number of points
; Uses an FFT of 2N points, so not as efficient
; as method below
;
; NOTE: Inverse not working!
FUNCTION DCTany, x, inverse=inverse
  j = DCOMPLEX(0.0D, 1.0D)
  
  n = N_ELEMENTS(x)
  
  data = x
  IF KEYWORD_SET(inverse) THEN BEGIN
    data[0] = data[0] / SQRT(2.D)
    PRINT, "SORRY, NOT WORKING YET"
    STOP
  ENDIF

  ; Double length
  x2 = [data, REVERSE(data)]
  
  ; Take FFT of 2N
  f = FFT(x2, /double)
  result = REAL_PART(f[0:(N-1)] * EXP(-j*DINDGEN(N)*!DPI / (2.D*N)) )
  
  IF NOT KEYWORD_SET(inverse) THEN result[0] = result[0] / SQRT(2.D)
  
  RETURN, result
END

; This method is more efficient, but only works
; for even N
FUNCTION DCT, x, inverse=inverse
  
  n = N_ELEMENTS(x)
  IF n MOD 2 NE 0 THEN BEGIN
    ; Odd number of points: Need to use more general method
    RETURN, DCTany(x, inverse=inverse)
  ENDIF

  j = DCOMPLEX(0.0D, 1.0D)

  IF NOT KEYWORD_SET(inverse) THEN BEGIN
    
    ; Re-order array
    y = DBLARR(n)
    FOR i=0, n/2 - 1 DO BEGIN
      y[i] = x[2*i]
      y[N-1-i] = x[2*i + 1]
    ENDFOR
    
    ; Calculate FFT
    yf = FFT(DOUBLE(y), /double)
    
    ; Multiply by phase
    result = REAL_PART(yf * EXP(-j*DINDGEN(n)*!DPI/(2.D*DOUBLE(n))) )
    
    result[0] = result[0] / SQRT(2.d)
    
    RETURN, result
  ENDIF ELSE BEGIN
    
    yf = DCOMPLEX(x) * EXP(j*DINDGEN(n)*!DPI/(2.D*DOUBLE(n)))
    yf[0] = yf[0] / SQRT(2.D)
    y = REAL_PART(FFT(yf, /inverse, /double))
    
    result = DBLARR(n)
    
    FOR i=0, n/2 - 1 DO BEGIN
      result[2*i] = y[i]
      result[2*i+1] = y[n-1-i]
    ENDFOR
    
    RETURN, 2.D*result
  ENDELSE
END
