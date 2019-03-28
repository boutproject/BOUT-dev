; Fast 2D Discrete Cosine Transform
;
; http://fourier.eng.hmc.edu/e161/lectures/dct/node2.html
;
; Ben Dudson, University of York, Feb 2010
; 
; NOTE: SOMETHING NOT QUITE RIGHT HERE

FUNCTION DCT2D, sig, inverse=inverse
  s = SIZE(sig, /dim)
  IF N_ELEMENTS(s) NE 2 THEN BEGIN
    PRINT, "ERROR: input to DCT2Dfast must be 2D"
    RETURN, 0
  ENDIF
  nx = s[0]
  ny = s[1]

  result = DBLARR(nx, ny)

  FOR i=0, ny-1 DO BEGIN
    result[*,i] = DCT(sig[*,i], inverse=inverse)
  ENDFOR
  
  FOR i=0, nx-1 DO BEGIN
    result[i,*] = DCT(result[i,*], inverse=inverse)
  ENDFOR
  
  IF NOT KEYWORD_SET(inverse) THEN BEGIN
    result = result * 2.D * SQRT(nx*ny)
  ENDIF ELSE BEGIN
    result = result / (2.D* SQRT(nx*ny))
  ENDELSE

  RETURN, result
END
