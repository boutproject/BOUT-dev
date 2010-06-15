;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find the gradient at a given point by fitting
; 
; Uses DCT to interpolate and get gradients
; 
; dctF   - the DCT of F
; ri, zi - R and Z indices
;

FUNCTION local_gradient, dctF, ri, zi, status=status
  s = SIZE(dctF, /DIMENSION)
  nr = s[0]
  nz = s[1]
  
  IF (ri LT 0) OR (ri GT nr-1) OR (zi LT 0) OR (zi GT nz-1) THEN BEGIN
    status = 1
    RETURN, 0
  ENDIF
  
  res = EvalCosPfast(dctF, x0=ri, y0=zi)
  
  status = 0

  RETURN, {f:res[0], dfdr:res[1], dfdz:res[2]}
END

