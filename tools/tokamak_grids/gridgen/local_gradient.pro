;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find the gradient at a given point by fitting
; 
; Uses DCT to interpolate and get gradients
; 
; dctF   - the DCT of F
; ri, zi - R and Z indices
;

FUNCTION local_gradient, dctF, ri, zi, status=status, f=f
  s = SIZE(dctF, /DIMENSION)
  nr = s[0]
  nz = s[1]
  
  IF (ri LT 0) OR (ri GT nr-1) OR (zi LT 0) OR (zi GT nz-1) THEN BEGIN
    status = 1
    RETURN, 0
  ENDIF
  
  IF KEYWORD_SET(f) THEN BEGIN
    ; Use a local approximation. Faster, but less accurate
    
    x = ROUND(ri)
    y = ROUND(zi)
    
    COMMON lgsvd, r, z, A, W, U, V
    
    IF NOT KEYWORD_SET(W) THEN BEGIN
      ; Calculate the SVD of the matrix once then re-use
      
      PRINT, "Calculating SVD for local gradient"

      ;r = [-1, -1, -1, 0, 0, 0, 1, 1, 1]
      ;z = [-1, 0,   1, -1, 0, 1, -1, 0, 1]
      
      r = [-1,  0, 0, 0, 1]
      z = [ 0, -1, 0, 1, 0]
      
      n = N_ELEMENTS(r)

      A = TRANSPOSE([[FLTARR(n)+1.], $
                     [r], $
                     [z], $
                     [r*r], $
                     [z*z]])
      
      SVDC, A,W,U,V
    ENDIF
    vals = f[ ((x+r) > 0) < (nr-1), ((y+z) > 0) < (nz-1) ]
    
    d = SVSOL(U,W,V,vals)
    
    status = 0
    
    dx = ri - x
    dy = zi - y

    RETURN, {f: d[0] + d[1]*dx + d[2]*dy + d[3]*dx*dx + d[4]*dy*dy , $
             dfdr: d[1] + d[3]*dx, $
             dfdz: d[2] + d[4]*dy}
  ENDIF ELSE BEGIN
    ; Use DCT method
    res = EvalCosPfast(dctF, x0=ri, y0=zi)
  ENDELSE

  status = 0

  RETURN, {f:res[0], dfdr:res[1], dfdz:res[2]}
END

