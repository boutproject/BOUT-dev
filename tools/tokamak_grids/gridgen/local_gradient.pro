;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find the gradient at a given point by fitting
; 
; This code is used a LOT so quite time-critical
; 
; interp_data  - Structure containing data for interpolation
; ri, zi - R and Z indices
;

PRO local_gradient, interp_data, ri, zi, status=status, $
                    f=f, dfdr=dfdr, dfdz=dfdz
  nr = interp_data.nx
  nz = interp_data.ny
  
  IF (ri LT 0) OR (ri GT nr-1) OR (zi LT 0) OR (zi GT nz-1) THEN BEGIN
    status = 1
    RETURN
  ENDIF
  
  CASE interp_data.method OF
    0: BEGIN ; DCT method
      res = EvalCosPfast(interp_data.dct, x0=ri, y0=zi)
      f    = res[0]
      dfdr = res[1]
      dfdz = res[2]
    END
    1: BEGIN ; SVD local derivatives
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
      vals = interp_data.f[ ((x+r) > 0) < (nr-1), ((y+z) > 0) < (nz-1) ]
      
      d = SVSOL(U,W,V,vals)
      
      status = 0
      
      dx = ri - x
      dy = zi - y

      f    = d[0] + d[1]*dx + d[2]*dy + d[3]*dx*dx + d[4]*dy*dy
      dfdr = d[1] + d[3]*dx
      dfdz = d[2] + d[4]*dy
    END
    ELSE: BEGIN
      PRINT, "ERROR: unknown method in local_gradient"
      STOP
    END
  ENDCASE
  status = 0
END

