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
      IF NOT in_struct(interp_data, "dct") THEN BEGIN
        PRINT, "Calculating DCT for method 0"
        
        dct = DCT2D(interp_data.F)
        interp_data = CREATE_STRUCT(interp_data, "dct", dct)
      ENDIF
      
      res = EvalCosPfast(interp_data.dct, x0=ri, y0=zi)
      f    = res[0]
      dfdr = res[1]
      dfdz = res[2]
    END
    1: BEGIN ; SVD local derivatives
      x = ROUND(ri)
      y = ROUND(zi)
      
      IF NOT in_struct(interp_data, "method1") THEN BEGIN
      
        ; Calculate the SVD of the matrix once then re-use
       
        PRINT, "Calculating SVD for local gradient (method 1)"

        ;r = [-1, -1, -1, 0, 0, 0, 1, 1, 1]
        ;z = [-1, 0,   1, -1, 0, 1, -1, 0, 1]
        
        r = [-1,  0, 0, 0, 1]
        z = [ 0, -1, 0, 1, 0]
        
        n = N_ELEMENTS(r)

        A = TRANSPOSE([[DBLARR(n)+1.D], $
                       [r], $
                       [z], $
                       [r*r], $
                       [z*z]])
        
        SVDC, A,W,U,V
        
        d = {r:r, z:z, A:A, W:W, U:U, V:V}
        
        interp_data = CREATE_STRUCT(interp_data, "method1", d)
      ENDIF ELSE d = interp_data.method1
      
      vals = interp_data.f[ ((x+d.r) > 0) < (nr-1), ((y+d.z) > 0) < (nz-1) ]
      
      d = SVSOL(d.U,d.W,d.V,vals)
      
      status = 0
      
      dx = ri - x
      dy = zi - y

      f    = d[0] + d[1]*dx + d[2]*dy + d[3]*dx*dx + d[4]*dy*dy
      dfdr = d[1] + d[3]*dx
      dfdz = d[2] + d[4]*dy
    END
    2: BEGIN
      IF NOT in_struct(interp_data, "method2") THEN BEGIN
        ; Calculate derivatives
        
        PRINT, "Calculating derivatives for local gradient (method 2)"

        ddr = DBLARR(nr, nz)
        ddz = ddr
        FOR i=0, nz-1 DO ddr[*,i] = DERIV(interp_data.f[*,i])
        FOR i=0, nr-1 DO ddz[i,*] = DERIV(interp_data.f[i,*])
        
        d = {ddr:ddr, ddz:ddz}
        
        interp_data = CREATE_STRUCT(interp_data, "method2", d)
      ENDIF ELSE d = interp_data.method2
      
      IF ARG_PRESENT(f)    THEN f    = INTERPOLATE(interp_data.f, ri, zi, cubic=-0.5D, /DOUBLE)
      IF ARG_PRESENT(dfdr) THEN dfdr = INTERPOLATE(d.ddr, ri, zi, cubic=-0.5D, /DOUBLE)
      IF ARG_PRESENT(dfdz) THEN dfdz = INTERPOLATE(d.ddz, ri, zi, cubic=-0.5D, /DOUBLE)
    END
    ELSE: BEGIN
      PRINT, "ERROR: unknown method in local_gradient"
      STOP
    END
  ENDCASE
  status = 0
END

