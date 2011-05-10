;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; X-point refinement using DCT


FUNCTION ref_xpt_newt, X
  ; Return gradients
  COMMON rfn_xpt_com, dctf
  
  a = local_gradient(dctf, X[0], X[1], status=status)

  RETURN, [a.dfdr, a.dfdz]
END

PRO refine_xpoint, dctF, ri0, zi0, ri, zi
  ; Find where Grad(F) = 0
  
  COMMON rfn_xpt_com, fdata
  fdata = dctf
  
  CATCH, theError
  IF theError NE 0 THEN BEGIN
    CATCH, /cancel
    PRINT, "** Error occurred whilst refining x-point location"
    ri = ri0
    zi = zi0
    RETURN
  ENDIF
  
  rz0 = [ri0, zi0]
  rz = NEWTON(rz0, 'ref_xpt_newt')
  
  ri = rz[0]
  zi = rz[1]
  CATCH, /cancel
END
