;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; X-point refinement 


FUNCTION ref_xpt_newt, X
  ; Return gradients
  COMMON rfn_xpt_com, interp_data
  
  local_gradient, interp_data, X[0], X[1], dfdr=dfdr, dfdz=dfdz
  
  RETURN, [dfdr, dfdz]
END

PRO refine_xpoint, interp_data, ri0, zi0, ri, zi
  ; Find where Grad(F) = 0
  
  COMMON rfn_xpt_com, fdata
  fdata = interp_data
  
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
