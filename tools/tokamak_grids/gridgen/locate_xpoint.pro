;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Locate an x-point 
; 
; F        F(R,Z) psi array
; ri, zi   1D array of indices into F describing a path
; xpt_ri, xpt_zi  Location index of x-point
; xpt_f    F at the x-point
;

FUNCTION loc_xpt_newt, mini
  COMMON loc_xpt_newt, dctF, fri, fzi, xpt_ri, xpt_zi, xpt_f
  
  ; Follow gradient towards the separatrix
  follow_gradient, dctF, R, Z, fft_interp(fri, mini), fft_interp(fzi, mini), $
      xpt_f, eri, ezi
  
  ; Work out the distance to the x-point, including
  ; a sign to determine which direction to go in
  ; Cross product of the unit gradient and vector to the x-point
 
  ;PRINT, mini, fft_interp(fri, mini), fft_interp(fzi, mini)
  
  a = local_gradient(dctF, eri, ezi)

  gmag = SQRT(a.dfdr^2 + a.dfdz^2)
  gr = a.dfdr / gmag
  gz = a.dfdz / gmag
 
  ;PRINT, "  -> ", gmag, gr, gz, gr*(ezi - xpt_zi) - gz*(eri - xpt_ri)

  RETURN, gr*(ezi - xpt_zi) - gz*(eri - xpt_ri)
END

FUNCTION locate_xpoint, dctF, ri, zi, xpt_ri, xpt_zi, xpt_f, pos=pos
  COMMON loc_xpt_newt, fdata, fri, fzi, xr, xz, xf
  
  fdata = dctF
  fri = FFT(ri) ; for interpolating
  fzi = FFT(zi)
  xr = xpt_ri
  xz = xpt_zi
  xf = xpt_f
  
  ; Get starting index
  mind = MIN((ri - xpt_ri)^2 + (zi - xpt_zi)^2, mini)
  
  PRINT, "Starting index: ", mini

  CATCH, theError
  IF theError EQ 0 THEN BEGIN
    ; Use Newton method to find x-point
    mini = (NEWTON(mini, 'loc_xpt_newt'))[0]
    CATCH, /cancel
  ENDIF ELSE BEGIN
    CATCH, /cancel
    ; Just keep original mini
  ENDELSE
  

  ; interpolate to this position
  pos = [fft_interp(fri,mini), fft_interp(fzi,mini)]
  
  RETURN, mini MOD N_ELEMENTS(ri)
END
