; Interpolates a quantity so that the derivatives are smooth
; Intended for interpolating pressure profiles

FUNCTION interp_smooth, f, x, x2
  ; First differentiate original, then interpolate
  d = interpol(DERIV(x, f), x, x2,/spline)
  ; Integrate up to get profile
  p = int_func(x2, d, /simple)

  ; Now p may not quite match f.
  
  ; match inner and outer values
  in = (interpol(f, x, x2[0], /spline))[0]

  m = MIN(ABS(x2 - x[N_ELEMENTS(x)-1]), oind) ; last index to use
  out = (interpol(f, x, x2[oind], /spline))[0]

  pin = in - p[0]
  pout = out - p[oind]

  p[0:oind] = p[0:oind] + pin + (pout - pin)*(x2 - x2[0])/(x2[N_ELEMENTS(x2)-1] - x2[0])
  IF oind LT (N_ELEMENTS(x2)-1) THEN p[(oind+1):*] = p[(oind+1):*] + pout
  
  RETURN, p
END
