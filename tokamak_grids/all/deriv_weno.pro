FUNCTION deriv_weno, x, f
  ; do the middle section

  n = N_ELEMENTS(f)

  df = FLTARR(n)

  eps = 1.0e-3

  FOR i=1, n-2 DO BEGIN
      ; derivatives using different stencils
      pr   = (f[i+1] - f[i]) / (x[i+1] - x[i])
      pl   = (f[i] - f[i-1]) / (x[i] - x[i-1])
      popt = (f[i+1] - f[i-1]) / (x[i+1] - x[i-1])
      pc   = 2.0*popt - 0.5*(pr + pl)

      ; smoothness indicators

      isl = (f[i] - f[i-1])^2
      isr = (f[i+1] - f[i])^2
      isc = (13./3.)*(f[i+1] - 2.*f[i] + f[i-1])^2 + 0.25*(f[i+1] - f[i-1])^2

      ; weights

      wl = 0.25 / (eps + isl)^2
      wr = 0.25 / (eps + isr)^2
      wc = 0.5  / (eps + isc)^2
      total = wl + wr + wc

      df[i] = (wl*pl + wr*pr + wc*pc)/total
      
  ENDFOR
  
  ; end points

  df[0] = (f[1] - f[0])/(x[1] - x[0])
  df[n-1] = (f[n-1] - f[n-2])/(x[n-1] - x[n-2])

  RETURN, df
END

