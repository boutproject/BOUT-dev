; Calculate the growth-rate and frequency

PRO gamma_diamag, t, signal, gamma=gamma, freq=freq
  ; Find local maxima where d/dx = 0 and d2/dx2 < 0

  ddx = DERIV(t, signal)
  d2dx2 = DERIV(t, ddx)
  
  n = N_ELEMENTS(signal)

  tp = WHERE( (ddx[1:*]*ddx[0:(n-2)] LE 0.0) AND (d2dx2[1:*] LT 0.0), count)

  plot, t, signal
  oplot, t[tp], signal[tp], psym=7
  
  ; calculate growth-rate and frequency
  gpeak = DERIV(t[tp], ALOG(signal[tp]))
  fpeak = 0.5 / DERIV(t[tp])
 
  ; interpolate onto original grid
  gamma = INTERPOL(gpeak, t[tp], t)
  freq = INTERPOL(fpeak, t[tp], t)
  
END
