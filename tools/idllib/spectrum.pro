; Calculate Fourier spectrum

FUNCTION periodogram, data, freq=freq
  n = N_ELEMENTS(data)
  ; Multiply by Hann window and FFT
  f = FFT(data*0.5*(1.0 - cos(2.0*!PI*DINDGEN(n)/DOUBLE(n))))
  
  freq = (1+FINDGEN(FIX(n/2))) / FLOAT(n)
  
  RETURN, ABS(f[1:FIX(n/2)])^2
END

FUNCTION spectrum, data, nwindows=nwindows, freq=freq
  n = N_ELEMENTS(data)
  
  IF NOT KEYWORD_SET(nwindows) THEN nwindows = SQRT(n)
  len = 2*FIX(n / (nwindows+1)) ; overlapping by 1/2
  
  f = periodogram(data[0:(len-1)], freq=freq)
  FOR i=0, nwindows-1 DO BEGIN
    s = i * len / 2
    f = f + periodogram(data[s:(s+len-1)])
  ENDFOR
  
  f = f / DOUBLE(nwindows)

  RETURN, f
END
