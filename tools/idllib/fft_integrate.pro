FUNCTION fft_integrate, y, loop=loop
  on_error,2  ; If an error occurs, return to caller
  
  n = N_ELEMENTS(y)

  F = FFT(y)
  imag = complex(0.0, 1.0)

  result = FINDGEN(n)*F[0]
  loop = FLOAT(n)*F[0] ;; return the loop integral

  F[0] = 0.0

  IF (n MOD 2) EQ 0 THEN BEGIN
      ; even number of points

      FOR i=1l, n/2-1 DO BEGIN
          a = imag*2.0*!PI*FLOAT(i)/FLOAT(n)
          F[i] = F[i] / a         ; positive frequencies
          F[n-i] = - F[n-i] / a   ; negative frequencies
      ENDFOR

      F[n/2] = F[n/2] / (imag*!PI)
  ENDIF ELSE BEGIN
      ; odd number of points

      FOR i=1l, (n-1)/2 DO BEGIN
          a = imag*2.0*!PI*FLOAT(i)/FLOAT(n)
          F[i] = F[i] / a
          F[n-i] = - F[n-i] / a 
      ENDFOR
  ENDELSE

  result = result + FFT(F, 1)

  result = result - result[0] ; just to make sure

  RETURN, result
END

PRO test_integrate
  n = 10

  dx = 2.0*!PI/FLOAT(n)
  x = dx * FINDGEN(n)

  y = 1 + COS(x) - 0.5*sin(2*x)
  iy = x + SIN(x) + 0.25*cos(2*x)

  result = fft_integrate(y)*dx

  plot, x, iy
  oplot, x, result, psym=1
  
END
