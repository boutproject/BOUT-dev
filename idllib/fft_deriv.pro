;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FFT_DERIV: Calculates the derivative of a variable on a         ;
; periodic domain.                                                ;
;                                                                 ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION fft_deriv, var
  on_error, 2
  
  n = N_ELEMENTS(var)

  F = FFT(var)

  imag = COMPLEX(0.0, 1.0)

  F[0] = 0.0

  IF (n MOD 2) EQ 0 THEN BEGIN
      ; even number
      FOR i=1, n/2-1 DO BEGIN
          a = imag*2.0*!PI*FLOAT(i)/FLOAT(n)
          F[i] = F[i] * a         ; positive frequencies
          F[n-i] = - F[n-i] * a   ; negative frequencies
      ENDFOR
      F[n/2] = F[n/2] * (imag*!PI)
  ENDIF ELSE BEGIN
      ; odd number
      FOR i=1, (n-1)/2 DO BEGIN
          a = imag*2.0*!PI*FLOAT(i)/FLOAT(n)
          F[i] = F[i] * a
          F[n-i] = - F[n-i] * a 
      ENDFOR
  ENDELSE

  result = FFT(F, 1)
  
  RETURN, result
END
