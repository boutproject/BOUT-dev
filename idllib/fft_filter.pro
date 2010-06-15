;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FFT_FILTER: Fourier filter a variable on a periodic domain      ;
; Supply number of Fourier components to keep                     ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION fft_filter, var, nf
  on_error, 2
  
  n = N_ELEMENTS(var)

  F = FFT(var)

  imag = COMPLEX(0.0, 1.0)
  
  IF (n MOD 2) EQ 0 THEN BEGIN
    ; even number
    FOR i=nf+1, n/2-1 DO BEGIN
      F[i] = 0.0
      F[n-i] = 0.0
    ENDFOR
    IF nf LT n/2 THEN F[n/2] = 0. ; The Nyquist frequency
  ENDIF ELSE BEGIN
    ; odd number
    FOR i=nf+1, (n-1)/2 DO BEGIN
      F[i] = 0.0
      F[n-i] = 0.0
    ENDFOR
  ENDELSE
  
  RETURN, FFT(F, 1)
END
