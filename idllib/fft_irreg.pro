; Calculate the Fourier Transform of an irregular set of points
;
; If the input is complex then computes both positive and negative
; frequencies, otherwise just the positive ones.
;
; Output is in the standard FFT layout:
; 
; DC, 1 ... (nf-1), nf(nyquist), -(nf-1) ... -1

FUNCTION fft_irreg, theta, data, nf=nf, sign=sign
  nd = N_ELEMENTS(data)
  IF NOT KEYWORD_SET(nf) THEN nf = FIX(nd/2)
  IF NOT KEYWORD_SET(sign) THEN sign = -1.0
  
  ; Check the type
  t = SIZE(data, /TYPE)
  cmplx = 0
  IF (t EQ 6) OR (t EQ 9) THEN BEGIN
    ; Complex or DComplex type
    ; Compute positive and negative frequencies
    cmplx = 1
    
    f = COMPLEXARR(2*nf)
  ENDIF ELSE BEGIN
    f = COMPLEXARR(nf+1)
  ENDELSE
  
  dth = ([theta[1:*], theta[0]] - theta + 2.*!PI) MOD 2*!PI
  
  FOR k=0, nf DO BEGIN
    
    ct = COS(FLOAT(k)*theta)
    st = SIN(FLOAT(k)*theta)
    
    f [k]= COMPLEX( TOTAL(ct * REAL_PART(data) - sign * st * IMAGINARY(data)), $
                    TOTAL(ct * IMAGINARY(data) + sign * st * REAL_PART(data)) )
  ENDFOR
  
  IF cmplx THEN BEGIN
    ; Calculate the negative frequencies
    FOR k=-(nf-1), -1 DO BEGIN
    
      ct = COS(FLOAT(k)*theta)
      st = SIN(FLOAT(k)*theta)
      
      f [2*nf+k]= COMPLEX( TOTAL(ct * REAL_PART(data) - sign * st * IMAGINARY(data)), $
                      TOTAL(ct * IMAGINARY(data) + sign * st * REAL_PART(data)) )
      
    ENDFOR
  ENDIF
  
  f = 2. * f / FLOAT(nd)
  
  RETURN, f
END
