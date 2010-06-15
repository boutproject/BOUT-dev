; Given the FFT of a function (from IDL's FFT routine), works out
; the value and derivative at any given x location

FUNCTION fft_interp, F, x, ddx=ddx
  n = N_ELEMENTS(F)

  result = F[0]
  ddx = COMPLEX(0.0,0.0)
  
  imag = complex(0.0, 1.0)

  IF (n MOD 2) EQ 0 THEN BEGIN
      ; even number of points

      ; positive and negative frequencies
      FOR i=1, n/2 - 1 DO BEGIN
          freq = FLOAT(i) / FLOAT(n)
          a = imag*2.0*!PI*freq
          e = EXP(a*x)

          result = result + 2.0*F[i] * e
          ddx = ddx + F[i]*a*e - F[n-i]*a*EXP(-1.0*a*x)
      ENDFOR

      result = result + F[n/2]*EXP(imag*!PI*x)
      ddx = ddx + F[n/2]*imag*!PI*EXP(imag*!PI*x)
      
  ENDIF ELSE BEGIN
      ; odd number of points

      FOR i=1, (n-1)/2 DO BEGIN
          freq = FLOAT(i) / FLOAT(n)
          a = imag*2.0*!PI*freq
          e = EXP(a*x)

          result = result + F[i]*e + F[n-i]*EXP(-1.0*a*x)
          ddx = ddx + F[i]*a*e - F[n-i]*a*EXP(-1.0*a*x)
      ENDFOR
  ENDELSE

  ddx = REAL_PART(ddx)

  RETURN, REAL_PART(result)
END


PRO test_interp
   n = 11
   seed = 123

   ; random test

   x = RANDOMN(seed, n)

   f = FFT(x)

   x2 = FLTARR(n)
   FOR I=0, n-1 DO BEGIN
       x2[i] = fft_interp(f, FLOAT(i))
       PRINT, i, x[i], x2[i]
   ENDFOR

   ; derivatives test
   x = 2.0*!PI*FINDGEN(n)/FLOAT(n)
   dx = 2.0*!PI / FLOAT(n)

   y = SIN(x)
   dydx = COS(x)

   f = FFT(y)
   
   y2 = FLTARR(n)
   dydx2 = FLTARR(n)
   FOR i=0, n-1 DO BEGIN
       y2[i] = fft_interp(f, FLOAT(i), d=d)
       dydx2[i] = d / dx

       PRINT, i, x[i], y[i], y2[i], dydx[i], dydx2[i]
   ENDFOR
   

   STOP
END
