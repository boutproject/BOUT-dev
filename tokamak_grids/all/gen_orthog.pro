; works out a transformation to orthogonal coordinates by changing
; theta. Returns the old theta coordinates of the new points

FUNCTION gen_orthog, r, z, inside_theta, tol=tol
  s = SIZE(r, /dimensions)
  nx = s[0]
  ny = s[1]

  IF NOT KEYWORD_SET(tol) THEN tol = 1.0e-3

  nlines = N_ELEMENTS(inside_theta) ; the number of resulting coordinates

  ; calculate FFT of r and z poloidally
  fr = COMPLEXARR(nx, ny)
  fz = COMPLEXARR(nx, ny)
  FOR i=0, nx-1 DO BEGIN
      fr[i,*] = FFT(r[i,*])
      fz[i,*] = FFT(z[i,*])
  ENDFOR

  ; result is theta locations for all new points
  result = FLTARR(nx, nlines)
  result[0,*] = inside_theta*FLOAT(ny)/(2.0*!PI) ; turn into indices

  FOR j=0, nlines-1 DO BEGIN ; loop over all output poloidal points
      
      ; FIRST PASS: Work from the inside out
      
      FOR i=0, nx-2 DO BEGIN
          ; get perpendicular vector from starting point
          r0 = fft_interp(fr[i,*], result[i,j], d=drdt0)
          z0 = fft_interp(fz[i,*], result[i,j], d=dzdt0)
          result[i+1,j] = result[i,j]

          dt = 0.0
          iter = 0
          REPEAT BEGIN
              olddt = dt
              ; find where this vector intersects next flux surface
              
              ; vector at current point
              r1 = fft_interp(fr[i+1,*], result[i+1,j], d=drdt1)
              z1 = fft_interp(fz[i+1,*], result[i+1,j], d=dzdt1)

              dt = (drdt0*(r1 - r0) + dzdt0*(z1 - z0)) / (dzdt0*dzdt1 + drdt0*drdt1)
              result[i+1,j] = result[i+1,j] - dt
              ;PRINT, i, dt
              ;cursor, x, y, /down

              iter = iter + 1
              IF (iter GT ny) AND (ABS(dt) GT ABS(olddt)) THEN BEGIN
                  PRINT, "Error: Gen_orthog diverging"
                  STOP
              ENDIF
              
          ENDREP UNTIL ABS(dt) LT tol
          ;PRINT, i
          ;STOP
      ENDFOR


      WRITEU, -1, 13, "Lines done:"+STRTRIM(STRING(j+1),2)+" of "+STRTRIM(STRING(nlines),2)
      ;PRINT, "Lines done: ", j+1
      ; could do higher-order corrections as a relaxation process here

  ENDFOR
  PRINT, " "

  ;STOP
  result = ((result MOD ny) + ny) MOD ny ; make sure indices between 0 and ny
  result = result*2.0*!PI / FLOAT(ny)    ; turn back into an angle

  RETURN, result
END
