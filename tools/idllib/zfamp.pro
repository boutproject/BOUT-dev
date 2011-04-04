;; Given a 4D variable [x,y,z,t], returns the Fourier amplitudes in
;; [x,y,f,t]

FUNCTION zfamp, var4d
  ON_ERROR, 2
  
  s = SIZE(var4d, /dim)

  IF N_ELEMENTS(s) NE 4 THEN BEGIN
    PRINT, "Error: zfamp takes 4D variables"
    RETURN, 0
  ENDIF

  nx = s[0]
  ny = s[1]
  nz = s[2]
  nt = s[3]


  IF (nz MOD 2) EQ 1 THEN nz = nz - 1

  IF NOT is_pow2(nz) THEN BEGIN
    PRINT, "Error: Expected a z length of power of 2 or power-of-2 + 1"
    RETURN, 0
  ENDIF
  
  result = FLTARR(nx, ny, (nz/2)+1, nt)
  
  FOR x=0, nx-1 DO BEGIN
    FOR y=0, ny-1 DO BEGIN
      FOR t=0, nt-1 DO BEGIN
        f = FFT(var4d[x,y,0:(nz-1), t])
        result[x,y,0:(nz/2), t] = 2.*ABS(f[0:(nz/2)])
      ENDFOR
    ENDFOR
  ENDFOR
  
  RETURN, result
END

