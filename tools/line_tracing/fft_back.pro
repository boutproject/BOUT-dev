;; Given a 3D variable [x,y,k] it returns the series
;; A(ix,jy,kz)=a_0(x,y)+SUM_OVER_K[ a_k(x,y)*cos(kz)+b_k(x,y)*sin(kz)]
;; a_0(ix,jy) is from var3d[ix,jy,0]
;; a_k(ix,jy) is from var3d[ix,jy,k]
;; b_k(ix,jy) is from var3d[ix,jy,NK-k] 
;; the amplitudes for the Nyquist frequency are set to zero
FUNCTION fft_back, var3d
  
  s = SIZE(var3d, /dim)

  IF N_ELEMENTS(s) NE 3 THEN BEGIN
    PRINT, "Error: fft_3d takes 3D variables"
    RETURN, 0
  ENDIF

  nx = s[0]
  ny = s[1]
  nz = s[2]

  IF (nz MOD 2) EQ 1 THEN nz = nz - 1

  IF NOT is_pow2(nz) THEN BEGIN
    PRINT, "Error: Expected a z length of power of 2 or power-of-2 + 1"
    RETURN, 0
  ENDIF
  
  result = FLTARR(nx, ny, nz)
  FOR ix=0, nx-1 DO BEGIN
    FOR jy=0, ny-1 DO BEGIN
      FOR kz=0, nz-1 DO BEGIN
        zval=DOUBLE(kz)/DOUBLE(nz)
        result[ix,jy,kz]=var3d[ix,jy,0]
        FOR fkz=1,(nz/2)-1 DO BEGIN
          result[ix,jy,kz]=result[ix,jy,kz]                                    $
                           +var3d[ix,jy,fkz]*COS(2*!DPI*fkz*zval)              $
                           +var3d[ix,jy,nz-fkz]*SIN(2*!DPI*fkz*zval)
        ENDFOR
      ENDFOR
    ENDFOR
  ENDFOR
  
  RETURN, result
END
