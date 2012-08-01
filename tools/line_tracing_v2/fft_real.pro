;; Given a 3D variable [x,y,z] it calculates the series
;; A(x,y,z)=a_0(x,y)+SUM_OVER_K[ a~_k(x,y)*cos(kz)+b~_k(x,y)*sin(kz)]
;; a~_0(ix,jy) is stored in result[ix,jy,0]
;; a~_k(ix,jy) is stored in result[ix,jy,k]
;; b~_k(ix,jy) is stored in result[ix,jy,NK-k] 
;; the amplitudes for the Nyquist frequency are set to zero
FUNCTION fft_real, var3d
  ON_ERROR, 2
  
  s = SIZE(var3d, /dim)

  IF N_ELEMENTS(s) NE 3 THEN BEGIN
    PRINT, "Error: fft_real takes 3D variables"
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
      f = FFT(var3d[ix,jy,0:(nz-1)],/DOUBLE)
   ;f is an array of nz complex numbers corresponding to (a+i*b)exp(ikx)
   ;  f[0]  a_0 + i*b_0 for freq.=0   (b_0 should be 0)
   ;  f[k]  a_k + i*b_k for freq.=k
   ;f[nz/2] a_(nz/2)+i*b_(nz/2) for crit. Nyquist freq. (b should be 0)
   ;f[nz-k] a_(-k) + i*b_(-k) for freq.=-k

   ;For a real signal, we should have a_k=a_(-k), b_k=-b_(-k)
      result[ix,jy,0]=REAL_PART(f[0])    ;n=0 component
      result[ix,jy,nz/2]=0.0             ;don't use nyquist freq
      FOR kz=1,nz/2-1 DO BEGIN
        IF REAL_PART(f[kz]) NE REAL_PART(f[nz-kz]) THEN STOP
        IF IMAGINARY(f[kz]) NE -1.*IMAGINARY(f[nz-kz]) THEN STOP
        result[ix,jy,kz]=2.*REAL_PART(f[kz])     ;==a~ cosine coefficient
        result[ix,jy,nz-kz]=-2.*IMAGINARY(f[kz]) ;==b~ sine coefficient
        ;storing a_k in result[k], b_k in result[nz-k] for k=1,nz/2
      ENDFOR
    ENDFOR
  ENDFOR
  
  RETURN, result
END

