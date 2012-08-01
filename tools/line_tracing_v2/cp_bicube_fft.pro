PRO cp_bicube_fft, apar=apar
  ; Convert apar to closed periodic form for bicubic interpolation in x,y, FFT in z

  COMMON griddata, g, deltaZtor, Ntor

  sz=SIZE(apar)
  nx=sz[1] ; should be equal to BOUT NX, g.nx
  ny=sz[2] ; should be equal to BOUT NY, g.ny
  nz=sz[3] ; should be equal to BOUT MZ and should be power of 2!!!

  ; (JPS) Again in BOUT the y-index runs from 0->g.ny-1
  ;  NOTE: dy=2*!PI/g.ny, y[0]=0, y[g.ny-1]=2*!PI-dy
  ;  For third-order finite differences we can more y-indices but must
  ;  take into account the twist-shift condition
  apar_yp=FLTARR(nx,ny+3,nz)
  dz=deltaZtor/nz
  zper=ROUND(2*!DPI/deltaZtor)
  ; Take Fourier transform along y=0 coordinate index
  FOR ix=0,nx-1 DO BEGIN
    twist=g.shiftangle[ix]    
    fft_0=FFT(REFORM(apar[ix,0,*]),/DOUBLE)
    fft_1=FFT(REFORM(apar[ix,1,*]),/DOUBLE)
    fft_n=FFT(REFORM(apar[ix,ny-1,*]),/DOUBLE)
    FOR kz=0,nz-1 DO BEGIN
      z_val=kz*dz   ;initial value of z; at y=2*PI it is this +/-? twist
      ;calculate the value at y=2*PI of this z index using Fourier coefficients
      val_0=REAL_PART(fft_0[0])   ;value at y=2*PI
      val_1=REAL_PART(fft_1[0])   ;value at y=2*PI+dy
      val_n=REAL_PART(fft_n[0])   ;value at y=0-dy
      FOR fkz=1,(nz/2)-1 DO BEGIN
        val_0=val_0+2.*REAL_PART(fft_0[fkz])*COS(fkz*zper*(z_val+twist)) $
                   -2.*IMAGINARY(fft_0[fkz])*SIN(fkz*zper*(z_val+twist))
        val_1=val_1+2.*REAL_PART(fft_1[fkz])*COS(fkz*zper*(z_val+twist)) $
                   -2.*IMAGINARY(fft_1[fkz])*SIN(fkz*zper*(z_val+twist))
        val_n=val_n+2.*REAL_PART(fft_n[fkz])*COS(fkz*zper*(z_val-twist)) $
                   -2.*IMAGINARY(fft_n[fkz])*SIN(fkz*zper*(z_val-twist))
      ENDFOR
      apar_yp[ix,*,kz]=[val_n,REFORM(apar[ix,*,kz]),val_0,val_1]
    ENDFOR 
  ENDFOR

  ; (JPS) Now take the Fourier transform in z, over all x,y
  apar_fft=FFT_REAL(apar_yp)

  ; return the result
  apar=apar_fft

END
