PRO cp_bilin_fft, apar=apar
  ; Convert apar to closed periodic form for bilinear interpolation in x,y
  ; with a Fourier transform in z

  COMMON griddata, g, deltaZtor, Ntor

  sz=SIZE(apar)
  nx=sz[1] ; should be equal to BOUT NX, g.nx
  ny=sz[2] ; should be equal to BOUT NY, g.ny
  nz=sz[3] ; should be equal to BOUT MZ and should be power of 2!!!

  ; (JPS) Again in BOUT the y-index runs from 0->g.ny-1
  ;  NOTE: dy=2*!PI/g.ny, y[0]=0, y[g.ny-1]=2*!PI-dy
  ;  For first-order finite differences we can add another y-index but must
  ;  take into account the twist-shift condition
  apar_yp=FLTARR(nx,ny+1,nz)
  dz=deltaZtor/nz
  zper=ROUND(2*!DPI/deltaZtor)
  ; Take Fourier transform along y=0 coordinate index
  FOR ix=0,nx-1 DO BEGIN
    fft_0=FFT(REFORM(apar[ix,0,*]))
    twist=g.shiftangle[ix]    
    FOR kz=0,nz-1 DO BEGIN
      z_val=kz*dz   ;initial value of z; at y=2*PI it is this +/-? twist
      ;calculate the value at y=2*PI of this z index using Fourier coefficients
      value=REAL_PART(fft_0[0])
      FOR fkz=1,(nz/2)-1 DO BEGIN
        value=value+2.*REAL_PART(fft_0[fkz])*COS(fkz*zper*(z_val+twist)) $
                   -2.*IMAGINARY(fft_0[fkz])*SIN(fkz*zper*(z_val+twist))
      ENDFOR
      apar_yp[ix,*,kz]=[REFORM(apar[ix,*,kz]),value]
    ENDFOR 
  ENDFOR

  ; (JPS) Now take the Fourier transform in z, over all x,y
  apar_fft=FFT_REAL(apar_yp)

  ; return the result
  apar=apar_fft

END
