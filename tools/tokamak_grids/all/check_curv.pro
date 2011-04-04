; Check the curvature calculation

PRO check_curv, file
  g = pd_import(EXPAND_PATH(file))


  brxy = FLTARR(g.nx, g.ny)
  bzxy = brxy

  FOR i=0, g.nx-1 DO BEGIN
      dr = fft_deriv(REFORM(g.rxy[i,*]))
      dz = fft_deriv(REFORM(g.zxy[i,*]))

      dl = sqrt(dr*dr + dz*dz)
      dr = dr / dl
      dz = dz / dl

      brxy[i,*] = g.bpxy[i,*]*dr
      bzxy[i,*] = g.bpxy[i,*]*dz
  ENDFOR
  
  thetaxy = FLTARR(g.nx, g.ny)
  thetaxy[0,*] = 2.0*!PI*findgen(g.ny)/g.ny
  FOR i=1, g.nx-1 DO thetaxy[i,*] = thetaxy[0,*]

  curvature, g.nx, g.ny, g.rxy, g.zxy, brxy, bzxy, g.btxy, g.psixy, thetaxy, $
    bxcv=bxcv

  qsafe = g.Btxy * g.hthe / (g.Bpxy * g.Rxy)

  bxcvx = FLTARR(g.nx, g.ny)
  bxcvy = bxcvx
  bxcvz = bxcvx

  FOR i=0, g.nx-1 DO BEGIN
      FOR j=0, g.ny-1 DO BEGIN
          bxcvx[i,j] = bxcv[i,j].psi
          bxcvy[i,j] = bxcv[i,j].theta
          bxcvz[i,j] = bxcv[i,j].phi - g.sinty[i,j]*bxcv[i,j].psi - qsafe[i,j]*bxcv[i,j].theta
      ENDFOR
  ENDFOR

  ; x borders
  bxcvx[0,*] = bxcvx[1,*]
  bxcvx[g.nx-1,*] = bxcvx[g.nx-2,*]

  bxcvy[0,*] = bxcvy[1,*]
  bxcvy[g.nx-1,*] = bxcvy[g.nx-2,*]

  bxcvz[0,*] = bxcvz[1,*]
  bxcvz[g.nx-1,*] = bxcvz[g.nx-2,*]

  STOP
END
