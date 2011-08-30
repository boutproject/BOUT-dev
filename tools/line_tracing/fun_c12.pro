FUNCTION fun_C12, x,y,z
; Right-hand side of equation d(ix)/d(iy)=C1(ix,iy,iz), d(ix)/d(iy)=C2(ix,iy,iz)

  COMMON griddata, g, deltaZtor, Ntor
  COMMON flags, flag, mc_flag

  ;; calculate RXY,BXY,BPXY,BTXY,SINTY,HTHE,BXCVX,BXCVZ
  CASE mc_flag OF
    'nearest'    : gvals=mc_close(x,y)
    'bilinear'   : gvals=mc_bilin(x,y)
    'bicubic'    : gvals=mc_bicube(x,y)
     else        : gvals=mc_close(x,y)
  ENDCASE
  RXY=gvals.rxy
  BXY=gvals.bxy
  BPXY=gvals.bpxy
  BTXY=gvals.btxy
  SINTY=gvals.sinty
  HTHE=gvals.hthe
  BXCVX=gvals.bxcvx
  BXCVZ=gvals.bxcvz

  ;;-derivatives of Apar in the index space
  CASE flag OF
    'trilin'      :   dApar=apar_trilin(x,y,z)
    'bilin_cube'  :   dApar=apar_bilin_cube(x,y,z)
    'bilin_fft'   :   dApar=apar_bilin_fft(x,y,z)
    'tricube'     :   dApar=apar_tricube(x,y,z)
    'bicube_fft'  :   dApar=apar_bicube_fft(x,y,z)
    'model'       :   dApar=apar_model(x,y,z)
     else         :   dApar=apar_trilin(x,y,z)
  ENDCASE
  Apar_val=dApar.f
  dApar_dx=dApar.x
  dApar_dy=dApar.y
  dApar_dz=dApar.z


  ;;-metric tensor components vs. radial & poloidal indices
  g_sup_11=(RXY*BPXY)^2
  g_sup_12=0.0d0
  g_sup_13=-SINTY*(RXY*BPXY)^2
  g_sup_22=1./HTHE^2
  g_sup_23=-(BTXY*HTHE/(BPXY*RXY))/(HTHE^2)
  g_sup_33=(SINTY*RXY*BPXY)^2 + (BXY/(RXY*BPXY))^2

  ;-subscripted coefs
  f_1 = - (1/BXY)*(g_sup_13*dApar_dx + g_sup_23*dApar_dy + g_sup_33*dApar_dz)
  f_3 =   (1/BXY)*(g_sup_11*dApar_dx + g_sup_13*dApar_dz)

  ;;-superscripted coefficients
  f_sup_1 = g_sup_11*f_1 + g_sup_13*f_3 + Apar_val*BXCVX ;(JPS) curvature term
  f_sup_3 = g_sup_13*f_1 + g_sup_33*f_3 + Apar_val*BXCVZ ;added so div(B~)=0

  ;-rhs of dynamic equations
  dx_dy = f_sup_1*HTHE/BPXY
  dz_dy = f_sup_3*HTHE/BPXY

  RETURN, [dx_dy, dz_dy]
END
