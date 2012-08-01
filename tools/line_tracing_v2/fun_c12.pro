FUNCTION fun_C12, x,y,z
; Right-hand side of equation dx/dy=C1(x,y,z), dz/dy=C2(x,y,z)

  COMMON griddata, g, deltaZtor, Ntor
  COMMON flags, flag

  ;; calculate RXY,BXY,BPXY,BTXY,SINTY,HTHE,BXCVX,BXCVZ
  gvals=mc_bicube(x,y)
  RXY=gvals.rxy
  BXY=gvals.bxy
  BPXY=gvals.bpxy
  BTXY=gvals.btxy
  SINTY=gvals.sinty
  HTHE=gvals.hthe
  BXCVX=gvals.bxcvx
  BXCVY=gvals.bxcvy
  BXCVZ=gvals.bxcvz
  JPAR0=gvals.jpar0

  MU0=4*!DPI*10.^(-7)

  ;;-derivatives of Apar in the index space
  CASE flag OF
    'tricube'     :   dApar=apar_tricube(x,y,z)
    'bicube_fft'  :   dApar=apar_bicube_fft(x,y,z)
    'model'       :   dApar=apar_model(x,y,z)
     else         :   dApar=apar_tricube(x,y,z)
  ENDCASE
  Apar_val=dApar.f
  dApar_dx=dApar.x
  dApar_dy=dApar.y
  dApar_dz=dApar.z

  ;;Don't rely on numerical cancellation of metric coefficient terms, just use
  ;;analytic expressions for area elements in the equations instead
  A1   = (RXY*BPXY*BTXY)/(HTHE)         ; (R*B_theta)^2*nu/(h_theta)^2
  A2   = (BXY)^2                        ; (B_0)^2
  A3   = SINTY*A1                       ; I*(R*B_theta)^2*nu/(h_theta)^2
  J_   = (BPXY/HTHE)*MU0/(BXY^2)*JPAR0  ; Parallel current term

  b0dgy = BPXY/HTHE
  bdgx  = (1/BXY)*( (-A1)*dApar_dy + (-A2)*dApar_dz ) + Apar_val* BXCVX
  bdgy  = (1/BXY)*( ( A1)*dApar_dx + (-A3)*dApar_dz ) + Apar_val*(BXCVY + J_)
  bdgz  = (1/BXY)*( ( A2)*dApar_dx + ( A3)*dApar_dy ) + Apar_val* BXCVZ

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  OLD METHOD IS BELOW, BUT QUALITATIVE DIFFERENCES ARE OBSERVED IN PLOTS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  ;;-metric tensor components vs. radial & poloidal indices
;  g_sup_11=(RXY*BPXY)^2
;  g_sup_12=0.0d
;  g_sup_13=-SINTY*(RXY*BPXY)^2
;  g_sup_22=1./HTHE^2
;  g_sup_23=-(BTXY*HTHE/(BPXY*RXY))/(HTHE^2)
;  g_sup_33=(SINTY*RXY*BPXY)^2 + (BXY/(RXY*BPXY))^2
;
;  ;-subscripted coefs
;  f_1 = - (1/BXY)*(g_sup_13*dApar_dx + g_sup_23*dApar_dy + g_sup_33*dApar_dz)
;  f_3 =   (1/BXY)*(g_sup_11*dApar_dx + g_sup_13*dApar_dz)
;
;  ;;-superscripted coefficients
;  b0dgy = BPXY/HTHE
;  bdgx = g_sup_11*f_1 + g_sup_13*f_3 + Apar_val*BXCVX
;  bdgy = g_sup_12*f_1 + g_sup_23*f_3 + Apar_val*BXCVY +                        $
;          Apar_val*(BPXY/HTHE)*MU0/(BXY^2)*JPAR0
;  bdgz = g_sup_13*f_1 + g_sup_33*f_3 + Apar_val*BXCVZ
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;-rhs of dynamic equations
  dx_dy = bdgx/(b0dgy+bdgy)
  dz_dy = bdgz/(b0dgy+bdgy)

  RETURN, [dx_dy, dz_dy]
END
