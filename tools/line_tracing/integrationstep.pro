PRO IntegrationStep, x0, y0, z0, x1, y1, z1
; Perform integration from y0 to y1
;============================================;

  COMMON griddata, g, deltaZtor, Ntor

  V1=[x0,z0] ;;-dependent variables
  y=y0       ;;-independent variable
  delta_y=y1-y0
  ;;-use RK4 with Nstep substeps
  Nstep=100
  dy=delta_y/Nstep

  FOR i=0,Nstep-1 DO BEGIN
    dvdy=DIFFERENTIAL(y,V1)
    V2=RK4(V1, dvdy, y, dy, 'differential', /DOUBLE)
    y = y + dy
    V1 = V2
    ;; (JPS) If we update z in fun_c12.pro, then we need to make sure it
    ;; remains in the range [0,deltaZtor) 
    V1[1]=FMODULO(V1[1],deltaZtor)
  ENDFOR
  x1=V2[0]
  z1=V2[1]
;
END
