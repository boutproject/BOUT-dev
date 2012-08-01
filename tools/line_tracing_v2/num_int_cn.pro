PRO num_int_cn, x0, y0, z0, x1, y1, z1
; Perform integration from y0 to y1
; integration is done using Crank-Nicolson method
;============================================;

  COMMON griddata, g, deltaZtor, Ntor

  COMMON newton_data, x_0, y_0, z_0, h ;must be passed to newtfunc
  COMMON int_data, int_meth, nsteps

  V1=[x0,z0] ;;-dependent variables
  y=y0       ;;-independent variable
  delta_y=y1-y0
  ;;-use Nstep substeps
  Nstep=nsteps
  dy=delta_y/Nstep ;;=h

;;pass in the starting data to newtfunc
  x_0=V1[0]
  y_0=y
  z_0=V1[1]
  h=dy

  FOR i=0,Nstep-1 DO BEGIN
    ;;determine x_(n+1),z_(n+1) by solving the non-linear system
    V2=NEWTON(V1,'newtfunc',/double)
    V1 = V2
    ;; make sure that z remains in the range [0,deltaZtor) 
    V1[1]=FMODULO(V1[1],deltaZtor)

    ;;update time
    y = y + dy
    ;;pass in new data to newtfunc
    x_0=V1[0]
    y_0=y
    z_0=V1[1]
  ENDFOR

  x1=V2[0]
  z1=V2[1]
;
END
