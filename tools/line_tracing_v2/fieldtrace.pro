; Starting from point (xIn,zIn) at outer midplane make a 
; full poloidal circle and return to midplane at (xout,zout)
;
;   Inputs:
;    x0 [weber]
;    z0 [radian]
;   Outputs:
;    x3 [weber]
;    z3 [radian]
PRO FieldTrace, xIn=xIn, zIn=zIn, xout=xout, zout=zout, debug=debug

  COMMON int_data, int_meth, nsteps

  format="('{x:',1e9.3,', y:',1e9.3,', z:',1e9.3,'}')"

  ;(JPS) By extending y-indices to cover 0->2*!PI instead of 0->2*!PI-delta_y
  ;we can integrate along y closer to the twist-shift location
  epsilon=1/128.
  y_imid=0.0d
  y_omid=!DPI
  y_bcut=2*!DPI-epsilon

  ;;-Step 1: integrate from outer midplane to last point below the cut
  x0=xIn
  z0=zIn
  y0=y_omid
  y1=y_bcut
  IF int_meth EQ 'rk4' THEN num_int_rk4, x0, y0, z0, x1, y1, z1
  IF int_meth EQ 'cn'  THEN num_int_cn , x0, y0, z0, x1, y1, z1
  if keyword_set(DEBUG) then begin
      str_In = STRING([x0,y0,z0], f=format)
      str_out = STRING([x1,y1,z1], f=format)
      print, "Step1 " + str_In + " -> " + str_out
  endif

  ;;-Step 2, twist-shift jump between y_bcut and y_imid
  ;; from (x1,y1,z1) go to (x2,y2,z)
  y2=y_imid
  TwistShift, x1,z1, x2,z2, debug=debug
  if keyword_set(DEBUG) then begin
      str_In = STRING([x1,y1,z1], f=format)
      str_out = STRING([x2,y2,z2], f=format) 
      print, "Step2 " + str_In + " -> " + str_out
  endif

  ;;-Step 3, integrate from y2 to y3
  y2=y_imid
  y3=y_omid
  IF int_meth EQ 'rk4' THEN num_int_rk4, x2, y2, z2, x3, y3, z3
  IF int_meth EQ 'cn'  THEN num_int_cn , x2, y2, z2, x3, y3, z3
  if keyword_set(DEBUG) then begin
    str_In = STRING([x2,y2,z2], f=format)
    str_out = STRING([x3,y3,z3], f=format)
    print, "Step3 " + str_In + " -> " + str_out
  endif

  xout=x3
  zout=z3

  IF KEYWORD_SET(debug) THEN STOP
END
