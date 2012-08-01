FUNCTION newtfunc, x_vec
;functions g_1(x,z) and g_2(x,z) to be solved to update x,z

;Inputs:
; x_vec is a vector [x,z] containing the location to evaluate g_1 and g_2
;
;Common:
; need information for x_0 and z_0 as well as y_val=y+h/2
;==============================================================================;
  COMMON newton_data, x_0, y_0, z_0, h

  x=x_vec[0]
  y=y_0+h/2.0
  z=x_vec[1] 

  f_12=fun_c12((x_0+x)/2.0,y,(z_0+z)/2.0)

  RETURN, [ x-x_0-h*f_12[0],z-z_0-h*f_12[1] ]

END 
