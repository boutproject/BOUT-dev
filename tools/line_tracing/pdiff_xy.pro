FUNCTION pdiff_xy, f, xi, yi
  ; fxy - 2x2 [x,y] 
  ; xi,yi - indices in range [0-1]
  
  ; f = f0 + fx * x + fy * y + fxy * x * y
  
  A = TRANSPOSE([[1,1,1,1], [0,1,0,1], [0,0,1,1], [0,0,0,1]])
  fvals = [f[0,0], f[1,0], f[0,1], f[1,1]]
  
  SVDC, A,W,U,V
  res=SVSOL(U,W,V,fvals)
  
  f0 = res[0]
  fx = res[1]
  fy = res[2]
  fxy = res[3]
  
  val = f0 + fx*xi + fy*yi + fxy*xi*yi
  dx  = fx + fxy*yi
  dy  = fy + fxy*xi
  
  RETURN, {f:val, dfdx:dx, dfdy:dy}
END





