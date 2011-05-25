; calculates x (psi) derivative for 2D variable
FUNCTION ddx, psi, var
  s = SIZE(var, /dimensions)
  nx = s[0]
  ny = s[1]

  dv = DBLARR(nx, ny)
  FOR i=0, ny-1 DO dv[*,i] = DERIV(psi[*,i], var[*,i])

  RETURN, dv
END
