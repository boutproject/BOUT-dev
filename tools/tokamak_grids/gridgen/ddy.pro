; Take derivative in y, taking into account branch-cuts
FUNCTION ddy, var, mesh
  f = var

  dtheta = 2.D*!DPI / DOUBLE(TOTAL(mesh.npol))

  status = gen_surface(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=period)
    IF period THEN BEGIN
       f[xi,yi] = fft_deriv(var[xi,yi])
    ENDIF ELSE f[xi,yi] = DERIV(var[xi,yi])
  ENDREP UNTIL last
  RETURN, f / dtheta
END
