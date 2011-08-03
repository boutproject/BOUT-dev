; Integrate a function over y
FUNCTION int_y, var, mesh, loop=loop, nosmooth=nosmooth, simple=simple
  f = var
  
  s = SIZE(var, /dim)
  nx = s[0]
  loop = FLTARR(nx)
  
  status = gen_surface(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=period)
    
    f[xi,yi] = int_func(var[xi,yi], simple=simple)
    IF NOT KEYWORD_SET(nosmooth) THEN BEGIN
      f[xi,yi] = SMOOTH(SMOOTH(f[xi,yi], 5, /edge_truncate), 5, /edge_truncate)
    ENDIF
    loop[xi] = f[xi,yi[N_ELEMENTS(yi)-1]] - f[xi,yi[0]]
  ENDREP UNTIL last
  
  RETURN, f
END
