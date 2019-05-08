; Integrate a function over y
FUNCTION int_y, var, mesh, loop=loop, nosmooth=nosmooth, simple=simple
  f = var
  
  s = SIZE(var, /dim)
  nx = s[0]
  loop = DBLARR(nx)
  
  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface_hypnotoad(last=last, xi=xi, period=period)

    IF period THEN BEGIN
       ; Add first point onto the end so wraps around for integration
       yi = [yi, yi[0]]
    ENDIF
    ; Integrate in 1D
    f1d = int_func(var[xi,yi], simple=simple)
    IF NOT KEYWORD_SET(nosmooth) THEN BEGIN
       f1d = SMOOTH(SMOOTH(f1d, 5, /edge_truncate), 5, /edge_truncate)
    ENDIF
    ny = N_ELEMENTS(f1d)
    loop[xi] = f1d[ny-1] - f1d[0]
    IF period THEN BEGIN
       ; Remove last point
       f1d = f1d[0:(ny-2)]
       yi = yi[0:(ny-2)]
    END
    f[xi,yi] = f1d
  ENDREP UNTIL last
  
  RETURN, f
END
