; Integrate a function over y
FUNCTION int_y, var, mesh, loop=loop, nosmooth=nosmooth
  f = var
  
  s = SIZE(var, /dim)
  nx = s[0]
  loop = FLTARR(nx)
  
  status = gen_surface(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=period)
    
    ;IF period THEN BEGIN
    ;  ; Periodic - use FFT
    ;  f[xi,yi] = REAL_PART(fft_integrate(var[xi,yi], loop=lo))
    ;  loop[xi] = lo
    ;ENDIF ELSE BEGIN
      f[xi,yi] = int_func(var[xi,yi])
      IF NOT KEYWORD_SET(nosmooth) THEN BEGIN
        f[xi,yi] = SMOOTH(SMOOTH(f[xi,yi], 5, /edge), 5, /edge)
      ENDIF
      loop[xi] = f[xi,yi[N_ELEMENTS(yi)-1]] - f[xi,yi[0]]
    ;ENDELSE
  ENDREP UNTIL last
  
  RETURN, f
END
