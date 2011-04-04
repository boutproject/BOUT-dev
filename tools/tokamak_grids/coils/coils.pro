; Calculate the Apar due to a given set of coils


; Convert polar to cartesian coordinates
FUNCTION polarToCart, vp
  RETURN, {x:vp.r*cos(vp.phi), y:vp.r*sin(vp.phi), z:vp.z}
END

; Convert cartesian to polar coordinates
FUNCTION cartToPolar, cp
  RETURN, {r:SQRT(cp.x^2 + cp.y^2), phi:atan(cp.y, cp.x), z:cp.z}
END

; Convert to cartesian
FUNCTION toCart, p
  CATCH, err
  IF err EQ 0 THEN BEGIN
    i = p.x
  ENDIF ELSE BEGIN
    CATCH, /cancel
    ; No 'x' component - convert 
    RETURN, polarToCart(p)
  ENDELSE
  ; Already in cartesian
  RETURN, p
END

; Convert to polar
FUNCTION toPolar, p
  CATCH, err
  IF err EQ 0 THEN BEGIN
    i = p.phi
  ENDIF ELSE BEGIN
    CATCH, /cancel
    ; No 'phi' component - convert 
    RETURN, cartToPolar(p)
  ENDELSE
  ; Already in polar
  RETURN, p
END

; Distance between two points
FUNCTION distance, p1, p2
  c1 = toCart(p1)
  c2 = toCart(p2)
  
  RETURN, SQRT((c1.x - c2.x)^2 + $
               (c1.y - c2.y)^2 + $
               (c1.z - c2.z)^2)
END

FUNCTION linediff, I, y
  COMMON linecom, c1, c2, c
  
  ; i between 0 and 1
  
  ; Current position
  ipos = {x:( i*c2.x + (1.-i)*c1.x ), $
          y:( i*c2.y + (1.-i)*c1.y ), $
          z:( i*c2.z + (1.-i)*c1.z )}
  
  d = distance(ipos, c)
  
  RETURN, [1./d]
END

; Wire from p1 to p2, carrying current 
; Get A at pos
FUNCTION AfromLine, p1, p2, current, pos, fast=fast 
  COMMON linecom, c1, c2, c
  
  ; Convert all coordinates to cartesian
  c1 = toCart(p1)
  c2 = toCart(p2)
  cin  = toCart(pos)
  
  len = distance(c1, c2) ; length of the wire
  
  Ivec = {x:current*(c2.x - c1.x)/len, $
          y:current*(c2.y - c1.y)/len, $
          z:current*(c2.z - c1.z)/len}
  
  IF KEYWORD_SET(fast) THEN BEGIN
    ; Use a fixed number of Simpson rule steps
    
    n = 2 ; Must be even
    h = len / FLOAT(n)
    FOR i=0, n DO BEGIN
      f = FLOAT(i) / FLOAT(n)
      ; Position along wire
      ipos = {x:( f*c2.x + (1.-f)*c1.x ), $
              y:( f*c2.y + (1.-f)*c1.y ), $
              z:( f*c2.z + (1.-f)*c1.z )}
      
      ; Distance
      d = distance(ipos, cin)
      
      IF i EQ 0 THEN BEGIN
        integral = 1. / d
      ENDIF ELSE IF i EQ n THEN BEGIN
        integral = integral + 1. / d
      ENDIF ELSE IF i MOD 2 EQ 1 THEN BEGIN
        integral = integral + 4. / d
      ENDIF ELSE BEGIN
        integral = integral + 2. / d
      ENDELSE
    ENDFOR
    integral = integral * (h / 3.) * 1.e-7
    result = {x:Ivec.x*integral, $
              y:Ivec.y*integral, $
              z:Ivec.z*integral}
  ENDIF ELSE BEGIN
    result = cin
    i = 0L
    REPEAT BEGIN
      c = {x:cin.x[i], y:cin.y[i], z:cin.z[i]}
      
      ; Integrate along the line
      a0 = [0.]
      a = LSODE(a0, 0., 1., 'linediff', lstat)
      a = a * 1.e-7 ; mu_0 / 4pi
      
      result.x[i] = Ivec.x*a[0]
      result.y[i] = Ivec.y*a[0]
      result.z[i] = Ivec.z*a[0]
      
      i = i + 1L
    ENDREP UNTIL i EQ N_ELEMENTS(cin.x)
  ENDELSE
  
  RETURN, result
END

FUNCTION AfromArc, p1, phi, current, pos, fast=fast 
  ; For now turn into a series of lines
  
  ps = p1
  pe = p1
  
  n = 10
  a = toCart(pos)
  
  dphi = phi / FLOAT(n)
  FOR i=0, n-1 DO BEGIN
    pe.phi = ps.phi + dphi
    
    a1 = AfromLine(ps, pe, current, pos, fast=fast)
    
    a.x = a.x + a1.x
    a.y = a.y + a1.y
    a.z = a.z + a1.z
    
    ps = pe
  ENDFOR
  
  RETURN, a
END

; Add cartesian vectors
FUNCTION addCart, a, b
  ac = toCart(a)
  bc = toCart(b)

  RETURN, {x:(ac.x+bc.x), y:(ac.y+bc.y), z:(ac.z+bc.z)}
END


; Add a coil, giving two corners
FUNCTION AfromCoil, p1, p2, current, pos, fast=fast
  
  ; Corners
  c0 = toPolar(p1)
  c2 = toPolar(p2)
  dphi = c2.phi - c0.phi
  
  c1 = {r:c0.r, z:c0.z, phi:c2.phi}
  c3 = {r:c2.r, z:c2.z, phi:c0.phi}
  
  A = AfromLine(c0, c1, current, pos, fast=fast)
  A = addCart(A, AfromLine(c1, c2, current, pos, fast=fast))
  A = addCart(A, AfromLine(c2, c3, current, pos, fast=fast))
  A = addCart(A, AfromLine(c3, c0, current, pos, fast=fast))
  
  RETURN, A
END

FUNCTION AfromCoilSet, set, current, pos, fast=fast
  
  shift = 2.*!PI / FLOAT(set.n)
  
  ; Go through the set of n coils
  FOR i=0, set.n-1 DO BEGIN
    c0 = {r:set.r1, z:set.z1, phi:(shift*i)}
    c1 = {r:set.r2, z:set.z2, phi:(c0.phi + set.dphi)}
    
    dA = AfromCoil(c0, c1, current*(-1.)^i, pos, fast=fast)
    IF i EQ 0 THEN A = dA ELSE A = addCart(A, dA)
  ENDFOR
  
  RETURN, A
END

; Take derivative in y, taking into account branch-cuts
FUNCTION ddy, var, mesh
  f = var

  dtheta = 2.*!PI / FLOAT(mesh.ny)

  status = gen_surface(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=period)
    IF period THEN BEGIN
       f[xi,yi] = fft_deriv(var[xi,yi])
    ENDIF ELSE f[xi,yi] = DERIV(var[xi,yi])
  ENDREP UNTIL last
  RETURN, f / dtheta
END

PRO coils, file, savefile=savefile, printps=printps, current=current, odd=odd

  IF NOT KEYWORD_SET(current) THEN current = 1.

  fast = 1
  
  ;;;; Define coil sets for MAST
  lower = {r1:1.311, z1:-0.791, r2:1.426, z2:-0.591, n:6, dphi:2.*!PI/12.}
  upper = {r1:1.311, z1:0.791, r2:1.426, z2:0.591, n:6, dphi:2.*!PI/12.}
  
  ;;;; Read in the grid file
  g = file_import(file)
  
  nz = 64
  dz = 2.*!PI / FLOAT(nz)
  
  ; Generate grid points
  r = FLTARR(g.nx, g.ny, nz)
  z = r
  phi = r
  FOR k=0,nz-1 DO BEGIN
    r[*,*,k] = g.Rxy
    z[*,*,k] = g.Zxy
    phi[*,*,k] = FLOAT(k)*dz
  ENDFOR
  
  pos = {r:r, z:z, phi:phi}
  A = AfromCoilSet(lower, current, pos, fast=fast)
  IF KEYWORD_SET(odd) THEN BEGIN
    A = addCart(A, AfromCoilSet(upper, current, pos, fast=fast))
  ENDIF ELSE BEGIN
    A = addCart(A, AfromCoilSet(upper, -current, pos, fast=fast))
  ENDELSE
  
  ; Convert to polar coordinates
  
  Ar   = A.x * COS(phi) + A.y * SIN(phi)
  Aphi = A.y * COS(phi) - A.x * SIN(phi)
  Az   = A.z
  
  ; Now need to convert to Apar. Get poloidal component
  Apol = FLTARR(g.nx, g.ny, nz)
  ; Loop over surfaces
  status = gen_surface(mesh=g) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=period)
    IF period THEN BEGIN
      ; periodic
      
      dr = fft_deriv(REFORM(g.Rxy[xi,yi]))
      dz = fft_deriv(REFORM(g.Zxy[xi,yi]))
    ENDIF ELSE BEGIN
      ; non-periodic
      dr = DERIV(REFORM(g.Rxy[xi,yi]))
      dz = DERIV(REFORM(g.Zxy[xi,yi]))
    ENDELSE
    
    dl = SQRT(dr^2 + dz^2)
    
    ; Get unit vector along surface
    dr = dr / dl
    dz = dz / dl
    
    ; Dot-product of dl with A
    FOR k=0, nz-1 DO BEGIN
      Apol[xi,yi,k] = dr * Ar[xi,yi,k] + dz * Az[xi,yi,k]
    ENDFOR
  ENDREP UNTIL last

  Apar = FLTARR(g.nx, g.ny, nz)
  
  FOR k=0, nz-1 DO BEGIN
    Apar[*,*,k] = (Apol[*,*,k] * g.Bpxy + Aphi[*,*,k] * g.Btxy) / g.Bxy
  ENDFOR
  
  nk = FIX(nz/2) - 1
  ; Take FFT of Apar
  Apar_k = FLTARR(g.nx, g.ny, 2*nk + 1)
  FOR x=0, g.nx-1 DO BEGIN
    FOR y=0, g.ny-1 DO BEGIN
      f = FFT(Apar[x,y,*])
      
      ; Put into BOUT++ input format
      Apar_k[x,y,0] = REAL_PART(f[0]) ; DC
      FOR k=0, nk-1 DO BEGIN
        ; Real then imaginary part of each 
        Apar_k[x,y,2*k+1] = REAL_PART(f[k+1])
        Apar_k[x,y,2*k+2] = IMAGINARY(f[k+1])
      ENDFOR
    ENDFOR
  ENDFOR
  
  ; Add this variable to the file
  
  f = file_open(file, /write)
  status = file_write(f, "rmp_A", Apar_k)
  file_close, f
  
  ; Plot the mode spectrum
  ergos_plot, apar, g, mode=3, /noshift
  
  ; Calculate psi component of B
  nx = g.nx
  ny = g.ny
  

  Bpsi = FLTARR(nx, ny, nz)
  
  dAdy = FLTARR(nx, ny, nz)
  dAdz = FLTARR(nx, ny, nz)
  
  ; Z derivative
  dz = 2.*!PI / FLOAT(nz)
  FOR i=0, nx-1 DO BEGIN
    FOR j=0, ny-1 DO BEGIN
      dAdz[i,j,*] = fft_deriv(apar[i,j,*]) / dz
    ENDFOR
  ENDFOR
  
  ; Theta derivative (branch-cuts)
  FOR k=0, nz-1 DO BEGIN
    pitch = g.hthe * g.Btxy / (g.Rxy * g.Bpxy) ; Field-line pitch
    dAdy[*,*,k] = DDY(REFORM(apar[*,*,k]), g) + pitch * dAdz[*,*,k]
  ENDFOR
  
  FOR k=0,nz-1 DO BEGIN  
    Bpsi[*,*,k] = (g.Btxy/g.hthe) * (dAdy[*,*,k] / g.Bxy) - ( g.Bxy/(g.Rxy*g.Bpxy) )*dAdz[*,*,k]
  ENDFOR
  
  ergos_plot, Bpsi, g, mode=3, /noshift

  ; Convert into contravariant component of B in PEST coordinates
  q = ABS(g.shiftangle)/(2.*!PI)
  f = g.Rxy * g.Btxy
  J = g.Rxy^2 / f
  FOR i=0,ny-1 DO J[*,i] = J[*,i] * q
  
  B1 = FLTARR(g.nx, g.ny, nz)
  FOR k=0,nz-1 DO B1[*,*,k] = Bpsi[*,*,k] * J

  ergos_plot, B1, g, mode=3, /noshift

  IF KEYWORD_SET(savefile) THEN SAVE, file=savefile
  
  IF KEYWORD_SET(printps) THEN BEGIN
    set_plot, 'PS'
    device, file=printps, /color
    
    ergos_plot, Ar, g, mode=3, /noshift, title="Ar", /rev
    ergos_plot, Az, g, mode=3, /noshift, title="Az", /rev
    ergos_plot, Aphi, g, mode=3, /noshift, title="Aphi", /rev
    ergos_plot, Apar, g, mode=3, /noshift, title="Apar", /rev
    ergos_plot, Bpsi, g, mode=3, /noshift, title="Bpsi from Apar", /rev
    ergos_plot, B1, g, mode=3, /noshift, title="B1", /rev
    device, /close
    set_plot, 'X'
  ENDIF
  
  STOP
END


