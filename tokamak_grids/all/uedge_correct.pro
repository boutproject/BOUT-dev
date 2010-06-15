; Correct a UEDGE-style BOUT input
; i.e. the inputs used for BOUT-06 which can be used for 2fluid models
;
; Does the following:
;  o Calculates the pressure balance and puts a delta N and T into
;    the grid file to correct for this (used by BOUT++ only)
;  o Checks the q profile is consistent, and optionally corrects
;
;  o Calculates Jpar0 to be consistent with the mesh
;
; B.Dudson, University of York. September 2009
;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Routines for calculating hthe
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION calc_hthe, Rxy, Zxy
  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]

  dtheta = 2.0*!PI / DOUBLE(ny)
  
  hthe = DBLARR(nx, ny)
  
  FOR i=0, nx-1 DO BEGIN
      ;FOR j=1, ny-1 DO BEGIN
      ;    dr = Rxy[i, j] - Rxy[i, j-1]
      ;    dz = Zxy[i, j] - Zxy[i, j-1]
      ;    
      ;    hthe[i,j] = SQRT(dr^2 + dz^2) / dtheta
      ;ENDFOR
      ;
      ;hthe[i,0] = (hthe[i,1] + hthe[i,ny-1])/2.0

      dr = DERIV([Rxy[i, ny-1], REFORM(Rxy[i,*]), Rxy[i,0]])
      dz = DERIV([Zxy[i, ny-1], REFORM(Zxy[i,*]), Zxy[i,0]])
      
      dr = dr[1:ny]
      dz = dz[1:ny]
      
      ;dr = fft_deriv(Rxy[i-1,*])
      ;dz = fft_deriv(Zxy[i-1,*])

      hthe[i,*] = SQRT(dr^2 + dz^2) / dtheta
  ENDFOR

  RETURN, hthe
END

; more resilient version of DERIV
FUNCTION new_deriv, x, f
  n = N_ELEMENTS(f)

  df = (f[1:*] - f[0:n-2]) / (x[1:*] - x[0:n-2])
  
  df = [df[0], df]

  RETURN, df
END

FUNCTION new_hfunc, h
  COMMON hvars, psi, bbp, btr, nu_h, mudp, h0, fixpos

  IF fixpos EQ 0 THEN BEGIN
    h2 = [h0, h]
  ENDIF ELSE IF fixpos EQ N_ELEMENTS(psi)-1 THEN BEGIN
    h2 = [h, h0]
  ENDIF ELSE BEGIN
    h2 = [h[0:(fixpos-1)], h0, h[fixpos:*]]
  ENDELSE

  a = DERIV(psi, bbp*h2)
  b = -btr*DERIV(psi, nu_h*h2)
  c = h2*mudp

  f = a+b+c

  w = WHERE(h2 LT 0.0, count)

  FOR i=0, count-1 DO f[w[i]] = -10.*h2[w[i]]
  
  RETURN, f[1:*]
END

FUNCTION robust_hfunc, h
  COMMON hvars, psi, bbp, btr, nu_h, mudp, h0, fixpos

  IF fixpos EQ 0 THEN BEGIN
    h2 = [h0, h]
  ENDIF ELSE IF fixpos EQ N_ELEMENTS(psi)-1 THEN BEGIN
    h2 = [h, h0]
  ENDIF ELSE BEGIN
    h2 = [h[0:(fixpos-1)], h0, h[fixpos:*]]
  ENDELSE

  a = new_DERIV(psi, bbp*h2)
  b = -btr*new_DERIV(psi, nu_h*h2)
  c = h2*mudp

  f = a+b+c

  w = WHERE(h2 LT 0.0, count)

  FOR i=0, count-1 DO f[w[i]] = -10.*h2[w[i]]
  
  RETURN, f[1:*]
END

;; Correct hthe using force balance
FUNCTION correct_hthe, Rxy, psixy, Btxy, Bpxy, Bxy, hthe, dpdpsi, fixhthe=fixhthe, robust=robust
  COMMON hvars, xarr, bbp, btr, nu_h, mudp, h0, fixpos
  
  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]

  MU = 4.e-7*!PI

  dbtr = DBLARR(nx, ny)
  FOR i=0, ny-1 DO dbtr[*,i] = DERIV(psixy[*,i], Btxy[*,i]*Rxy[*,i])

  IF NOT KEYWORD_SET(fixhthe) THEN fixhthe = 0
  IF fixhthe LT 0 THEN fixhthe = 0
  IF fixhthe GT nx-1 THEN fixhthe = nx-1

  fixpos = fixhthe
  
  nh = DBLARR(nx, ny)
  nh[fixhthe,*] = hthe[fixhthe,*]
  FOR i=0, ny-1 DO BEGIN
      xarr = psixy[*,i]
      bbp = REFORM(Bxy[*,i]^2 / Bpxy[*,i])
      btr = REFORM(Btxy[*,i]*Rxy[*,i])
      nu_h = REFORM(Btxy[*,i]/(Rxy[*,i]*Bpxy[*,i]))
      mudp = MU*dpdpsi / REFORM(Bpxy[*,i])
      
      IF fixhthe EQ 0 THEN BEGIN
        ; use the original way and fix hthe on the inside
        
        h0 = hthe[0,i]
        
        IF KEYWORD_SET(robust) THEN BEGIN
          nh[1:*,i] = NEWTON(REFORM(hthe[1:*,i]), "robust_hfunc")
        ENDIF ELSE BEGIN
          nh[1:*,i] = NEWTON(REFORM(hthe[1:*,i]), "new_hfunc")
        ENDELSE
      ENDIF ELSE BEGIN
        ; Fix hthe in the middle of the grid
        
        h0 = hthe[fixhthe,i]
        
        IF fixhthe GE nx-1 THEN BEGIN
          ; fix last point
          htmp = hthe[0:(nx-2),i]
          fixhthe = nx-1
        ENDIF ELSE BEGIN
          ; fix somewhere in the middle  
          htmp = [hthe[0:(fixhthe-1),i], hthe[(fixhthe+1):*,i]]
        ENDELSE

        IF KEYWORD_SET(robust) THEN BEGIN
          htmp = NEWTON(htmp, "robust_hfunc")
        ENDIF ELSE BEGIN
          htmp = NEWTON(htmp, "new_hfunc")
        ENDELSE

        IF fixhthe GE nx-1 THEN BEGIN
          nh[0:(nx-2), i] = htmp
        ENDIF ELSE BEGIN
          nh[0:(fixhthe-1), i] = htmp[0:(fixhthe-1)]
          nh[(fixhthe+1):*, i]  = htmp[fixhthe:*]
        ENDELSE
      ENDELSE

      ; CHECK IF CONVERGED

      ;f = new_hfunc(REFORM(nh[1:*,i]))
      ;fmax = MAX(f)
      
      ;PRINT, i, fmax

      w = WHERE(nh[*,i] LT 0.0, count)
      IF count GT 0 THEN BEGIN
          PRINT, "Error in hthe solver: Negative solution at y = ", i
          STOP
      ENDIF
  ENDFOR

  ; NOTE: WEIRD ERROR IN FINAL POINT
  nh[0,ny-1] = 2.*nh[1,ny-1] - nh[2,ny-1] ; extrapolate from next two

  PRINT, "Maximum change in hthe = ", max(abs(nh - hthe))
  PRINT, "Maximum percentage change = ", 100*max(abs(nh - hthe) / hthe)

  RETURN, nh
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate f = R * Bt
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION f_newt, f
  COMMON f_newt_com, h, r, bp, rhs, psi

  bt = f / r
  
  b = sqrt(Bt^2 + Bp^2)

  RETURN, f * DERIV(psi, bt*h / (r*bp)) $
    - DERIV(psi, (B^2)*h / bp) $
    - rhs
END

;; Correct f = RBt by Newton iteration of radial force balance
PRO correct_f, Rxy, Zxy, dpdpsi, hthe, Btxy, Bpxy, Bxy, psixy
  COMMON f_newt_com, h, r, bp, rhs, pxy

  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]

  MU = 4.e-7*!PI

  PRINT, "Correcting f = R*Bt using radial force balance"

  hthe = calc_hthe(Rxy, Zxy)
  
  dpdpsi2d = DBLARR(nx, ny)
  FOR i=0,ny-1 DO dpdpsi2d[*,i] = dpdpsi
  
  rhs_all =  MU*hthe*dpdpsi2d / Bpxy
  
  f_all = Btxy * Rxy
  
  FOR y=0, ny-1 DO BEGIN
    h = hthe[*,y]
    r = Rxy[*,y]
    bp = Bpxy[*,y]
    rhs = rhs_all[*,y]
    pxy = psixy[*,y]
    
    old_max_f = MAX(ABS(f_newt(f_all[*,y])))
    
    f_all[*,y] = NEWTON(f_all[*,y], "f_newt", tolf=1e-6, itmax=1000)
    
    PRINT, y, "force imbalance: ", old_max_f, "-> ", MAX(ABS(f_newt(f_all[*,y])))
  ENDFOR
  
  Btold = Btxy
  
  FOR x=0, nx-1 DO BEGIN
    f = MEAN(f_all[x, *])
    
    Btxy[x,*] = f / Rxy[x,*]
  ENDFOR
  
  ;IF MIN(Btxy) LT 0.0 THEN Btxy = -Btxy
  ; Since sign doesn't matter, try to match input
  IF MEAN(ABS(Btxy + Btold)) LT MEAN(ABS(Btxy - Btold)) THEN Btxy = -Btxy
  
  
  PRINT, "Maximum change in Bt = ", max(abs(Btxy - Btold))
  PRINT, "Maximum percentage change = ", 100.*max(abs(Btxy - Btold) / Btold)
  Bxy = SQRT(Btxy^2 + Bpxy^2)
  
END

;; Calculate pressure
FUNCTION calc_pressure, ingrid, Bxy, Btxy, inhthe, title=title
  MU0 = 4.e-7*!PI

  nx = ingrid.nx
  ny = ingrid.ny

  j1_1 = ingrid.jyseps1_1
  j1_2 = ingrid.jyseps1_2
  j2_1 = ingrid.jyseps2_1
  j2_2 = ingrid.jyseps2_2
  
  pres = FLTARR(nx, ny)

  cfirst = 1
  FOR yind=0, ny-1 DO BEGIN
    ; Solve radial force-balance equation on each y slice
    
    psi = REFORM(ingrid.psixy[*,yind])
    
    B    = Bxy[*,yind]
    Bp   = ingrid.Bpxy[*,yind]
    Bt   = Btxy[*,yind]
    hthe = inhthe[*,yind]
    R    = ingrid.Rxy[*,yind]
    
    ; Calculate pressure gradient in field-aligned coordinates
    
    dpdx = ( Bt*R*DERIV(psi, Bt*hthe/(R*Bp)) - DERIV(psi, B*B*hthe/Bp) )*Bp/(MU0*hthe)
    
    ; Integrate to get pressure in Pascals
    
    p = int_func(psi, dpdx)
    
    ; If in the core, set outer edge to zero pressure. If in PF, set
    ; inner to zero.
    
    IF ((yind GT j1_1) AND (yind LE j1_2)) OR ((yind GT j2_1) AND (yind LE j2_2)) THEN BEGIN
      ; In the core
      p = p - p[nx-1]
      
    ENDIF ELSE BEGIN
      ; In PF / legs
      p = p - p[0]
    ENDELSE
    
    pres[*,yind] = p
    
  ENDFOR
  
  IF KEYWORD_SET(title) THEN BEGIN
    
    PLOT, pres[*,0], xtitle="X index", ytitle="Core pressure profiles at each y location [Pa]", chars=1.5, title=title, yr=[MIN(pres), MAX(pres)]
    FOR i=1, ny-1 DO OPLOT, pres[*,i]
  ENDIF

  RETURN, pres
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate Jpar, curvature
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; calculates x (psi) derivative for 2D variable
FUNCTION ddx, psi, var
  s = SIZE(var, /dimensions)
  nx = s[0]
  ny = s[1]

  dv = DBLARR(nx, ny)
  FOR i=0, ny-1 DO dv[*,i] = DERIV(psi[*,i], var[*,i])

  RETURN, dv
END

; calculates y (poloidal) derivatives.
; NOTE: For periodic domain only (NO X-POINTS)
FUNCTION ddy, var
  s = SIZE(var, /dimensions)
  nx = s[0]
  ny = s[1]
  
  dy = 2.0*!PI / ny

  result = DBLARR(nx, ny)
  
  FOR i=0, nx-1 DO BEGIN
      result[i,*] = fft_deriv(REFORM(var[i,*])) / dy
  ENDFOR

  RETURN, result
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; MAIN CODE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO uedge_correct, input, output, x0=x0, y0=y0, printps=printps
  ; Read the input grid
  ingrid = file_import(input)
  
  MU0 = 4.e-7*!PI

  IF KEYWORD_SET(printps) THEN BEGIN
    SET_PLOT, 'PS'
    DEVICE, file=printps, /landscape
  ENDIF

  ; Put topology info into more convenient variables
  
  nx = ingrid.nx
  ny = ingrid.ny

  j1_1 = ingrid.jyseps1_1
  j1_2 = ingrid.jyseps1_2
  j2_1 = ingrid.jyseps2_1
  j2_2 = ingrid.jyseps2_2

  ; Calculate the pressure according to the input grid
  pingrid = calc_pressure(ingrid, ingrid.Bxy, ingrid.Btxy, ingrid.hthe, title="Input grid")

  IF NOT KEYWORD_SET(printps) THEN CURSOR, x, y, /down

  ; Use the profile at y0 as the "correct" value
  IF NOT KEYWORD_SET(x0) THEN x0 = FIX(ingrid.nx/2)
  IF NOT KEYWORD_SET(y0) THEN y0 = ingrid.gjy0

  
  ; calculate profile pressure in Pa
  Pprof = 1.602e-19*1.e20*ingrid.Ni0*(ingrid.Te0 + ingrid.Ti0)
  oplot, Pprof[*,y0], lines=2
  
  IF get_yesno("Use fluid pressure? (no = magnetic) ") THEN BEGIN
    core_pres = Pprof[*,y0]
  ENDIF ELSE core_pres = pingrid[*,y0] ; Core pressure profile
  
  dpdpsi = DERIV(ingrid.psixy[*,y0], core_pres)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Correct f = RBt using force balance

  PRINT, "Correcting force balance by adjusting toroidal field"
  Btxy = ingrid.Btxy
  Bxy=ingrid.Bxy
  correct_f, ingrid.Rxy, ingrid.Zxy, dpdpsi, ingrid.hthe, Btxy, ingrid.Bpxy, Bxy, ingrid.psixy
  
  plot, Btxy[nx-1,*], lines=2, title="Toroidal field. Solid=input, dashed=new"
  oplot, ingrid.Btxy[nx-1,*]
  
  oplot, Btxy[0,*], lines=2
  oplot, ingrid.Btxy[0,*]

  IF NOT KEYWORD_SET(printps) THEN CURSOR, x, y, /down
 
  p1 = calc_pressure(ingrid, Bxy, Btxy, ingrid.hthe, title="Corrected Bt")
  
  IF NOT KEYWORD_SET(printps) THEN CURSOR, x, y, /down
  
  PRINT, "Correcting force balance by adjusting hthe"

  hthe = correct_hthe(ingrid.Rxy, ingrid.psixy, Btxy, ingrid.Bpxy, Bxy, $
                    ingrid.hthe, dpdpsi, fixhthe=x0, /robust)

  plot, ingrid.hthe[nx-1,*], lines=2, title="hthe. Solid=input, symbols=new"
  oplot, hthe[nx-1,*], psym=1
  
  oplot, ingrid.hthe[0,*]
  oplot, hthe[0,*], psym=1

  IF NOT KEYWORD_SET(printps) THEN CURSOR, x, y, /down

  p2 = calc_pressure(ingrid, Bxy, Btxy, hthe, title="Corrected Bt then hthe")

  IF NOT KEYWORD_SET(printps) THEN CURSOR, x, y, /down

  !P.Multi=[0,2,0,0,0]
  zr=[0.0,MAX([pingrid,p2])]
  surface, pingrid, az=-45, zr=zr, chars=1.5, xtitle="X", ytitle="Y", ztitle="Pressure [Pa]", title="Input grid"
  surface, p2, az=-45, zr=zr, chars=1.5, xtitle="X", ytitle="Y", ztitle="Pressure [Pa]", title="Corrected grid"
  !P.multi=[0,0,1,0,0]

  pressure = p2
  
  IF NOT KEYWORD_SET(printps) THEN CURSOR, x, y, /down
  
  Bpxy = ingrid.Bpxy
  psixy = ingrid.psixy
  Rxy = ingrid.Rxy

  ;; Need to re-calculate jpar0, qinty, sinty, qsafe, bxk etc.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Calculate Jpar

  jpar0 = ((Bpxy*Btxy*Rxy/(Bxy*hthe))*( DDX(psixy, Bxy^2*hthe/Bpxy) - Btxy*Rxy*DDX(psixy,Btxy*hthe/(Rxy*Bpxy)) ) $
           - Bxy*DDX(psixy, Btxy*Rxy)) / MU0

  ; interpolate last point and first point
  FOR i=0,nx-1 DO BEGIN
    jpar0[i,ny-1] = (2.*jpar0[i,ny-2] + jpar0[i,1])/3.
    jpar0[i,0] = (jpar0[i,ny-2] + 2.*jpar0[i,1])/3.
  ENDFOR

  !P.Multi=[0,2,0,0,0]
  zr=[MIN([jpar0,ingrid.jpar0]),MAX([jpar0,ingrid.jpar0])]
  surface, jpar0, zr=zr, az=-45, chars=1.5, ztitle="Parallel current [A/m^2]", xtitle="X", ytitle="Y", title="Input file"
  surface, ingrid.jpar0, zr=zr, az=-45, chars=1.5, ztitle="Parallel current [A/m^2]", xtitle="X", ytitle="Y", title="Corrected grid"
  !P.multi=[0,0,1,0,0]
  
  IF NOT KEYWORD_SET(printps) THEN CURSOR, x, y, /down

  dtheta = 2.0*!PI / DOUBLE(ny)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Calculate Safety factor
  
  ; calculate q
  pitch=hthe*Btxy/(Rxy*Bpxy)
  
  ; derivative with psi
  dqdpsi = DDX(psixy, pitch)

  ; perform loop integrals to get qinty
  PRINT, "Performing integrals"
  qinty = DBLARR(nx,ny)
  sinty = DBLARR(nx,ny)
  qloop = DBLARR(nx)
  pol_angle = DBLARR(nx, ny) ;; xi = qinty / q
  
  ; Core region: Use FFT to do loop integrals
  ; NOTE: NEED TO COPE WITH PRIVATE FLUX REGIONS ETC.
  FOR i=0, ingrid.ixseps1-1 DO BEGIN
    qinty[i,*] = fft_integrate(pitch[i,*], loop=loop)*dtheta
    qloop[i] = loop * dtheta  ;; the loop integral
    sinty[i,*] = fft_integrate(dqdpsi[i,*])*dtheta
    pol_angle[i, *] = 2.0*!PI * qinty[i,*] / qloop[i]

    ; Set qinty, sinty to zero at y0
    qinty[i,*] = qinty[i,*] - qinty[i,y0]
    sinty[i,*] = sinty[i,*] - sinty[i,y0]
  ENDFOR

  qsafe = qloop / (2.*!PI)

  plot, ingrid.qint1d/(2.*!PI), xtitle="X", ytitle="Safety factor", title="q solid=input grid, dashed=corrected", chars=1.5
  oplot, qsafe, lines=2

  IF NOT KEYWORD_SET(printps) THEN CURSOR, x, y, /down

  ;;;;;;;;;;;;;;;;;;;; CURVATURE ;;;;;;;;;;;;;;;;;;;;;;;
  ; Calculating b x kappa in field-aligned coordinates
  
  dpb = DBLARR(nx, ny)      ; quantity used for y and z components
      
  FOR i=0, ny-1 DO BEGIN
    dpb[*,i] = MU0*dpdpsi/Bxy[*,i]
  ENDFOR
  dpb = dpb + DDX(psixy, Bxy)
      
  bxcvx = Bpxy * DDY(Btxy*Rxy / Bxy) / hthe
  bxcvy = Bpxy*Btxy*Rxy*dpb / (hthe*Bxy^2)
  bxcvz = -dpb - sinty*bxcvx 

  ;;;;;;;;;;;;;;;;;;;;
  ;; Initial profiles
  
  ; calculate profile pressure in Pa
  Pprof = 1.602e-19*1.e20*ingrid.Ni0*(ingrid.Te0 + ingrid.Ti0)

  ; Set edge pressure
  pressure = pressure - MIN(pressure) +  MIN(Pprof)

  yr=[0,MAX([pressure, Pprof])]
  plot, Pprof[*,y0], ytitle="Pressure [Pa]", xtitle="X", title="Initial profiles=solid, JxB=dashed", yr=yr
  oplot, pressure[*,y0], lines=2
  
  pratio = pressure / Pprof

  print, "Ratio of pressures: ", MIN(pratio), MAX(pratio)

  ; Multiply all quantities by sqrt(pratio)
  Ni0 = ingrid.Ni0*sqrt(pratio)
  Ti0 = ingrid.Ti0*sqrt(pratio)
  Te0 = ingrid.Te0*sqrt(pratio)
  
  Ni_x = max(Ni0)
  Te_x = max(Te0)
  Ti_x = max(Ti0)

  IF KEYWORD_SET(printps) THEN BEGIN
    DEVICE, /close
    SET_PLOT, 'X'
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;
  ;; Write output grid

  fd = file_open(output, /create, /write)
  
  err = file_write(fd, "Bpxy", Bpxy)
  err = file_write(fd, "Btxy", Btxy)
  err = file_write(fd, "Bxy",  Bxy)
  
  err = file_write(fd, "Jpar0", Jpar0)
  err = file_write(fd, "Ni0", Ni0)
  err = file_write(fd, "Ni_x", Ni_x)
  err = file_write(fd, "R0", ingrid.R0)
  err = file_write(fd, "Rpsi", ingrid.Rpsi)
  err = file_write(fd, "Rthe", ingrid.Rthe)
  err = file_write(fd, "Rxy", Rxy)
  err = file_write(fd, "Te0", Te0)
  err = file_write(fd, "Te_x", Te_x)
  err = file_write(fd, "Theta", ingrid.Theta)
  err = file_write(fd, "Ti0", Ti0)
  err = file_write(fd, "Ti_x", Ti_x)
  err = file_write(fd, "VE0", ingrid.VE0)
  err = file_write(fd, "Vi0", ingrid.Vi0)
  err = file_write(fd, "Vi_x", ingrid.Vi_x)
  err = file_write(fd, "Zpsi", ingrid.Zpsi)
  err = file_write(fd, "Zthe", ingrid.Zthe)
  err = file_write(fd, "Zxy", ingrid.Zxy)
  err = file_write(fd, "bmag", ingrid.bmag)
  err = file_write(fd, "bxcvx", bxcvx)
  err = file_write(fd, "bxcvy", bxcvy)
  err = file_write(fd, "bxcvz", bxcvz)
  err = file_write(fd, "dlthe", ingrid.dlthe) ; Should be dy * hthe?
  err = file_write(fd, "dpsi", ingrid.dpsi)
  err = file_write(fd, "dx", ingrid.dx)
  err = file_write(fd, "dy", ingrid.dy)
  err = file_write(fd, "gjy0", ingrid.gjy0)
  err = file_write(fd, "hthe", hthe)
  err = file_write(fd, "hthe0", ingrid.hthe0)
  err = file_write(fd, "iNixnorm", ingrid.iNixnorm)
  err = file_write(fd, "ixlb2", ingrid.ixlb2)
  err = file_write(fd, "ixseps1", ingrid.ixseps1)
  err = file_write(fd, "ixseps2", ingrid.ixseps2)
  err = file_write(fd, "jyseps1_1", ingrid.jyseps1_1)
  err = file_write(fd, "jyseps2_1", ingrid.jyseps2_1)
  err = file_write(fd, "jyseps1_2", ingrid.jyseps1_2)
  err = file_write(fd, "jyseps2_2", ingrid.jyseps2_2)
  
  ; some kappa_ terms which aren't used (at least in BOUT++)
  
  err = file_write(fd, "nx", nx)
  err = file_write(fd, "ny", ny)
  err = file_write(fd, "phi0", ingrid.phi0)
  err = file_write(fd, "psixy", ingrid.psixy)
  
  err = file_write(fd, "q_safe", pitch)
  err = file_write(fd, "qint1d", qloop)
  err = file_write(fd, "qinty", qinty)
  err = file_write(fd, "sibdryg", ingrid.sibdryg)
  err = file_write(fd, "simagxg", ingrid.simagxg)
  err = file_write(fd, "sinty", sinty)
  err = file_write(fd, "x_array", ingrid.x_array)
  
  file_close, fd
  
  STOP
END



