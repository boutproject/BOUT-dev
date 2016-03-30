;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; BOUT / BOUT++ input grid file generator                         ;
;                                                                 ;
; NOTE: This code assumes periodic in Y - NO X-POINTS!            ;
;
; Input units are:
; Ni      - 10^20 m^-3
; Te, Te  - eV
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; concatenate two 2D arrays in the X direction
FUNCTION xconcat, f1, f2
  s1 = SIZE(f1, /dim)
  s2 = SIZE(f2, /dim)

  nx1 = s1[0]
  nx2 = s2[0]
  ny = s1[1]
  IF s2[1] NE ny THEN BEGIN
      PRINT, "Error in xconcat: Y size of the arrays does not match"
      STOP
  ENDIF
  
  f = DBLARR(nx1 + nx2, ny)
  
  FOR i=0, nx1-1 DO f[i,*] = f1[i,*]
  FOR i=0, nx2-1 DO f[i+nx1, *] = f2[i,*]

  RETURN, f
END

; returns the value with the maximum absolute value (+ve or -ve)
FUNCTION maxabs, var
  n = LONG(N_ELEMENTS(var))
  v = REFORM(var, n) ; turn into 1D vector

  ma = v[0]
  
  i = 1L
  REPEAT BEGIN
      IF ABS(v[i]) GT ABS(ma) THEN ma = v[i]
      i = i + 1L
  ENDREP UNTIL i EQ n

  RETURN, ma
END

; smooth a 2D function using SAVGOL in the x direction, FFT in Y
FUNCTION smooth_savgolfft, var, degree=degree, nmodes=nmodes, width=width, noedge=noedge
  IF NOT KEYWORD_SET(degree) THEN degree = 4
  IF NOT KEYWORD_SET(width) THEN width = degree
  
  s = SIZE(var, /dim)
  nx = s[0]
  ny = s[1]
  
  IF NOT KEYWORD_SET(nmodes) THEN nmodes = FIX( ny / 4 )

  farr = COMPLEXARR(nx, ny)

  ; take FFT in Y

  FOR x=0, nx-1 DO BEGIN
    farr[x, *] = FFT(var[x,*])
  ENDFOR

  ; smooth with savgol

  fres = COMPLEXARR(nx, ny)

  filter = SAVGOL(width, width, 0, degree)
  FOR y=0, nmodes DO BEGIN
    fres[*,y] = CONVOL(farr[*,y], filter, /edge_truncate)
    
    ;r = POLY_FIT(findgen(nx), REAL_PART(farr[*,y]), degree, yfit=rfit)
    ;i = POLY_FIT(findgen(nx), IMAGINARY(farr[*,y]), degree, yfit=ifit)
    ;fres[*,y] = COMPLEX(rfit, ifit)

    IF y NE 0 THEN fres[*,ny-y] = CONJ(fres[*,y])
  ENDFOR

  ; inverse FFT

  res = FLTARR(nx, ny)

  FOR x=0, nx-1 DO BEGIN
    res[x,*] = REAL_PART(FFT(fres[x,*], /inv))
  ENDFOR

  IF KEYWORD_SET(noedge) THEN BEGIN
    res = res[width:(nx-1-width), *]
  ENDIF

  RETURN, res
END

; rotate a 2D function in the y direction
FUNCTION roty, var, n
  s = size(var, /dimensions)
  nx = s[0]
  ny = s[1]

  v = DBLARR(nx, ny)

  FOR y=0, ny-1 DO BEGIN
      yp = (((y + n) MOD ny) + ny) MOD ny
      v[*,y] = var[*,yp]
  ENDFOR

  RETURN, v
END

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Routines used in flux-coordinate calculations
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

; calculate Grad P - JxB in flux coordinates
; NOTE: +ve = inward force
FUNCTION calc_force, psixy, Bpxy, Btxy, hthe, Rxy, dpdpsi
  MU = !PI * 4.e-7

  Bxy = SQRT(Bpxy^2 + Btxy^2)

  RETURN, DDX(psixy, Bxy^2 * hthe / Bpxy) - Btxy*Rxy*DDX(psixy, Btxy*hthe / (Rxy*Bpxy)) $
    + MU*hthe*dpdpsi / Bpxy
END

; Calculate poloidal field
FUNCTION calc_bp, psixy, Rxy, Zxy
  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]

  Bpxy = DBLARR(nx, ny)
  
  FOR y=0, ny-1 DO BEGIN
      dr = DERIV(Rxy[*,y])
      dz = DERIV(Zxy[*,y])
      dl = SQRT(dr^2 + dz^2)
      
      Bpxy[*,y] = DERIV(psixy[*,y]) / (dl * Rxy[*,y])
  ENDFOR

  RETURN, Bpxy
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
    
    f_all[*,y] = NEWTON(f_all[*,y], "f_newt")
    
    PRINT, y, "force imbalance: ", old_max_f, "-> ", MAX(ABS(f_newt(f_all[*,y])))
  ENDFOR
  
  Btold = Btxy
  
  FOR x=0, nx-1 DO BEGIN
    f = MEAN(f_all[x, *])
    
    Btxy[x,*] = f / Rxy[x,*]
  ENDFOR
  
  IF MIN(Btxy) LT 0.0 THEN Btxy = -Btxy
  
  PRINT, "Maximum change in Bt = ", max(abs(Btxy - Btold))
  PRINT, "Maximum percentage change = ", 100.*max(abs(Btxy - Btold) / Btold)
  Bxy = SQRT(Btxy^2 + Bpxy^2)
  
END

;; Calculate q
FUNCTION calc_q, Rxy, hthe, Btxy, Bpxy
  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]

  dtheta = 2.0*!PI / DOUBLE(ny)

  ; get local pitch
  pitch=hthe*Btxy/(Rxy*Bpxy)

  q = FLTARR(nx)
  ; perform loop integral
  FOR x=0, nx-1 DO BEGIN
    qinty = fft_integrate(pitch[x,*], loop=loop)
    q[x] = loop * dtheta
  ENDFOR

  RETURN, q / (2.0*!PI)
END

;; Correct hthe using force balance
FUNCTION correct_hthe, Rxy, psixy, Btxy, Bpxy, Bxy, hthe, dpdpsi, fixhthe=fixhthe
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
          CATCH, err
          IF err EQ 0 THEN BEGIN
            htmp = NEWTON(htmp, "new_hfunc")
          ENDIF
          CATCH, /CANCEL
            
          IF err GT 0 THEN BEGIN
            robust = 1   ; Next time 
            CATCH, err2
            IF err2 EQ 0 THEN htmp = NEWTON(htmp, "robust_hfunc")
            CATCH, /cancel
          ENDIF
          
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
          RETURN, hthe
      ENDIF
  ENDFOR

  ; NOTE: WEIRD ERROR IN FINAL POINT
  nh[0,ny-1] = 2.*nh[1,ny-1] - nh[2,ny-1] ; extrapolate from next two

  RETURN, nh
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Correct RBt and hthe using force balance and jpar
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION fj_newt, f
  COMMON fj_com, fixpos, Bt0, h0, r, Bp, mu0j, psi, mu0dp
  
  nx = N_ELEMENTS(f)/2 + 1 

  h1 = f[0:(nx-2)]
  Bt1 = f[(nx-1):*]

  IF fixpos EQ 0 THEN BEGIN
    h = [h0, h1]
    Bt = [Bt0, Bt1]
  ENDIF ELSE IF fixpos EQ nx-1 THEN BEGIN
    h = [h1, h0]
    Bt = [Bt1, Bt0]
  ENDIF ELSE BEGIN
    h = [h1[0:(fixpos-1)], h0, h1[fixpos:*]]
    Bt = [Bt1[0:(fixpos-1)], Bt0, Bt1[fixpos:*]]
  ENDELSE
  
  
  B = sqrt(Bt^2 + Bp^2)

  ; calculate force imbalance

  perr = (Bp*Bt*r/h)*DERIV(psi, Bt*h/(r*Bp)) - (Bp/h)*DERIV(psi, B*B*h/Bp) - mu0dp
  
  ; calculate difference in MU0*j

  jerr = (Bp*Bt*r / (B * h)) * ( DERIV(psi, B*B*h/Bp) - Bt*r*DERIV(h*Bt / (r*Bp)) ) - B*DERIV(Bt * r) - mu0j

  err = [perr[1:*], jerr[1:*]]

  RETURN, err
END

PRO correct_fh, Rxy, Btxy, Bpxy, hthe, psixy, dpdpsi, jpar, fix=fix
  COMMON fj_com, fixpos, Bt0, h0, r, Bp, mu0j, psi, mu0dp

  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]

  MU = 4.e-7*!PI

  IF NOT KEYWORD_SET(fix) THEN fix = 0
  fixpos = fix
  IF fixpos LT 0 THEN fixpos = 0
  IF fixpos GE nx THEN fixpos = nx-1
  
  mu0dp = dpdpsi * MU
  psi = psixy[*,0]

  ;!P.multi=[0,0,2,0,0]
  
  FOR y=0, ny-1 DO BEGIN
    
    Bt0 = Btxy[fixpos,y]
    h0  = hthe[fixpos,y]

    r = Rxy[*,y]
    Bp = Bpxy[*,y]
    mu0j = MU*jpar[*,y]
    
    IF fixpos EQ 0 THEN BEGIN
      ; fix the first point
      h1 = hthe[1:*,y]
      bt1 = Btxy[1:*,y]
    ENDIF ELSE IF fixpos EQ nx-1 THEN BEGIN
      ; fix last point
      h1 = hthe[0:(nx-2), y]
      bt1 = Btxy[0:(nx-2), y]
    ENDIF ELSE BEGIN
      ; fix somewhere in the middle  
      h1 = [hthe[0:(fixpos-1), y], hthe[(fixpos+1):*, y]]
      bt1 = [Btxy[0:(fixpos-1), y], Btxy[(fixpos+1):*, y]]
    ENDELSE

    f0 = [h1,Bt1]

    f1 = NEWTON(f0, "fj_newt")
    
    h1 = f1[0:(nx-2)]
    Bt1 = f1[(nx-1):*]
    
    IF fixpos EQ 0 THEN BEGIN
      h = [h0, h1]
      Bt = [Bt0, Bt1]
    ENDIF ELSE IF fixpos EQ nx-1 THEN BEGIN
      h = [h1, h0]
      Bt = [Bt1, Bt0]
    ENDIF ELSE BEGIN
      h = [h1[0:(fixpos-1)], h0, h1[fixpos:*]]
      Bt = [Bt1[0:(fixpos-1)], Bt0, Bt1[fixpos:*]]
    ENDELSE

    ;plot, hthe[*,y]
    ;oplot, h, psym=1
    ;plot, Btxy[*,y]
    ;oplot, Bt, psym=1
    
    ;cursor, i, j, /down

    hthe[*,y] = h
    Btxy[*,y] = Bt
  ENDFOR
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Relax equilibrium
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION grid_newt, data

  COMMON grid_newt_com, nx, ny, psixy, gs_f, dpdpsi, R0, Z0, min_f, yind
  
  IF yind LT 0 THEN BEGIN
      n = nx*ny
      pos = REFORM(data, nx, ny)
      
      Rxy = FLTARR(nx, ny)
      Zxy = Rxy
      
      FOR y=0, ny-1 DO BEGIN
          Rxy[*,y] = INTERPOL(R0[*,y], FINDGEN(nx), pos[*,y])
          Zxy[*,y] = INTERPOL(Z0[*,y], FINDGEN(nx), pos[*,y])
      ENDFOR
  ENDIF ELSE BEGIN
      ; just solving for one y index
      n = nx
      
      Rxy = R0
      Zxy = Z0

      Rxy[*,yind] = INTERPOL(R0[*,yind], FINDGEN(nx), data)
      Zxy[*,yind] = INTERPOL(Z0[*,yind], FINDGEN(nx), data)
  ENDELSE

  ; calculate Bpxy, Btxy and hthe

  Btxy = DBLARR(nx, ny)
  FOR x=0, nx-1 DO Btxy[x,*] = gs_f[x] / Rxy[x,*]
  hthe = calc_hthe(Rxy, Zxy)
  Bpxy = calc_bp(psixy, Rxy, Zxy)
  
  F = -1.0*calc_force(psixy, Bpxy, Btxy, hthe, Rxy, dpdpsi)

  fm = MAX(ABS(F))

  IF (fm LT min_f) OR (min_f LT 0.0) THEN BEGIN
      min_f = fm
      PRINT, MAX(ABS(Rxy - R0)), MAX(ABS(Zxy - Z0)), MAX(ABS(F))
  ENDIF

  ;!P.multi=[0,0,2,0,0]
  ;surface, Bpxy, chars=3
  ;surface, F, chars=3

  IF yind LT 0 THEN val = REFORM(F, n) ELSE val = F[*,yind]
  
  RETURN, val
END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Main BOUT_OUTPUT routine
;
; specifying /default asks no questions (non-interactive)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO bout_output, data, output=output, input=input,$
                 smooth=smooth, width=width, degree=degree, $
                 vacsmooth=vacsmooth, oldcurv=oldcurv, $
                 default=default, reverse=reverse, fixvals=fixvals

  safe_colors, /first

  IF KEYWORD_SET(output) THEN output=EXPAND_PATH(output)

  IF KEYWORD_SET(smooth) THEN BEGIN
    IF NOT KEYWORD_SET(width) THEN width = 10
    IF NOT KEYWORD_SET(degree) THEN degree = 3
  ENDIF

  IF NOT KEYWORD_SET(output) THEN output="bout.grd.pdb"

  MU = 4.0*!PI*1.0e-7

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Get the needed data from the structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  list = TAG_NAMES(data)

  IF in_list(list, "NX", /need) THEN nx    = data.nx
  IF in_list(list, "NY", /need) THEN ny    = data.ny
  IF in_list(list, "RXY", /need) THEN Rxy   = data.Rxy
  IF in_list(list, "ZXY", /need) THEN Zxy   = data.Zxy
  IF in_list(list, "BPXY", /need) THEN Bpxy  = data.Bpxy
  IF in_list(list, "BTXY", /need) THEN Btxy  = data.Btxy
  IF in_list(list, "PSIXY", /need) THEN psixy = data.psixy
  
  IF in_list(list, "PSI_AXIS") THEN psi_axis = data.psi_axis ELSE psi_axis = 0.0
  IF in_list(list, "PSI_BNDRY") THEN psi_bndry = data.psi_bndry ELSE psi_bndry = psixy[nx-1,0]

  dtheta = 2.0*!PI / DOUBLE(ny)

  got_qsafe = 1
  IF in_list(list, "QSAFE") THEN qsafe = data.qsafe ELSE BEGIN
      got_qsafe = 0
      PRINT, "WARNING: NO QSAFE SUPPLIED. WILL CALCULATE"
  ENDELSE
  
  blank = FLTARR(nx, ny)
  got_jpar = 1
  IF in_list(list, "JPAR") THEN Jpar  = data.Jpar ELSE BEGIN
      Jpar = blank
      got_jpar = 0
      PRINT, "Setting Jpar to zero"
  ENDELSE
  IF in_list(list, "TE") THEN Te    = data.Te ELSE BEGIN
      Te = blank
      PRINT, "Setting Te to zero"
  ENDELSE
  IF in_list(list, "TI") THEN Ti   = data.Ti ELSE BEGIN
      Ti = blank
      PRINT, "Setting Ti to zero"
  ENDELSE
  IF in_list(list, "NI") THEN Ni   = data.Ni ELSE BEGIN
      Ni = blank
      PRINT, "Setting Ni to zero"
  ENDELSE

  ;; Get pressure in Pa
  P = 1.602e-19*1.e20*Ni*(Te + Ti)

  PRINT, "Density "+STRTRIM(STRING(MIN(Ni)),2)+" -> "+STRTRIM(STRING(MAX(Ni)),2)
  PRINT, "Maximum pressure [Pa]: ", max(P)

  IF(max(P) LT 0.1) THEN BEGIN
      PRINT, "ERROR: NO PRESSURE!"
      RETURN
  ENDIF
  
  IF in_list(list, "DPDPSI") THEN dpdpsi = data.dpdpsi ELSE BEGIN
      PRINT, "WARNING: dpdpsi missing. Calculating from p and psixy"

      dpdpsi =DERIV(psixy[*,0], P[*,0])
  ENDELSE

  Bxy = SQRT(Btxy^2 + Bpxy^2)

  IF in_list(list, "RMAG") THEN rmag = data.rmag ELSE BEGIN
      rmag = MAX(ABS(Rxy))
      PRINT, "Setting rmag = ", rmag
  ENDELSE
  IF in_list(list, "BMAG") THEN bmag = data.bmag ELSE BEGIN
      bmag = MAX(ABS(Bxy))
      PRINT, "Setting bmag = ", bmag
  ENDELSE

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Add grid-points in the vacuum (optional)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT KEYWORD_SET(default) THEN BEGIN  ;; by default skip this

      if keyword_set(INPUT) then add_vac=input.add_vac else $
        add_vac=get_yesno("Add vacuum region?")

      IF add_vac THEN BEGIN
          done = 0
          REPEAT BEGIN
             ; get dpsi at the edge
          
              dpsi = psixy[nx-1, 0] - psixy[nx-2, 0]
              
              REPEAT BEGIN
                  nvac = get_integer("Number of vacuum surfaces to add:")
              ENDREP UNTIL nvac GT 0
              
              !P.multi=[0,0,1,0,0]
              
              vg = vacuum(REFORM(Rxy[nx-1, *]), $
                          REFORM(Zxy[nx-1, *]), $
                          REFORM(Bpxy[nx-1,*]), $
                          REFORM(Btxy[nx-1,*]), $
                          ABS(dpsi), n = nvac, sm=vacsmooth)    
              
              done = get_yesno("Is this ok?")
          ENDREP UNTIL done EQ 1
          
          Rxy = xconcat(Rxy, vg.r[1:*,*])
          Zxy = xconcat(Zxy, vg.z[1:*,*])
          Bpxy = xconcat(Bpxy, vg.bp[1:*,*])
          Btxy = xconcat(Btxy, vg.bt[1:*,*])
          
          nularr = FLTARR(nvac, ny)
          Jpar = xconcat(Jpar, nularr )
          Te = xconcat(Te, nularr+Te[nx-1, 0])
          Ti = xconcat(Ti, nularr+Ti[nx-1, 0])
          Ni = xconcat(Ni, nularr+Ni[nx-1, 0])
          
          dpdpsi = [dpdpsi, FLTARR(nvac)]
      
          parr = DBLARR(nvac, ny)
          FOR i=0, nvac-1 DO parr[i,*] = psixy[nx-1,0] + FLOAT(i+1)*dpsi
          psixy = xconcat(psixy, parr)
          
          nx = nx + nvac
      
      ENDIF
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; SMOOTHING
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF KEYWORD_SET(smooth) THEN BEGIN
    Rxy = smooth_savgolfft(Rxy, degree=degree, width=width, /noedge)
    Zxy = smooth_savgolfft(Zxy, degree=degree, width=width, /noedge)
    
    Btxy = smooth_savgolfft(Btxy, degree=degree, width=width, /noedge)
    Bpxy = smooth_savgolfft(Bpxy, degree=degree, width=width, /noedge)
    
    Jpar = smooth_savgolfft(Jpar, degree=degree, width=width, /noedge)
    
    psixy = psixy[width:(nx-1-width), *]
    dpdpsi = dpdpsi[width:(nx-1-width)]
    
    Te = Te[width:(nx-1-width), *]
    Ti = Ti[width:(nx-1-width), *]
    Ni = Ni[width:(nx-1-width), *]

    nx = nx - 2*width
  ENDIF

  ; make sure that R*Bt = const on flux surfaces
      
  FOR i=0, nx-1 DO BEGIN
      tmp = MEAN(Rxy[i,*]*Btxy[i,*])
      Btxy[i,*] = tmp / Rxy[i,*]
  ENDFOR

  Bxy = SQRT(Btxy^2 + Bpxy^2)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; SIGNS: SET EVERYTHING +VE FOR NOW
  ; 
  ; Will reverse Bt if needed later. Hence Bp is always
  ; positive
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  dp = DERIV(psixy[1:*,0])
  IF MIN(dp) LT 0.0 THEN BEGIN  ; make sure psi always increasing
      psixy = -psixy
  ENDIF

  dpsi = DBLARR(nx, ny)
  FOR i=0, ny-1 DO dpsi[*,i] = DERIV(psixy[*,i])
  
  ; 29/8/08 non-uniform mesh correction
  d2x = DBLARR(nx, ny)
  FOR i=0, ny-1 DO d2x[*,i] = DERIV(dpsi[*,i])

  IF MIN(Bpxy) LT 0.0 THEN Bpxy = -Bpxy
  IF MIN(Btxy) LT 0.0 THEN Btxy = -Btxy

  IF MAX(dpdpsi) GT ABS(MIN(dpdpsi)) THEN BEGIN ; dpdpsi must be negative
      dpdpsi = -dpdpsi
  ENDIF

  hthe = calc_hthe(Rxy, Zxy)

  IF NOT KEYWORD_SET(default) THEN BEGIN ; Skip this by default
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Correct values to make JxB, Jpar and q closer to original
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    PRINT, "Equilibrium correction options:"
    PRINT, "  0  No correction"
    PRINT, "  1  RBt using force balance"
    PRINT, "  2  hthe and RBt using force balance and q (FAILS)"
    PRINT, "  3  hthe and RBt using force balance and jpar"
    
    if keyword_set(INPUT) then copt=input.correct_opt else $
      copt = get_integer("Enter option:")
    
    IF copt EQ 1 THEN BEGIN
      ;; Correct f = RBt using force balance
      correct_f, Rxy, Zxy, dpdpsi, hthe, Btxy, Bpxy, Bxy, psixy
    ENDIF ELSE IF copt EQ 2 THEN BEGIN
      ;; Correct RBt and hthe using force balance and q

      IF got_qsafe THEN BEGIN
        
        horig = hthe
        btorig = Btxy
        
        !P.MULTI=[0,0,2,0,0]
        REPEAT BEGIN
          ; Use q to set RBt
          qnew = calc_q(Rxy, hthe, Btxy, Bpxy) 
          
          qdiff = MAX(ABS(qnew - qsafe))
          PRINT, "Maximum difference in q: ", qdiff
          
          plot, qsafe
          oplot, qnew, psym=1
          
          ; calculate hthe using force balance
          nh = correct_hthe(Rxy, psixy, Btxy, Bpxy, Bxy, hthe, dpdpsi, fixhthe=fixhthe)
          
          ; Use this difference to calculate f = RBt
          
          FOR x=0, nx-1 DO BEGIN
            Btxy[x,*] = Btxy[x,*] * qsafe[x] / qnew[x]
          ENDFOR
          
          hdiff = MAX(ABS(nh - hthe))
          
          PRINT, "Maximum change in hthe: ", hdiff
          hthe = nh
          
          plot, horig[*,FIX(ny/2)]
          oplot, nh[*,FIX(ny/2)], psym=1
          
          cursor, x, y, /down
        ENDREP UNTIL 0
      ENDIF ELSE BEGIN
        PRINT, " => Don't have qsafe. Can't do this correction"
      ENDELSE
    ENDIF ELSE IF copt EQ 3 THEN BEGIN
      ;; Correct hthe and RBt using force balance and jpar
    
      IF KEYWORD_SET(fixvals) AND got_qsafe THEN BEGIN
        ;; Calculate Bt based on q
        
        qnew = calc_q(Rxy, hthe, Btxy, Bpxy) 
        
        FOR x=0, nx-1 DO BEGIN
          Btxy[x,*] = Btxy[x,*] * qsafe[x] / qnew[x]
        ENDFOR
      ENDIF
      
      correct_fh, Rxy, Btxy, Bpxy, hthe, psixy, dpdpsi, jpar, fix=fixvals
    ENDIF
  ENDIF ; End of default

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; CALCULATE HTHE
  ; calculate an initial guess, then use force-balance
  ; equation to find hthe. Does not depend on signs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  PRINT, "Calculating poloidal arc length dthe"
 
  IF KEYWORD_SET(fixvals) AND got_qsafe THEN BEGIN
    ; Calculate hthe based on q
    
    qnew = calc_q(Rxy, hthe, Btxy, Bpxy) 
      
    IF MIN(qnew * qsafe) LT 0.0 THEN BEGIN
      qnew = -1.*qnew
    ENDIF

    FOR x=0, nx-1 DO BEGIN
      hthe[x,*] = hthe[x,*] * qsafe[x] / qnew[x]
    ENDFOR
  ENDIF

  ; calculate hthe such that in equilibrium
  nh = correct_hthe(Rxy, psixy, Btxy, Bpxy, Bxy, hthe, dpdpsi, fixhthe=fixvals)

  !P.multi=[0,0,1,0,0]
  
  PRINT, "Maximum difference in hthe: ", MAX(ABS(hthe - nh))
  PRINT, "Maximum percentage difference: ", 100.*MAX(ABS((hthe - nh)/hthe))
      
  PLOT, hthe[*,0], title="Poloidal arc length at midplane. line is initial estimate", color=1
  OPLOT, nh[*,0], psym=1, color=2
  OPLOT, nh[*,0], color=2
  
  ;PRINT, "CLICK TO CONTINUE"
  ;cursor, x, y, /down

  old_hthe = hthe

  IF NOT KEYWORD_SET(default) THEN BEGIN
    if keyword_set(INPUT) then use_new_hthe=input.use_new_hthe else $
      use_new_hthe=get_yesno("Use new hthe?")
    
    IF use_new_hthe THEN hthe = nh
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; CALCULATE PARALLEL CURRENT
  ; Provides a way to check if Btor should be reversed
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  PRINT, "Checking parallel current"

  j0 = ((Bpxy*Btxy*Rxy/(Bxy*hthe))*( DDX(psixy, Bxy^2*hthe/Bpxy) - Btxy*Rxy*DDX(psixy,Btxy*hthe/(Rxy*Bpxy)) ) $
        - Bxy*DDX(psixy, Btxy*Rxy)) / MU
  
  IF (MEAN(ABS(j0 + jpar)) LT MEAN(ABS(j0 - jpar))) OR KEYWORD_SET(reverse) THEN BEGIN
      PRINT, "****Equilibrium has -ve toroidal field"
      Btxy = -Btxy
      j0 = -j0
  ENDIF

  PRINT, "Maximum difference in jpar0: ", MAX(ABS(j0[2:*,*] - jpar[2:*,*]))
  PRINT, "Maximum percentage difference: ", 200.*MAX(ABS((jpar[2:*,*] - j0[2:*,*])/(jpar[2:*,*]+j0[2:*,*])))
      
  plot, jpar[*,0], yr=[MIN([MIN(jpar[*,0]), MIN(j0[*,0])]), MAX([MAX(jpar[*,0]), MAX(j0[*,0])])], $
    title="Parallel current at midplane. Black line is from equilibrium", color=1
  oplot, j0[*,0], psym=1, color=2
  oplot, j0[*,0], color=2
  
  IF got_jpar THEN BEGIN
      IF NOT KEYWORD_SET(default) THEN BEGIN ;; by default keep existing jpar

          if keyword_set(INPUT) then use_new_jpar=input.use_new_jpar else $
            use_new_jpar=get_yesno("Use new Jpar?")

          IF use_new_jpar THEN jpar = j0

      ENDIF
  ENDIF ELSE BEGIN
      PRINT, "Using this new Jpar"
      ;PRINT, "CLICK TO CONTINUE"
      ;cursor, x, y, /down
      jpar = j0
  ENDELSE

  ; ROTATE SO THAT y=0 is on inboard side

  rmin = MIN(Rxy[0,*], mid) ; get inboard midplane index

  Rxy = roty(Rxy, mid)
  Zxy = roty(Zxy, mid)
  Bpxy = roty(Bpxy, mid)
  Btxy = roty(Btxy, mid)
  Bxy = roty(Bxy, mid)
  dpsi = roty(dpsi, mid)
  hthe = roty(hthe, mid)
  Jpar = roty(jpar, mid)
  
  dy = DBLARR(nx, ny) + dtheta

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
  FOR i=0, nx-1 DO BEGIN
      qinty[i,*] = fft_integrate(pitch[i,*], loop=loop)*dtheta
      qloop[i] = loop * dtheta  ;; the loop integral
      sinty[i,*] = fft_integrate(dqdpsi[i,*])*dtheta
      pol_angle[i, *] = 2.0*!PI * qinty[i,*] / qloop[i]
  ENDFOR

  IF got_qsafe THEN BEGIN

      ShiftAngle = qsafe * 2.0*!PI ; Twist-shift angle
      ; zShift = qinty
      IF maxabs(qinty)*maxabs(ShiftAngle) LT 0.0 THEN BEGIN
          ; need to reverse ShiftAngle to be consistent
          PRINT, "q is negative. Reversing values from equilibrium"
          ShiftAngle = -ShiftAngle
      ENDIF

      plot, ShiftAngle/(2.0*!PI), title="Safety factor. Solid is from input q", color=1
      oplot, qloop/(2.0*!PI), psym=1, color=2
      
;       if keyword_set(INPUT) then begin
;           wait, 3
;       endif else begin
;           PRINT, "CLICK TO CONTINUE"
;           cursor, x, y, /down
;       endelse

      IF NOT KEYWORD_SET(default) THEN BEGIN ;; by default keep existing qsafe

          if keyword_set(INPUT) then use_new_qsafe=input.use_new_qsafe else $
            use_new_qsafe=get_yesno("Use new qsafe?")
          
            IF use_new_qsafe THEN ShiftAngle = qloop

      ENDIF
  ENDIF ELSE BEGIN
      qsafe = qloop
      PRINT, "Using loop integral for qsafe"
      
      plot, qsafe / (2.0*!PI), title="Safety factor", color=2

;       if keyword_set(INPUT) then begin
;           wait, 3
;       endif else begin
;           PRINT, "CLICK TO CONTINUE"
;           cursor, x, y, /down
;       endelse

      ShiftAngle = qsafe

      PRINT, "Using new qsafe"
  ENDELSE

  ;;;;;;;;;;;;;;;;;;;; THETA_ZERO ;;;;;;;;;;;;;;;;;;;;;;
  ; re-set zshift to be zero at the outboard midplane

  m = MAX(Rxy[0,*], mind) ; outboard midplane index

  qm = qinty[*,mind]
  sm = sinty[*,mind]

  FOR i=0, ny-1 DO BEGIN
     qinty[*,i] = qinty[*,i] - qm
     sinty[*,i] = sinty[*,i] - sm
  ENDFOR

  ;;;;;;;;;;;;;;;;;;;; CURVATURE ;;;;;;;;;;;;;;;;;;;;;;;
  ; Calculating b x kappa

  IF KEYWORD_SET(oldcurv) THEN BEGIN
      thetaxy = FLTARR(nx, ny)
      thetaxy[0,*] = 2.0*!PI*findgen(ny)/ny
      FOR i=1, nx-1 DO thetaxy[i,*] = thetaxy[0,*]

      brxy = FLTARR(nx, ny)
      bzxy = brxy

      FOR i=0, nx-1 DO BEGIN
          dr = fft_deriv(REFORM(Rxy[i,*]))
          dz = fft_deriv(REFORM(Zxy[i,*]))
          
          dl = sqrt(dr*dr + dz*dz)
          dr = dr / dl
          dz = dz / dl
          
          brxy[i,*] = bpxy[i,*]*dr
          bzxy[i,*] = bpxy[i,*]*dz
      ENDFOR
      
      ;STOP
      
      curvature, nx, ny, FLOAT(Rxy), FLOAT(Zxy), FLOAT(brxy), FLOAT(bzxy), FLOAT(btxy), FLOAT(psixy), FLOAT(thetaxy), $
        bxcv=bxcv

      bxcvx = bxcv.psi 
      bxcvy = bxcv.theta
      bxcvz = bxcv.phi - sinty*bxcv.psi - pitch*bxcv.theta

      ; x borders
      bxcvx[0,*] = bxcvx[1,*]
      bxcvx[nx-1,*] = bxcvx[nx-2,*]
      
      bxcvy[0,*] = bxcvy[1,*]
      bxcvy[nx-1,*] = bxcvy[nx-2,*]
      
      bxcvz[0,*] = bxcvz[1,*]
      bxcvz[nx-1,*] = bxcvz[nx-2,*]
      
  ENDIF ELSE BEGIN
      ; calculate in flux coordinates.
      
      dpb = DBLARR(nx, ny)      ; quantity used for y and z components
      
      FOR i=0, ny-1 DO BEGIN
          dpb[*,i] = MU*dpdpsi/Bxy[*,i]
      ENDFOR
      dpb = dpb + DDX(psixy, Bxy)
      
      bxcvx = Bpxy * DDY(Btxy*Rxy / Bxy) / hthe
      ;bxcvx = -Bpxy*Btxy*Rxy*DDY(Bxy) / (hthe*Bxy^2)
      bxcvy = Bpxy*Btxy*Rxy*dpb / (hthe*Bxy^2)
      bxcvz = -dpb - sinty*bxcvx 

  ENDELSE

  ;;;;;;;;;;;;;;;;;;; PLASMA PROFILES ;;;;;;;;;;;;;;;;;;;
  
  P = Ni * (Te + Ti)*1.602e-19*1.0e20

  IF MIN(P) LT 1.0e-4*MAX(P) THEN BEGIN
      PRINT, "****Minimum pressure is very small:", MIN(P)
      PRINT, "****Setting minimum pressure to 1e-4 of maximum"

      pmin = 1e-4*MAX(P)
      P = P + pmin

      nimin = 1e-2*MAX(Ni)
      Ni = Ni + nimin
      
      tnew = P / (Ni*1.602e-19*1.0e20) ; Te + Ti
      told = Te + Ti
      tmin = 0.5*pmin / (nimin*1.602e-19*1.0e20)
      
      Te = tnew * (Te + tmin) / (told + 2.*tmin)
      Ti = tnew * (Ti + tmin) / (told + 2.*tmin)
  ENDIF

  Ni_x = MAX(Ni)
  Ti_x = MAX(Ti)
  Te_x = MAX(Te)

  ; set indices
  ixseps1 = nx                  ; all inside separatrix
  ixseps2 = nx

  jyseps1_1 = -1                ; all in core
  jyseps2_2 = ny-1
  jyseps1_2 = ny/2
  jyseps2_1 = ny/2
  
  ; save to file

  handle = file_open(output, /CREATE)
  
  s = file_write(handle, "nx", nx)
  s = file_write(handle, "ny", ny)
      
  ;topology

  s = file_write(handle, "ixseps1", ixseps1)
  s = file_write(handle, "ixseps2", ixseps2)
  s = file_write(handle, "jyseps1_1", jyseps1_1)
  s = file_write(handle, "jyseps1_2", jyseps1_2)
  s = file_write(handle, "jyseps2_1", jyseps2_1)
  s = file_write(handle, "jyseps2_2", jyseps2_2)

  ; geometry
  
  s = file_write(handle, "dx", dpsi)
  s = file_write(handle, "dy", dy)

  s = file_write(handle, "d2x", d2x) ; non-uniform correction

  s = file_write(handle, "ShiftAngle", ShiftAngle)
  s = file_write(handle, "zShift", qinty)
  s = file_write(handle, "pol_angle", pol_angle)
  s = file_write(handle, "ShiftTorsion", dqdpsi)

  s = file_write(handle, "Rxy",  Rxy)
  s = file_write(handle, "Zxy",  Zxy)
  s = file_write(handle, "Bpxy", Bpxy)
  s = file_write(handle, "Btxy", Btxy)
  s = file_write(handle, "Bxy",  Bxy)
  s = file_write(handle, "hthe", hthe)
  s = file_write(handle, "sinty", sinty)
  s = file_write(handle, "psixy", psixy)
      
  ; plasma profiles

  s = file_write(handle, "pressure", P)
  s = file_write(handle, "Ni0", Ni)
  s = file_write(handle, "Te0", Te)
  s = file_write(handle, "Ti0", Ti)
  s = file_write(handle, "Ni_x", Ni_x)
  s = file_write(handle, "Te_x", Te_x)
  s = file_write(handle, "Ti_x", Ti_x)
  s = file_write(handle, "Jpar0", Jpar)
  s = file_write(handle, "bmag", bmag)
  s = file_write(handle, "rmag", rmag)
  
  s = file_write(handle, "bxcvx", bxcvx)
  s = file_write(handle, "bxcvy", bxcvy)
  s = file_write(handle, "bxcvz", bxcvz)
  
  s = file_write(handle, "psi_axis", psi_axis)
  s = file_write(handle, "psi_bndry", psi_bndry)

  file_close, handle

  PRINT, "DONE"
END
