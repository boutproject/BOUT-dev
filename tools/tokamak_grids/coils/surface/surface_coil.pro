; Implements an idealised surface current to model external field
; perturbations

FUNCTION loop_integrate, f
  RETURN, TOTAL(f)
END

; return the 2D rms of a 3D variable
FUNCTION zrms, f
  s = SIZE(f, /dim)

  nx = s[0]
  ny = s[1]
  nz = s[2]

  r = FLTARR(nx, ny)

  FOR x=0, nx-1 DO BEGIN
    FOR y=0, ny-1 DO BEGIN
      r[x,y] = SQRT(MEAN(f[x,y,*]^2))
    ENDFOR
  ENDFOR

  RETURN, r
END

PRO surface_coil, grid
  ; read in the grid file
  g = pd_import(EXPAND_PATH(grid)) 
  
  s = SIZE(g.Rxy, /dim)

  nx = s[0]
  ny = s[1]
  PRINT, "Size of grid: ", nx, ny

  y0 = FIX(0.3*ny)
  y1 = FIX(0.7*ny)

  safe_colors, /first

  PRINT, "Setting poloidal location of coils"
  REPEAT BEGIN
    PRINT, "Coil index range ", y0, y1
    
    plot, g.Rxy, g.Zxy, psym=3, /iso, title="Grid", color=1
   
    oplot, g.Rxy[nx-1,y0:y1], g.Zxy[nx-1,y0:y1], thick=3, color=2
    
    done = get_yesno("Is this ok?")
    
    IF NOT done THEN BEGIN
      y0 = get_float("Enter y0: ")
      y1 = get_float("Enter y1: ")

      y0 = y0 * ny
      y1 = y1 * ny

      IF y0 LT 0 THEN y0 = 0
      IF y0 GE ny-1 THEN y0 = ny-2
      IF y1 LE y0 THEN y1 = y0+1
      IF y1 GE ny THEN y1 = ny-1
    ENDIF
  ENDREP UNTIL done
  
  REPEAT BEGIN
    nnodes = get_integer("Number of nodes in theta: ")
    n_mode = get_integer("Toroidal mode number n  : ")
    
    IF nnodes LT 1 THEN nnodes = 1

    m_mode = FLOAT(nnodes + 1) / 2.0

    Jzeta = FLTARR(ny)
    Jthe  = FLTARR(ny)

    tmp = FLTARR(ny)

    FOR y=y0, y1 DO BEGIN
      Jzeta[y] = sin(2.0*!PI*m_mode*(y-y0)/(y1-y0)) / (g.hthe[nx-1,y]*g.Rxy[nx-1,y])
      tmp[y] = (g.hthe[nx-1,y] / (g.Rxy[nx-1,y] * g.Bpxy[nx-1,y])) * Jzeta[y]
    ENDFOR

    ; Make sure there's no net toroidal current in the coils
    ; this needs to be done for an even number of nodes

    w = WHERE(jzeta GT 0.0)
    factor = 1.0
    REPEAT BEGIN
      Jzeta[w] = Jzeta[w] * factor
      FOR y=y0, y1 DO tmp[y] = (g.hthe[nx-1,y] / (g.Rxy[nx-1,y] * g.Bpxy[nx-1,y])) * Jzeta[y]

      Jthe = REAL_PART(fft_integrate(tmp, loop=jtloop) * 2.0*!PI / FLOAT(ny))
      Jthe = Jthe - Jthe[y0]
      
      tot = fft_integrate(ABS(tmp), loop=jtnorm) * 2.0*!PI / FLOAT(ny)
      
      jtloop = REAL_PART(jtloop)
      jtnorm = REAL_PART(jtnorm)

      factor = 1.0 - jtloop / jtnorm

      ;PRINT, "loop = ", jtloop, " factor = ", factor
      ;plot, jzeta, color=1
      ;cursor, a,b, /down
    ENDREP UNTIL ABS(jtloop) LT 1.0e-5

    IF y0 NE 0 THEN Jthe[0:y0] = 0.0
    IF y1 NE ny-1 THEN Jthe[y1+1:*] = 0.0

    plot, Jzeta, color=1, yr=[min([jzeta, jthe]), max([jzeta, jthe])]
    oplot, Jthe, color=2
    cursor, a, b, /down
    
    ; Multiply Jthe by the other normalisation quantities
    Jthe = -1.* Jthe * n_mode * REFORM(g.Bpxy[nx-1,*])
    
    ; Get Jr and Jz from Jthe (convert to cylindrical coords)
    
    dr = fft_deriv(REFORM(g.Rxy[nx-1,*]))
    dz = fft_deriv(REFORM(g.Zxy[nx-1,*]))
    
    dl = SQRT(dr*dr + dz * dz) ; get unit vector
    dr = dr / dl
    dz = dz / dl

    Jr = Jthe * dr / REFORM(g.hthe[nx-1,*])
    Jz = Jthe * dz / REFORM(g.hthe[nx-1,*])

    plot, Jzeta, color=1, yr=[min([jzeta, jthe]), max([jzeta, jthe])], $
      title="Surface current: Zeta (black), R (red), Z (blue)"
    oplot, Jr, color=2
    oplot, Jz, color=4
    cursor, a, b, /down

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Calculate A. For now a brute-force method in cartesian coords

    nz = 8*n_mode ; number of points to use toroidally
    dzeta = 2.*!PI / FLOAT(nz)

    ; display surface current

    Jd_z = FLTARR(nz, ny)
    Jd_t = Jd_z

    FOR z=0, nz-1 DO BEGIN
      zeta = 2.0*!PI * FLOAT(z) / FLOAT(nz)
      
      Jd_z[z,*] = Jzeta * sin(zeta)
      Jd_t[z,*] = Jthe * cos(zeta)
    ENDFOR
    
    !P.color = 1
    PLOT_FIELD, Jd_z, Jd_t, $
      title="Vessel surface current (1/n th of torus)", n=500
    
    cursor, x, y, /down

    ; cartesian components on the surface
    Jx2 = FLTARR(ny, nz)
    Jy2 = Jx2
    Jz2 = Jx2

    ; position (x,y,z)
    xs = Jx2
    ys = xs
    zs = xs

    FOR z=0, nz-1 DO BEGIN
      zeta = 2.0*!PI * FLOAT(z) / FLOAT(nz)

      rs = (2.*g.Rxy[nx-1,*] - g.Rxy[nx-2,*])

      Jx2[*,z] = Jr * cos(n_mode * zeta) * cos(zeta) $
        - Jzeta * sin(n_mode * zeta) * sin(zeta) / rs

      Jy2[*,z] = Jr * cos(n_mode * zeta) * sin(zeta) $
        + Jzeta * sin(n_mode * zeta) * cos(zeta) / rs

      Jz2[*,z] = Jz * cos(n_mode * zeta)

      ; put the coils just outside the domain
      xs[*,z] = rs * cos(zeta)
      ys[*,z] = rs * sin(zeta)
      zs[*,z] = 2.*g.Zxy[nx-1,*] - g.Zxy[nx-2,*]
    ENDFOR
    
    
    PRINT, "Calculating vector potential"

    Ax = FLTARR(nx, ny, nz)
    Ay = Ax
    Az = Ax

    PRINT, ""

    FOR x=0, nx-1 DO BEGIN
      FOR y=0, ny-1 DO BEGIN
        FOR z=0, nz-1 DO BEGIN
          ; Calculate A at this point

          zeta = 2.0*!PI * FLOAT(z) / FLOAT(nz)
          
          ; get cartesian position of this point
          xp = g.Rxy[x,y] * cos(zeta)
          yp = g.Rxy[x,y] * sin(zeta)
          zp = g.Zxy[x,y]
          
          dx = xs - xp
          dy = ys - yp
          dz = zs - zp

          ; distance to coils
          r = SQRT( dx*dx + dy*dy + dz*dz )
          
          ; Need to do a double integral over theta and zeta
          
          tmpx = Jx2 / r
          tmpy = Jy2 / r
          tmpz = Jz2 / r
          
          tmpx2 = FLTARR(nz)
          tmpy2 = tmpx2
          tmpz2 = tmpx2
          
          FOR z2=0, nz-1 DO BEGIN
            s = REFORM(g.hthe[nx-1,*]*g.Rxy[nx-1,*]) ; surface element
            tmpx[*,z2] = tmpx[*,z2] * s
            tmpy[*,z2] = tmpy[*,z2] * s
            tmpz[*,z2] = tmpz[*,z2] * s
            
            ; do poloidal integral

            ;a = fft_integrate(tmpx[*,z2], loop=loop)
            ;tmpx2[z2] = loop
            tmpx2[z2] = loop_integrate(tmpx[*,z2])

            ;a = fft_integrate(tmpy[*,z2], loop=loop)
            ;tmpy2[z2] = loop
            tmpy2[z2] = loop_integrate(tmpy[*,z2])
            
            ;a = fft_integrate(tmpz[*,z2], loop=loop)
            ;tmpz2[z2] = loop
            tmpz2[z2] = loop_integrate(tmpz[*,z2])
          ENDFOR
          
          ; toroidal integral

          ;a = fft_integrate(tmpx2, loop=loop)
          ;Ax[x,y,z] = loop
          
          Ax[x,y,z] = loop_integrate(tmpx2)

          ;a = fft_integrate(tmpy2, loop=loop)
          ;Ay[x,y,z] = loop

          Ay[x,y,z] = loop_integrate(tmpy2)
          
          ;a = fft_integrate(tmpz2, loop=loop)
          ;Az[x,y,z] = loop

          Az[x,y,z] = loop_integrate(tmpz2)
        ENDFOR
      ENDFOR
       WRITEU, -1, 13, "Calculated x = "+STRING(x)+" of "+STRING(nx)
    ENDFOR
    
    ; Calculate Ar and Athe
    
    Ar = FLTARR(nx, ny, nz)
    A_zeta = Ar
    
    FOR z=0, nz-1 DO BEGIN
      zeta = 2.0*!PI * FLOAT(z) / FLOAT(nz)
      
      Ar[*,*,z] = Ax[*,*,z] * cos(zeta) + Ay[*,*,z] * sin(zeta)
      A_zeta[*,*,z] = Ay[*,*,z] * cos(zeta) - Ax[*,*,z]*sin(zeta)
    ENDFOR

    FOR y=0, ny-1 DO A_zeta[*,y,*] = A_zeta[*,y,*] * g.Rxy[nx-1,y]

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Have A in cylindrical coords (Ar, Az, Azeta)
    ; Put into toroidal coords (Apsi, Athe, Azeta)
    
    A_the = Ar
    A_psi = Ar

    FOR x=0, nx-2 DO BEGIN
      dr = fft_deriv(REFORM(g.Rxy[x,*]))
      dz = fft_deriv(REFORM(g.Zxy[x,*]))
      
      dl = SQRT(dr*dr + dz * dz) ; get unit vector
      dr = dr / dl
      dz = dz / dl

      ; (dr, dz) is unit vector in flux surface

      FOR y=0, ny-1 DO BEGIN
        A_the[x,y,*] = (Ar[x,y,*] * dr[y] + Az[x,y,*] * dz[y]) * g.hthe[x,y]
        
        A_psi[x,y,*] = (Ar[x,y,*] * dz[y] - Az[x,y,*] * dr[y]) / (g.Rxy[x,y] * g.Bpxy[x,y])
      ENDFOR
      
    ENDFOR

    ; calculate A^psi, A^the, A^zeta
    Apsi = A_psi
    Athe = A_the
    Azeta = A_zeta
    FOR z=0, nz-1 DO BEGIN
      Apsi[*,*,z] = Apsi[*,*,z] * (g.Rxy * g.Bpxy)^2
      Athe[*,*,z] = Athe[*,*,z] / (g.hthe^2)
      Azeta[*,*,z] = Azeta[*,*,z] / (g.Rxy^2)
    ENDFOR
    

    ; Display A just inside vessel (should be parallel to J)
    
    Ad_z = FLTARR(nz, ny)
    Ad_t = Ad_z

    FOR z=0, nz-1 DO BEGIN
      Ad_z[z,*] = A_zeta[nx-4,*,z]
      Ad_t[z,*] = A_the[nx-4,*,z]
    ENDFOR
    
    !P.color = 1
    PLOT_FIELD, Ad_z, Ad_t, title="A", n=500
    cursor, x, y, /down

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Have A in toroidal coords (A_psi, A_the, A_zeta)
    ; Calculate B field from B = Curl(A)

    PRINT, "Calculating B field"

    ; psi derivatives
    dAtdp = Ar
    dAzdp = Ar

    FOR y=0, ny-1 DO BEGIN
      FOR z=0, nz-1 DO BEGIN
        dAtdp[*,y,z] = DERIV(g.psixy[*,y], A_the[*,y,z])
        dAzdp[*,y,z] = DERIV(g.psixy[*,y], A_zeta[*,y,z])
      ENDFOR
    ENDFOR
    
    ; theta derivatives

    dApdt = Ar
    dAzdt = Ar

    FOR x=0, nx-1 DO BEGIN
      FOR z=0, nz-1 DO BEGIN
        dApdt[x,*,z] = fft_deriv(A_psi[x,*,z]) / g.dy[x,0]
        dAzdt[x,*,z] = fft_deriv(A_zeta[x,*,z]) / g.dy[x,0]
      ENDFOR
    ENDFOR
    
    ; zeta derivatives

    dApdz = Ar
    dAtdz = Ar

    FOR x=0, nx-1 DO BEGIN
      FOR y=0, ny-1 DO BEGIN
        dApdz[x,y,*] = fft_deriv(A_psi[x,y,*]) / dzeta
        dAtdz[x,y,*] = fft_deriv(A_the[x,y,*]) / dzeta
      ENDFOR
    ENDFOR
    
    Bpsi = (dAzdt - dAtdz)
    Bthe = (dApdz - dAzdp)
    Bzeta = (dAtdp - dApdt)

    B_psi = Bpsi
    B_the = B_psi
    B_zeta = B_psi

    ; multiply by 1 / sqrt(g)

    Bmag = Bpsi

    FOR x=0, nx-1 DO BEGIN
      FOR y=0, ny-1 DO BEGIN
        Bpsi[x,y,*] = Bpsi[x,y,*] * (g.Bpxy[x,y] / g.hthe[x,y])
        Bthe[x,y,*] = Bthe[x,y,*] * (g.Bpxy[x,y] / g.hthe[x,y])
        Bzeta[x,y,*] = Bzeta[x,y,*] * (g.Bpxy[x,y] / g.hthe[x,y])

        B_psi[x,y,*] = Bpsi[x,y,*] / (g.Rxy[x,y] * g.Bpxy[x,y])^2
        B_the[x,y,*] = Bthe[x,y,*] * g.hthe[x,y]^2
        B_zeta[x,y,*] = Bzeta[x,y,*] * g.Rxy[x,y]^2

        Bmag[x,y,*] = SQRT(Bpsi[x,y,*]*B_psi[x,y,*] + $
                           Bthe[x,y,*]*B_the[x,y,*] + $
                           Bzeta[x,y,*]*B_zeta[x,y,*])
      ENDFOR
    ENDFOR

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Get magnitude of B components

    Bprms = zrms(Bpsi)
    Btrms = zrms(Bthe)
    Bzrms = zrms(Bzeta)
    
    Brms = zrms(Bmag)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Calculate Div(B) as a check
    ; Div(B) = ( d/dx(J*Bx) + d/dy(J*By) + d/dz(J*Bz) )/J

    tmp_p = Bpsi
    tmp_t = Bthe
    tmp_z = Bzeta

    FOR z=0, nz-1 DO BEGIN
      tmp_p[*,*,z] = tmp_p[*,*,z] * g.hthe / g.Bpxy
      tmp_t[*,*,z] = tmp_t[*,*,z] * g.hthe / g.Bpxy
      tmp_z[*,*,z] = tmp_z[*,*,z] * g.hthe / g.Bpxy
    ENDFOR

    ; psi derivative
    FOR y=0, ny-1 DO BEGIN
      FOR z=0, nz-1 DO BEGIN
        tmp_p[*,y,z] = DERIV(g.psixy[*,y], tmp_p[*,y,z])
      ENDFOR
    ENDFOR
    
    ; theta derivative
    FOR x=0, nx-1 DO BEGIN
      FOR z=0, nz-1 DO BEGIN
        tmp_t[x,*,z] = fft_deriv(tmp_t[x,*,z]) / g.dy[x,0]
      ENDFOR
    ENDFOR

    ; zeta derivative
    FOR x=0, nx-1 DO BEGIN
      FOR y=0, ny-1 DO BEGIN
        tmp_z[x,y,*] = fft_deriv(tmp_z[x,y,*]) / dzeta
      ENDFOR
    ENDFOR
    
    divB = (tmp_p + tmp_t + tmp_z)
    FOR z=0,nz-1 DO divB[*,*,z] = divB[*,*,z] * g.Bpxy / g.hthe

    dBrms = zrms(divB)

    PRINT, "Maximum RMS(B) = ", MAX(Brms)
    PRINT, "Maximum RMS(divB) = ", MAX(dBrms)

    plot, g.pressure[*,0], chars=2, title="Pressure", xtitle="X index"

    xnorm = get_integer("Enter X index for normalisation: ")

    bnorm = MAX(Brms[xnorm, *])
    PRINT, "Bnorm = ", bnorm
    
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Calculate J from B as a check

    PRINT, "Calculating current in domain (Should be zero)"

    ; psi derivatives
    dBtdp = Ar
    dBzdp = Ar

    FOR y=0, ny-1 DO BEGIN
      FOR z=0, nz-1 DO BEGIN
        dBtdp[*,y,z] = DERIV(g.psixy[*,y], B_the[*,y,z])
        dBzdp[*,y,z] = DERIV(g.psixy[*,y], B_zeta[*,y,z])
      ENDFOR
    ENDFOR
    
    ; theta derivatives

    dBpdt = Ar
    dBzdt = Ar

    FOR x=0, nx-1 DO BEGIN
      FOR z=0, nz-1 DO BEGIN
        dBpdt[x,*,z] = fft_deriv(B_psi[x,*,z]) / g.dy[x,0]
        dBzdt[x,*,z] = fft_deriv(B_zeta[x,*,z]) / g.dy[x,0]
      ENDFOR
    ENDFOR
    
    ; zeta derivatives

    dBpdz = Ar
    dBtdz = Ar

    FOR x=0, nx-1 DO BEGIN
      FOR y=0, ny-1 DO BEGIN
        dBpdz[x,y,*] = fft_deriv(B_psi[x,y,*]) / dzeta
        dBtdz[x,y,*] = fft_deriv(B_the[x,y,*]) / dzeta
      ENDFOR
    ENDFOR
    
    MU0 = 4.e-7 * !PI

    Jpsi = (dBzdt - dBtdz) / MU0
    Jthe = (dBpdz - dBzdp) / MU0
    Jzeta = (dBtdp - dBpdt) / MU0

    J_psi = Jpsi
    J_the = J_psi
    J_zeta = J_psi

    ; multiply by 1 / sqrt(g)

    Jmag = Jpsi

    FOR x=0, nx-1 DO BEGIN
      FOR y=0, ny-1 DO BEGIN
        Jpsi[x,y,*] = Jpsi[x,y,*] * (g.Bpxy[x,y] / g.hthe[x,y])
        Jthe[x,y,*] = Jthe[x,y,*] * (g.Bpxy[x,y] / g.hthe[x,y])
        Jzeta[x,y,*] = Jzeta[x,y,*] * (g.Bpxy[x,y] / g.hthe[x,y])

        J_psi[x,y,*] = Jpsi[x,y,*] / (g.Rxy[x,y] * g.Bpxy[x,y])^2
        J_the[x,y,*] = Jthe[x,y,*] * g.hthe[x,y]^2
        J_zeta[x,y,*] = Jzeta[x,y,*] * g.Rxy[x,y]^2

        Jmag[x,y,*] = SQRT(Jpsi[x,y,*]*J_psi[x,y,*] + $
                           Jthe[x,y,*]*J_the[x,y,*] + $
                           Jzeta[x,y,*]*J_zeta[x,y,*])
      ENDFOR
    ENDFOR

    Jrms = zrms(Jmag)

    PRINT, "Maximum RMS(J) = ", MAX(Jrms)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Cleaning process
    ; could try to refine A so that J = 0
    
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Put into field-aligned coordinates

    PRINT, "Outputing perturbation of 1T at x = ", xnorm
    
    Ax = Apsi / bnorm
    Ay = Athe / bnorm
        
    pitch=g.hthe*g.Btxy/(g.Rxy*g.Bpxy)
    Az = Ax
    FOR z=0, nz-1 DO Az[*,*,z] = Azeta - g.sinty*Ax - pitch*Ay
    Az = Az / bnorm

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Get A parallel

    
    Apar = Ax
    FOR z=0, nz-1 DO Apar[*,*,z] = (Ax[*,*,z] * g.Btxy * g.sinty * g.Rxy + $
                                    Ay[*,*,z] * (g.Bxy^2) * g.hthe / g.Bpxy + $
                                    Az[*,*,z] * g.Btxy * g.Rxy) / g.Bxy

    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Calculate and plot resonance
    
    mode_structure, apar, g, period=1, famp=famp, /quiet
    
    qprof = ABS(g.shiftangle) / (2.0*!PI) ; q profile
    
    mind = MIN(famp)
    maxd = MAX(famp)
    
    nlev = 50
    lev=mind + (maxd-mind)*indgen(nLev)/(nLev-1)
    col=2+253*indgen(nLev)/(nLev-1)

  
    ; Define red-blue color table
    common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
    
    red = BYTARR(256)
    green = red
    blue = red
    
    ; need to keep color[0] = white, color[1] = black
    red[0] = 255
    green[0] = 255
    blue[0] = 255
    
    ; now create scale
    
    FOR i=2, 255 DO BEGIN
      green[i] = 256 - 2*ABS(i - 128.5)
      
      IF i GT 129 THEN blue[i] = 256 - 2*ABS(i - 128.5) ELSE blue[i] = 255
      IF i LE 129 THEN red[i] = 256 - 2*ABS(i - 128.5) ELSE red[i] = 255
      
    ENDFOR
    
    tvlct,red,green,blue

    qm = n_mode * qprof ; poloidal mode-number

    max_m = FIX(1.2*max(qm))

    CONTOUR, famp[*,0:(max_m-1)], g.psixy[*,0], indgen(max_m)+1, nlev=50, /fill, lev=lev, c_col=col
    
    oplot, g.psixy[*,0], qm

    STOP

  ENDREP UNTIL 0
  
END

