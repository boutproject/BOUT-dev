; Calculates mode structure from BOUT++ output
; for comparison to ELITE
;
; April 2009 - Added ERGOS flag. This is intended
;              for producing plots similar to the ERGOS
;              vacuum RMP code

; interpolates a 1D periodic function
FUNCTION zinterp, v, zind
  v = REFORM(v)
  nz = N_ELEMENTS(v)
  z0 = ROUND(zind)

  p = zind - FLOAT(z0)          ; between -0.5 and 0.5
  
  IF p LT 0.0 THEN BEGIN
      z0 = z0 - 1
      p = p + 1.0
  ENDIF
  
  z0 = ((z0 MOD (nz-1)) + (nz-1)) MOD (nz-1)
  
  ; for now 3-point interpolation
          
  zp = (z0 + 1) MOD (nz - 1)
  zm = (z0 - 1 + (nz-1)) MOD (nz - 1)
  
  result = 0.5*p*(p-1.0)*v[zm] $
    + (1.0 - p*p)*v[z0] $
    + 0.5*p*(p+1.0)*v[zp]

  RETURN, result
END

PRO mode_structure, var_in, grid_in, period=period, $
                    zangle=zangle, n=n, addq=addq, output=output, $
                    xq=xq, xpsi=xpsi, slow=slow, subset=subset, $
                    filter=filter, famp=famp, quiet=quiet, $
                    ergos=ergos, title=title, $
                    xrange=xrange, yrange=yrange, rational=rational, pmodes=pmodes, $
                    _extra=_extra


  ON_ERROR, 2
  
  IF NOT KEYWORD_SET(period) THEN period = 1 ; default = full torus
  IF NOT KEYWORD_SET(zangle) THEN zangle = 0.0
  IF NOT KEYWORD_SET(n) THEN BEGIN
    IF KEYWORD_SET(filter) THEN n = filter*period ELSE n = period
  ENDIF

  IF (grid_in.JYSEPS1_1 GE 0) OR (grid_in.JYSEPS1_2 NE grid_in.JYSEPS2_1) OR (grid_in.JYSEPS2_2 NE grid_in.ny-1) THEN BEGIN
    PRINT, "Mesh contains branch-cuts. Keeping only core"
    
    grid = core_mesh(grid_in)
    var = core_mesh(var_in, grid_in)
  ENDIF ELSE BEGIN
    grid = grid_in
    var = var_in
  ENDELSE

  IF KEYWORD_SET(filter) THEN BEGIN
    var = zfilter(var, filter)
  ENDIF

  nx = grid.nx
  ny = grid.ny

  s = SIZE(var, /dimensions)
  IF N_ELEMENTS(s) NE 3 THEN BEGIN
      PRINT, "Error: Variable must be 3 dimensional"
      RETURN
  ENDIF
  IF (s[0] NE nx) OR (s[1] NE ny) THEN BEGIN
      PRINT, "Error: Size of variable doesn't match grid"
      
      RETURN
  ENDIF
  nz = s[2]

  dz = 2.0*!PI / FLOAT(period*(nz-1))
  
  ; GET THE TOROIDAL SHIFT
  tn = TAG_NAMES(grid)
  tn = STRUPCASE(tn)
  w = WHERE(tn EQ "QINTY", count)
  IF count GT 0 THEN BEGIN
      PRINT, "Using qinty as toroidal shift angle"
      zShift = grid.qinty
  ENDIF ELSE BEGIN
      w = WHERE(tn EQ "ZSHIFT", count)
      IF count GT 0 THEN BEGIN
          PRINT, "Using zShift as toroidal shift angle"
          zShift = grid.zShift
      ENDIF ELSE BEGIN
          PRINT, "ERROR: Can't find qinty or zShift variable"
          RETURN
      ENDELSE
  ENDELSE

  np = 4*ny

  nf = (np - 2) / 2
  famp = FLTARR(nx, nf)

  FOR x=0, nx-1 DO BEGIN
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; transform data into fixed poloidal angle
      
      ; get number of poloidal points
      nskip = FLTARR(ny-1)
      FOR y=0, ny-2 DO BEGIN
          yp = y + 1
          
          nskip[y] = ABS(zshift[x,yp] - zshift[x,y]) / dz - 1
      ENDFOR
      nskip = LONG(ROUND(nskip)) > 0L
      ny2 = LONG(ny + TOTAL(nskip)) ; number of poloidal points

      IF NOT KEYWORD_SET(quiet) THEN PRINT, x, ny2

      f = FLTARR(ny2) ; array for values
      R = f      ; Rxy
      Z = f      ; Zxy
      BtBp = f   ; Bt / Bp
      
      ; interpolate values onto points
      ypos = 0l
      FOR y=0, ny-2 DO BEGIN
          ; original points
          zind = (zangle - zshift[x,y])/dz
          IF N_ELEMENTS(zind) NE 1 THEN STOP
          f[ypos] = zinterp(var[x,y,*], zind)
          R[ypos] = grid.rxy[x,y]
          Z[ypos] = grid.zxy[x,y]
          BtBp[ypos] = grid.Btxy[x,y] / grid.Bpxy[x,y]

          ypos = ypos + 1

          ; add the extra points
          
          zi0 = (zangle - zshift[x,y])/dz
          zip = (zangle - zshift[x,y+1])/dz

          dzi = (zip - zi0) / (nskip[y] + 1)

          FOR i=0l, nskip[y]-1 DO BEGIN
              zi = zi0 + FLOAT(i+1)*dzi ; zindex 
              w = FLOAT(i+1)/FLOAT(nskip[y]+1) ; weighting
              
              f[ypos+i] = w*zinterp(var[x,y+1,*], zi) + (1.0-w)*zinterp(var[x,y,*], zi)
              
              R[ypos+i] = w*grid.rxy[x,y+1] + (1.0-w)*grid.rxy[x,y]
              Z[ypos+i] = w*grid.zxy[x,y+1] + (1.0-w)*grid.zxy[x,y]
              BtBp[ypos+i] = (w*grid.Btxy[x,y+1] + (1.0-w)*grid.Btxy[x,y]) / (w*grid.Bpxy[x,y+1] + (1.0-w)*grid.Bpxy[x,y])
          ENDFOR
          ypos = ypos + nskip[y]
          
          ; final point

          zind = (zangle - zShift[x,ny-1])/dz
          
          f[ypos] = zinterp(var[x,ny-1,*], zind)
          R[ypos] = grid.rxy[x,ny-1]
          Z[ypos] = grid.zxy[x,ny-1]
          BtBp[ypos] = grid.Btxy[x,ny-1] / grid.Bpxy[x,ny-1]
      ENDFOR

      ;STOP

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; calculate poloidal angle
      
      drxy = DERIV(r)
      dzxy = DERIV(z)
      dl = SQRT(drxy*drxy + dzxy*dzxy)
      
      nu = dl * BtBp / R ; field-line pitch
      theta = REAL_PART(fft_integrate(nu)) / grid.shiftangle[x]
      
      IF MAX(theta) GE 1.0 THEN BEGIN
          ; mis-match between q and nu (integration error?)
          IF NOT KEYWORD_SET(quiet) THEN PRINT, "Mismatch  ", x, MAX(theta)
          theta = theta / (MAX(theta) + ABS(theta[1] - theta[0]))
      ENDIF
      
      theta = 2.0*!PI * theta

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; take Fourier transform in theta angle

      tarr = 2.0*!PI*FINDGEN(np) / FLOAT(np) ; regular array in theta

      farr = INTERPOL(f, theta, tarr)

      ;STOP

      ff = FFT(farr)

      FOR i=0, nf-1 DO BEGIN
          famp[x, i] = 2.0*ABS(ff[i+1])
      ENDFOR
      
  ENDFOR
  
  ; sort modes by maximum size

  fmax = FLTARR(nf)
  FOR i=0, nf-1 DO BEGIN
      fmax[i] = MAX(famp[*,i])
  ENDFOR

  inds = REVERSE(SORT(fmax))

  IF NOT KEYWORD_SET(pmodes) THEN pmodes = 10

  qprof = ABS(grid.shiftangle) / (2.0*!PI)

  xarr = FINDGEN(nx)
  xtitle="Radial index"
  IF KEYWORD_SET(xq) THEN BEGIN
      ; show as a function of q*n
      xarr = qprof*FLOAT(n)

      xtitle="q * n"
  ENDIF ELSE IF KEYWORD_SET(xpsi) THEN BEGIN
      ; show as a function of psi. Should be normalised psi
      xarr = REFORM(grid.psixy[*,0])

      ; Check if the grid includes psi axis and boundary
      w = WHERE(tn EQ "PSI_AXIS", count1)
      w = WHERE(tn EQ "PSI_BNDRY", count2)
      
      IF (count1 GT 0) AND (count2 GT 0) THEN BEGIN
        xarr = (xarr - grid.psi_axis) / (grid.psi_bndry - grid.psi_axis)
      ENDIF ELSE BEGIN
        ; Use hard-wired values
        PRINT, "WARNING: Using hard-wired psi normalisation"
        ; for circular case
        ;xarr = (xarr + 0.1937) / (0.25044 + 0.1937)
        ; for ellipse case
        ;xarr = xarr / 0.74156
        
        ; cbm18_dens8
        xarr = (xarr + 0.854856) / (0.854856 + 0.0760856)
      ENDELSE
      
      xtitle="Psi normalised"
  ENDIF

  
  IF KEYWORD_SET(slow) THEN BEGIN
    ; plot modes slowly for examination
    safe_colors, /first
    
    ; go through and plot each mode
    FOR i=0, nf-1 DO BEGIN
      IF max(famp[*,i]) GT 0.05*max(famp) THEN BEGIN
        PRINT, "Mode m = ", i+1, " of ", nf
        plot, xarr, famp[*,i], yr=[0,MAX(famp)], color=1, xtitle=xtitle, chars=1.5, xrange=xrange, _extra=_extra
        q = FLOAT(i+1) / FLOAT(n)
        
        pos = INTERPOL(xarr, qprof, q)
        
        oplot, [pos, pos], [0, 2.*MAX(fmax)], lines=2, color=1
        cursor, a, b, /down
      ENDIF
    ENDFOR
    
  ENDIF ELSE IF KEYWORD_SET(ergos) THEN BEGIN
    ; ERGOS - style output
    
    IF KEYWORD_SET(output) AND NOT KEYWORD_SET(slow) THEN BEGIN
      SET_PLOT, 'PS'
      DEVICE, file=output+".ps", /color
    ENDIF

    contour2, famp, xarr, indgen(nf)+1, $
              xtitle=xtitle, xrange=xrange, yrange=yrange, _extra=_extra

    ; overplot the q profile

    oplot, xarr, qprof * n, color=1, thick=2
  
    IF KEYWORD_SET(rational) THEN BEGIN
      maxm = FIX(MAX(qprof)) * n

      qreson = (FINDGEN(maxm)+1) / FLOAT(n)

      ; get x location for each of these resonances
      qloc = INTERPOL(xarr, qprof, qreson)

      oplot, qloc, findgen(maxm)+1., psym=4, color=1
    ENDIF
    
    IF KEYWORD_SET(output) THEN BEGIN
      ; output data to save file
      SAVE, xarr, qprof, famp, file=output+".idl"
      
      DEVICE, /close
      SET_PLOT, 'X'
    ENDIF

  ENDIF ELSE BEGIN
    IF KEYWORD_SET(output) AND NOT KEYWORD_SET(slow) THEN BEGIN
      SET_PLOT, 'PS'
      DEVICE, file=output+".ps", /color
    ENDIF
    
    safe_colors, /first
    
    IF KEYWORD_SET(subset) THEN BEGIN
      
      ; get number of modes larger than 5% of the maximum
      w = WHERE(fmax GT 0.10*MAX(fmax), count)
      
      minind = MIN(inds[0:(count-1)])
      maxind = MAX(inds[0:(count-1)])
      
      PRINT, "Mode number range: ", minind, maxind
      
      plot, xarr, famp[*,0], yr=[0,MAX(famp)], color=1, xtitle=xtitle, chars=1.5, xrange=xrange, /nodata,title=title, _extra=_extra
      color = 2
      FOR i=minind, maxind, subset DO BEGIN
        oplot, xarr, famp[*,i], color=color
        
        q = FLOAT(i+1) / FLOAT(n)
        pos = INTERPOL(xarr, qprof, q)
        
        oplot, [pos, pos], [0, 2.*MAX(fmax)], lines=2, color=color
        
        color = (color MOD 2) + 1
      ENDFOR
      
    ENDIF ELSE BEGIN
      ; default - just plot everything
      plot, xarr, famp[*,0], yr=[0,MAX(famp)], color=1, xtitle=xtitle, chars=1.5, xrange=xrange,title=title, _extra=_extra
      
      color = 2
      FOR i=0, nf-1 DO BEGIN
        oplot, xarr, famp[*,i], color=color
        
        color = (color MOD 2) + 1
      ENDFOR
          
      IF KEYWORD_SET(addq) THEN BEGIN
        
        FOR i=0, pmodes-1 DO BEGIN
          PRINT, "m = "+STRTRIM(STRING(inds[i]+1), 2)+" amp = "+STRTRIM(STRING(fmax[inds[i]]),2)
          q = FLOAT(inds[i]+1) / FLOAT(n)
          
          pos = INTERPOL(xarr, qprof, q)
          
          oplot, [pos, pos], [0, 2.*MAX(fmax)], lines=2, color=1
        ENDFOR
      ENDIF
      
    ENDELSE
    IF KEYWORD_SET(output) THEN BEGIN
      ; output data to save file
      SAVE, xarr, qprof, famp, file=output+".idl"
      
      DEVICE, /close
      SET_PLOT, 'X'
    ENDIF
  ENDELSE
  
  ;STOP
END
