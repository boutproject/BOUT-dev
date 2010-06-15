;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Solver for vacuum region, given B at plasma boundary.   ;
; Works in cylindrical geometry (i.e. tokamaks)           ;
;                                                         ;
; NOTE: Assumes closed flux-surfaces. NO X-POINTS         ;
;                                                         ;
; Inputs:                                                 ;
;   Rs     - 1D array of major radius for surface points  ;
;   Zs     - 1D array of height (Z)                       ;
;   Bps    - Poloidal field at each point                 ;
;   Bts    - Toroidal field                               ;
;   dpsi   - Change in psi between output flux-surfaces   ;
;                                                         ;
; Optional inputs:                                        ;
;   nrad   - number of surfaces to add (default 5)        ;
;   sm     - Poloidal smoothing length. Not really needed ;
;   output - Name of a PS file to plot to                 ;
;                                                         ;
; Output: Structure containing                            ;
;   r      - 2D array of R values (including surface)     ;
;   z      - 2D array of Z value                          ;
;   bp     - Poloidal field                               ;
;   bt     - Toroidal field                               ;
;  These are all of the form r[nrad + 1, npoloidal]       ;
;  so that r[0,*] = rs etc.                               ;
;                                                         ;
; B.Dudson, University of York, April 2008                ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION int_func, x, f
   n = N_ELEMENTS(f)

   n2 = FIX(n/2)

   g = FLTARR(n)

   g[0] = 0.0
   FOR i=n2, n-1 DO BEGIN
       g[i] = INT_TABULATED(x[0:i], f[0:i])
   ENDFOR

   FOR i=1, n2-1 DO BEGIN
       g[i] = g[n-1] - INT_TABULATED(x[i:*], f[i:*])
   ENDFOR

   RETURN, g
END

FUNCTION loop_smooth, f, n
  IF n LT 3 THEN n = 3
  IF n MOD 2 EQ 0 THEN n = n + 1 ; n must be odd

  len = N_ELEMENTS(f)
  nedge = (n-1)/2

  IF nedge GT len THEN BEGIN
      PRINT, "Error: Array to be smoothed is too small"
      RETURN, f
  ENDIF

  arr = [REFORM(f[(len-nedge):*]),REFORM(f), REFORM(f[0:nedge-1])]

  arr = SMOOTH(arr, n)
  arr = arr[nedge:nedge+len-1]
  
  RETURN, arr
END

FUNCTION vacuum, Rs, Zs, Bps, Bts, dpsi, nrad=nrad, $
                 rev=rev, sm=sm, $
                 output=output, max_iter=max_iter

  IF NOT KEYWORD_SET(max_iter) THEN max_iter = 20  

  safe_colors, /first

  npol = N_ELEMENTS(Rs)

  rbt = MEAN(Rs*Bts)

start:

  IF KEYWORD_SET(rev) THEN BEGIN
      r = REVERSE(rs)
      z = REVERSE(zs)
      Bp = REVERSE(bps)
      Bt = REVERSE(bts)
      rev = 1
  ENDIF ELSE BEGIN
      r = rs
      z = zs
      bp = bps
      bt = bts
      rev = 0
  ENDELSE
  
  ; check input
  IF (N_ELEMENTS(Z) NE npol) OR (SIZE(R, /n_dim) NE 1) OR (SIZE(Z, /n_dim) NE 1) THEN BEGIN
      PRINT, "Error: R and Z must be 1D vectors of the same length"
      RETURN, 0
  ENDIF
  
  ; get settings
  IF NOT KEYWORD_SET(nrad) THEN nrad = 5

  psiarr = FINDGEN(nrad+1) * dpsi 

  PRINT, "Creating starting mesh"
  
  Rxy = FLTARR(nrad+1, npol)
  Zxy = Rxy
  Bpxy = Rxy
  hthe = Rxy
  psi = Rxy

  ; at each point take a vector orthogonal to initial surface
  ; doesn't have to be particularly accurate (refined later)
  ; works as long as shape is convex

  Rxy[0,*] = R
  Zxy[0,*] = Z
  Bpxy[0,*] = Bp
  
  ; get maximum distance necessary for this psi
  ; dpsi = R*Bp*dr (NO 2*pi!)
      
  distance = dpsi*FLOAT(nrad) / (MIN(R)*MIN(Bp))

  plot, R, Z, psym=1, /iso, color=1, $
    xr=[MIN(R)-distance*2, MAX(R)+distance*2], $
    yr=[MIN(Z)-distance*2, MAX(Z)+distance*2]

  ; get the local flux-surface pitch
  dr = DERIV([Rxy[0, npol-1], REFORM(Rxy[0,*]), Rxy[0,0]])
  dz = DERIV([Zxy[0, npol-1], REFORM(Zxy[0,*]), Zxy[0,0]])
  dr = dr[1:npol]
  dz = dz[1:npol]

  ;dr = fft_deriv(Rxy[0,*])
  ;dz = fft_deriv(Zxy[0,*])

  dl0 = SQRT(dr^2 + dz^2)

  hthe[0,*] = dl0

  ang = ATAN(dz, dr)        ; local pitch of the surface

  dl = dl0
  FOR i=1, nrad DO BEGIN

      IF KEYWORD_SET(dpsi) THEN distance = dpsi*FLOAT(nrad)*dl/(Rxy[i-1,*]*Bp*dl0)
      
      Rxy[i,*] = Rxy[i-1,*] - (distance/FLOAT(nrad)) * sin(!PI - ang)
      Zxy[i,*] = Zxy[i-1,*] - (distance/FLOAT(nrad)) * cos(!PI - ang)

      dr = DERIV([Rxy[i, npol-1], REFORM(Rxy[i,*]), Rxy[i,0]])
      dz = DERIV([Zxy[i, npol-1], REFORM(Zxy[i,*]), Zxy[i,0]])
      
      dr = dr[1:npol]
      dz = dz[1:npol]
      
      ;dr = fft_deriv(Rxy[i-1,*])
      ;dz = fft_deriv(Zxy[i-1,*])

      dl = SQRT(dr^2 + dz^2)

      IF i EQ 1 THEN BEGIN
          ; Check that the direction is correct
          ; inner surface should be smaller than outer

          len0 = int_tabulated(findgen(npol), dl0)
          len1 = int_tabulated(findgen(npol), dl)

          IF len1 LT len0 THEN BEGIN
              ; going inwards not outwards
              
              PRINT, "Reversing direction"
              IF rev THEN rev = 0 ELSE rev = 1
              GOTO, start
          ENDIF
      ENDIF

      hthe[i,*] = dl

      ; set initial guess for Bp
      Bpxy[i,*] = Bp * (dl0 / dl)

      oplot, Rxy[i,*], Zxy[i,*], psym=3, color=1

  ENDFOR


  rold = rxy
  zold = zxy
  bold = Bpxy
  ;;;;;;;;;; MAIN LOOP ;;;;;;;;;;;;

  iter = 0

  REPEAT BEGIN
      iter = iter + 1

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; calculate psi at each point ;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      f = Rxy * Bpxy

      ; integrate f along radial coordinate to get psi
      
      FOR i=0, npol-1 DO BEGIN
          dr = DERIV(REFORM(rxy[*,i]))
          dz = DERIV(REFORM(zxy[*,i]))
          
          dl = SQRT(dr^2 + dz^2)
          
          l = int_func(FINDGEN(nrad+1), dl)
          
          psi[*,i] = int_func(l, REFORM(f[*,i]))
      ENDFOR

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Create new poloidal mesh     ;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      r2 = Rxy
      z2 = Zxy
      bp2 = Bpxy
      FOR i=0, npol-1 DO BEGIN
          ; get indices for psi iso-surfaces
          inds = INTERPOL(FINDGEN(nrad+1), REFORM(psi[*, i]), psiarr)

          r2[*,i] = INTERPOL(REFORM(Rxy[*,i]), FINDGEN(nrad+1), inds)
          z2[*,i] = INTERPOL(REFORM(Zxy[*,i]), FINDGEN(nrad+1), inds)
          bp2[*,i] = INTERPOL(REFORM(Bpxy[*,i]), FINDGEN(nrad+1), inds)
      ENDFOR

      IF KEYWORD_SET(sm) THEN BEGIN
          ; smooth the grid
          
          FOR i=1, nrad DO BEGIN
              r2[i,*] = loop_smooth(r2[i,*], sm)
              z2[i,*] = loop_smooth(z2[i,*], sm)
              ;bp2[i,*] = loop_smooth(bp2[i,*], sm)
          ENDFOR
      ENDIF

      FOR i=1, nrad DO BEGIN
          OPLOT, r2[i,*], z2[i,*], color=2
      ENDFOR

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Create new orthogonal radial mesh ;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      orthog = 1
      
      IF orthog EQ 1 THEN BEGIN
          theta_inner = 2.0*!PI*DINDGEN(npol)/DOUBLE(npol)
          transform = gen_orthog(r2, z2, theta_inner, tol=1.0e-2)

          ; turn into an index
          transform = transform * FLOAT(npol) / (2.0*!PI)

          ; interpolate to new grid-points

          FOR i=1, nrad DO BEGIN
              fr = FFT(REFORM(r2[i,*]))
              fz = FFT(REFORM(z2[i,*]))
              fb = FFT(REFORM(Bp2[i,*]))
              
              FOR j=0, npol-1 DO BEGIN
                  Rxy[i,j] = fft_interp(fr, transform[i,j])
                  Zxy[i,j] = fft_interp(fz, transform[i,j])
                  Bpxy[i,j] = fft_interp(fb, transform[i,j])
              ENDFOR
          ENDFOR
      ENDIF ELSE BEGIN
          Rxy = r2
          Zxy = z2
          Bpxy = Bp2
      ENDELSE

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Re-calculate hthe and Bpol  ;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      FOR i=0, nrad DO BEGIN
          dr = DERIV([Rxy[i, npol-1], REFORM(Rxy[i,*]), Rxy[i,0]])
          dz = DERIV([Zxy[i, npol-1], REFORM(Zxy[i,*]), Zxy[i,0]])
      
          dr = dr[1:npol]
          dz = dz[1:npol]
      
          ;dr = fft_deriv(Rxy[i-1,*])
          ;dz = fft_deriv(Zxy[i-1,*])

          hthe[i,*] = SQRT(dr^2 + dz^2)

          
          IF i NE 0 THEN Bpxy[i,*] = Bp * REFORM(hthe[0,*] / hthe[i,*])
          
          IF KEYWORD_SET(sm) THEN Bpxy[i,*] = loop_smooth(Bpxy[i,*], sm)
          
      ENDFOR

      ;;;;;;;;;;;;;;;;;
      ; Plot new grid ;
      ;;;;;;;;;;;;;;;;;

      distance = MAX(rxy) - MAX(rxy[0,*])
      
      plot, R, Z, psym=1, /iso, color=1, $
        xr=[MIN(R)-distance*2, MAX(R)+distance*2], $
        yr=[MIN(Z)-distance*2, MAX(Z)+distance*2]
      FOR i=1, nrad DO BEGIN
          oplot, rxy[i,*], zxy[i,*], color=1
      ENDFOR
      
      FOR i=0, npol-1 DO BEGIN
          oplot, rxy[*,i], zxy[*,i], color=2
      ENDFOR

      dl = SQRT((rxy - rold)^2 + (zxy - zold)^2)
      db = MAX(ABS(bpxy - bold))

      rold = rxy
      zold = zxy
      bold = bpxy

      PRINT, "Change in location:", iter, MEAN(dl), MAX(dl), db

      ;cursor, a, b, /down

      ;STOP

  ENDREP UNTIL (db LT 1.0e-4) OR (iter EQ max_iter) ;(MEAN(dl) LT 1.0e-3) OR (iter EQ 5)

  IF KEYWORD_SET(output) THEN BEGIN
      SET_PLOT, 'PS'
      DEVICE, file=output, /color

      safe_colors, /first

      plot, R, Z, psym=1, /iso, color=1, $
        xr=[MIN(R)-distance*2, MAX(R)+distance*2], $
        yr=[MIN(Z)-distance*2, MAX(Z)+distance*2]
      FOR i=1, nrad DO BEGIN
          oplot, rxy[i,*], zxy[i,*], color=1
      ENDFOR
      
      FOR i=0, npol-1 DO BEGIN
          oplot, rxy[*,i], zxy[*,i], color=2
      ENDFOR
      
      DEVICE, /close
      SET_PLOT, 'X'
  ENDIF

  IF rev THEN BEGIN
      ; need to reverse the poloidal direction
      FOR i=0, nrad DO BEGIN
          rxy[i,*] = REVERSE(REFORM(rxy[i,*]))
          zxy[i,*] = REVERSE(REFORM(zxy[i,*]))
          bpxy[i,*] = REVERSE(REFORM(bpxy[i,*]))
      ENDFOR
  ENDIF
  
  ; calculate toroidal field
  btxy = rbt / rxy

  result = {r:rxy, z:zxy, bp:bpxy, bt:btxy}

  RETURN, result
END
