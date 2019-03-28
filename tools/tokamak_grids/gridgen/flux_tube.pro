; Generate a grid file for a tokamak flux-tube simulation
;
; For example, a DIII-D case. 
; Not a particularly good case, as pressure goes negative
;
; IDL> flux_tube, "efit/neqdsk", 0.8D
; 
; NOTE: NOT WORKING PROPERLY YET. MAGNETIC SHEAR MUCH TOO LARGE
;       
;                 DO NOT USE FOR RUNS!!
;
;

PRO flux_tube, gfile, psinorm, output=output
  
  IF NOT KEYWORD_SET(output) THEN output="fluxtube"+STR(psinorm)+".grd.nc"

  IF (psinorm LT 0.D) OR (psinorm GE 1.D) THEN BEGIN
     PRINT, "ERROR: input psinorm must be between 0 and 1"
     RETURN
  ENDIF
  
  ; Read the EFIT G-EQDSK file
  g = read_neqdsk(gfile)
  
  ; Check the return
  IF SIZE(g, /TYPE) NE 8 THEN BEGIN
     PRINT, "ERROR: Couldn't read G-EQDSK file '"+gfile+"'"
     RETURN
  ENDIF
  
  ; Extract needed data from g-file struct
  
  rzgrid = {nr:g.nx, nz:g.ny, $ ; Number of grid points
            r:REFORM(g.r[*,0]), z:REFORM(g.z[0,*]), $ ; R and Z as 1D arrays
            simagx:g.simagx, sibdry:g.sibdry, $       ; Range of psi
            psi:g.psi, $       ; Poloidal flux in Weber/rad on grid points
            pres:g.pres, $     ; Plasma pressure in nt/m^2 on uniform flux grid
            qpsi:g.qpsi, $     ; q values on uniform flux grid
            nlim:g.nlim, rlim:g.xlim, zlim:g.ylim} ; Wall boundary
  
  ; Plot psi
  nlev = 100
  minf = MIN(rzgrid.psi)
  maxf = MAX(rzgrid.psi)
  levels = findgen(nlev)*(maxf-minf)/DOUBLE(nlev-1) + minf
  
  safe_colors, /first
  CONTOUR, rzgrid.psi, rzgrid.R, rzgrid.Z, $
           levels=levels, color=1, /iso, xstyl=1, ysty=1
  
  ; Find O- and X- points
  critical = analyse_equil(rzgrid.psi, rzgrid.R, rzgrid.Z)
  
  ; Overplot the separatrices, O-points
  oplot_critical, rzgrid.psi, rzgrid.R, rzgrid.Z, critical
  
  ; Use primary O-point and X-point to convert normalised psi into psi
  
  psi_axis = critical.opt_f[critical.primary_opt]
  psi_sep = critical.xpt_f[critical.inner_sep]
  
  ; See how this compares to grid values
  PRINT, "Psi Axis : ", psi_axis, g.simagx
  PRINT, "Psi Bndry: ", psi_sep, g.sibdry
  
  psi_axis = g.simagx
  psi_sep = g.sibdry

  psi = psi_axis + psinorm*(psi_sep - psi_axis)
  
  ;;;;;;;;;;;;;;; Calculate DCT ;;;;;;;;;;;;;;
  
  PRINT, "Calculating DCT..."
  DCT2Dslow, rzgrid.psi, dctpsi
  PRINT, "Finished DCT"
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Use contour to get a flux surface
  
  contour_lines, rzgrid.psi, findgen(rzgrid.nr), findgen(rzgrid.nz), $
                 levels=[psi], $
                 path_info=info, path_xy=xy
  
  IF N_ELEMENTS(info) GT 1 THEN BEGIN
     ; Find the surface closest to the o-point
      
     ind = closest_line(info, xy, critical.opt_ri[critical.primary_opt], critical.opt_zi[critical.primary_opt])
     info = info[ind]
  ENDIF ELSE info = info[0]
  
  ; Extract r and z indices
  ri = REFORM(xy[0,info.offset:(info.offset+info.n-1)])
  zi = REFORM(xy[1,info.offset:(info.offset+info.n-1)])
  np = N_ELEMENTS(ri)

  ; Re-order so starts on the inboard midplane
  ; so that twist-shift matching will be on inboard side
  mri = MIN(ri, rmini)
  IF (rmini NE 0) AND (rmini NE (np-1)) THEN BEGIN
    ri = [ri[rmini:*], ri[0:(rmini-1)]]
    zi = [zi[rmini:*], zi[0:(rmini-1)]]
  ENDIF

  ; Refine flux surface location using DCT
  FOR i=0, np-1 DO BEGIN
    follow_gradient, dctpsi, rzgrid.R, rzgrid.Z, ri[i], zi[i], psi, $
                     rinew, zinew, status=status
    IF status THEN BEGIN
      PRINT, "ERROR refining surface location"
      RETURN
    ENDIF
    
    ri[i] = rinew
    zi[i] = zinew
  ENDFOR

  ; Fourier filter to smooth surface
  ;nfreq = 20
  ;ri = REAL_PART(fft_filter(ri, nfreq))
  ;zi = REAL_PART(fft_filter(zi, nfreq))

  ; Plot the flux surface
  OPLOT, INTERPOLATE(rzgrid.R, ri, /DOUBLE), INTERPOLATE(rzgrid.Z, zi, /DOUBLE), color=4, thick=2
  
  ; Get quantities from g-eqdsk
  ngrid = N_ELEMENTS(g.fpol)
  gpos = psinorm * ngrid ; Index into the psi grid. CHECK THIS
  psigrid = psi_axis + (psi_sep - psi_axis)*FINDGEN(ngrid)/DOUBLE(ngrid)
  fpol = INTERPOLATE(g.fpol, gpos, /DOUBLE) ; Poloidal current function
  pres = INTERPOLATE(g.pres, gpos, /DOUBLE) ; Pressure [Pascals]
  dpdpsi = INTERPOLATE(DERIV(psigrid, g.pres), gpos, /DOUBLE)
  dfdpsi = INTERPOLATE(DERIV(psigrid, g.fpol), gpos, /DOUBLE)
  qsafe = INTERPOLATE(g.qpsi, gpos, /DOUBLE) ; q
  
  PRINT, "Pressure [Pa]: " + STR(pres)
  PRINT, "Safety factor: " + STR(qsafe)
  
  Rmaj = INTERPOLATE(rzgrid.R, ri, /DOUBLE) ; Major radius

  ; Use DCT to get local gradients of psi for Bp
  np = N_ELEMENTS(ri)
  Bpol = DBLARR(np)
  FOR i=0, np-1 DO BEGIN
    grad = local_gradient(dctpsi, ri[i], zi[i], status=status)
    IF status THEN BEGIN
      PRINT, "ERROR: Problems evaluating derivatives at "+STR([ri,zi])
      RETURN
    ENDIF
    ; dfd* are derivatives wrt the indices. Need to multiply by dr/di etc
    Bpol[i] = SQRT((grad.dfdr/(DERIV(rzgrid.r))[ri[i]])^2 + $
                   (grad.dfdz/(DERIV(rzgrid.z))[zi[i]])^2) / Rmaj[i] ; NOTE: Magnitude only
    
  ENDFOR
  
  ; Get toroidal field
  Btor = fpol / Rmaj

  B = SQRT(Btor^2 + Bpol^2)
  
  ; Need to work out direction of B and poloidal direction
  
  ; Now integrate along this surface to get toroidal angle and
  ; distance along field-line
  
  ; Poloidal distance dl/di (i = index)
  drdi = REAL_PART(fft_deriv(INTERPOLATE(rzgrid.R, ri, /DOUBLE)))
  dzdi = REAL_PART(fft_deriv(INTERPOLATE(rzgrid.Z, zi, /DOUBLE)))
  dldi = SQRT(drdi^2 + dzdi^2)

  ; Integrate to get distance in poloidal direction
  l = REAL_PART(fft_integrate(dldi, loop=poldist))
  poldist = REAL_PART(poldist)

  ; Change in toroidal angle, following a field-line
  dtdi = dldi * Btor / (Bpol * Rmaj)
  
  ; Integrate this over the poloidal direction to get toroidal angle
  qinty = REAL_PART(fft_integrate(dtdi, loop=qloop))
  qloop = REAL_PART(qloop)
  
  ; Get parallel distance
  dsdi = dldi * B / Bpol
  s = REAL_PART(fft_integrate(dsdi, loop=pardist))
  pardist = REAL_PART(pardist)
  
  PRINT, "Safety factor = "+STR(qloop/(2.D*!DPI))+" ("+STR(qsafe)+" in g-file)"
  PRINT, "Poloidal distance = "+STR(poldist)+"m"
  PRINT, "Parallel distance = "+STR(pardist)+"m"

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Choose parallel grid spacing

  REPEAT BEGIN
    npar = get_integer("Number of parallel grid points (>3): ")
  ENDREP UNTIL npar GT 3
  
  PRINT, "Parallel grid spacing options"
  PRINT, "   1) Equal poloidal distance"
  PRINT, "   2) Equal toroidal angle"
  PRINT, "   3) Equal parallel distance"
  REPEAT BEGIN
    opt = get_integer("Enter option [1-3]: ")
  ENDREP UNTIL (opt GE 1) AND (opt LE 3)
  
  CASE opt OF
    1: BEGIN ; Equal poloidal
      dist = l
      tot = poldist
    END
    2: BEGIN ; Equal toroidal angle
      dist = qinty
      tot = qloop
    END
    3: BEGIN ; Equal parallel
      dist = s
      tot = pardist
    END
  ENDCASE
  
  pos = tot * FINDGEN(npar)/DOUBLE(npar)
  inds = INTERPOL(FINDGEN(np), dist, pos)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Interpolate everything onto grid
  
  ri   = INTERPOLATE(ri, inds, /DOUBLE)
  zi   = INTERPOLATE(zi, inds, /DOUBLE)
  Bpol = INTERPOLATE(Bpol, inds, /DOUBLE)
  Btor = INTERPOLATE(Btor, inds, /DOUBLE)
  B    = SQRT(Bpol^2 + Btor^2)
  Rmaj = INTERPOLATE(rzgrid.R, ri, /DOUBLE)
  Zpos = INTERPOLATE(rzgrid.Z, zi, /DOUBLE)
  qinty = INTERPOLATE(qinty, inds, /DOUBLE)
  s = INTERPOLATE(s, inds, /DOUBLE) ; parallel distance
  l = INTERPOLATE(l, inds, /DOUBLE) ; poloidal distance
  hthe = DERIV(l) / (2.D*!DPI / DOUBLE(npar))
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Calculate psi derivatives
  
  dpsidR  = DBLARR(npar)
  dpsidZ  = DBLARR(npar)
  dRdpsi  = DBLARR(npar) ; dR / dpsi
  dZdpsi  = DBLARR(npar) ; dZ / dpsi
  dBpdR   = DBLARR(npar)
  dBpdZ   = DBLARR(npar)
  FOR i=0, npar-1 DO BEGIN
    ; Get gradients of psi
    grad = EvalCosP(dctpsi, x0=ri[i], y0=zi[i])
    dpsidR[i]    = grad[1]
    dpsidZ[i]    = grad[2]
    d2psidR2  = grad[3]
    d2psidZ2  = grad[4]
    d2psidRdZ = grad[5]

    ; NOTE: THESE MAY BE DERIVATIVES W.R.T INDEX, NOT R/Z
    
    ; dR / dpsi
    dRdpsi[i] = dpsidR[i] / (Rmaj[i] * Bpol[i])^2
    ; dZ / dpsi
    dZdpsi[i] = dpsidZ[i] / (Rmaj[i] * Bpol[i])^2
    
    dBpdR[i] = -Bpol[i]/Rmaj[i] + (dpsidR[i]*d2psidR2 + dpsidZ[i]*d2psidRdZ)/(Bpol[i]*Rmaj[i]^2)
    
    dBpdZ[i] = (dpsidZ[i]*d2psidZ2 + dpsidR[i]*d2psidRdZ)/(Bpol[i]*Rmaj[i]^2)
   
  ENDFOR
  
  ; Psi derivatives of B field components
  dBpdpsi = dBpdR*dRdpsi + dBpdZ*dZdpsi
  dBtdpsi = dfdpsi/Rmaj - (fpol / (Rmaj^2)) * dRdpsi
  dBdpsi = (Btor*dBtdpsi + Bpol*dBpdpsi)/B
  
  ;;;;;;;;;;;;;;;;;; MAGNETIC SHEAR ;;;;;;;;;;;;;;;;;;;;
  ; Calculate from radial force balance equation
  
  MU0 = 4.d-7 * !DPI
  
  pitch = hthe * Btor / (Bpol * Rmaj)
  dnudpsi = - ( (MU0*hthe*dpdpsi/Bpol) + pitch*( 2.D*Rmaj*B*dBdpsi/Bpol + B^2*dRdpsi/Bpol - B^2*Rmaj*dBpdpsi/(Bpol^2) ) ) / (Rmaj*Bpol^2 / Btor)
  
  ; Integrate this to get the integrated shear sinty
  sinty = REAL_PART(fft_integrate(dnudpsi, loop=sloop)) * 2.D*!DPI/DOUBLE(npar) 
  sloop = REAL_PART(sloop) * 2.D*!DPI/DOUBLE(npar)

  ; Want this shift to be zero at the outboard midplane,
  ; and matching location on inboard side
  ; sloop needed for taylor expansion of q at twist-shift
  
  mri = MAX(Rmaj, rmaxi)
  sinty = sinty - sinty[rmaxi]

  ;;;;;;;;;;;;;;;;;;;; CURVATURE ;;;;;;;;;;;;;;;;;;;;;;;
  ; Calculate b x kappa using coordinates of a single field-line
  
  ; First derivatives along the line to get tangent vector
  dr = DERIV(s, Rmaj)  ; R position
  dz = DERIV(s, Zpos)  ; Z position 
  dp = DERIV(s, qinty) ; Toroidal angle
  
  ; Second derivatives
  d2r = DERIV(s, dr)
  d2z = DERIV(s, dz)
  d2p = DERIV(s, dp)
  
  ; Components of b (tangent vector)
  br = dr
  bz = dz
  bp = Rmaj*dp

  ; Components of curvature (unit vectors)
  kr = d2r - Rmaj*dp^2
  kz = d2z
  kp = 2.D*dr*dp + Rmaj*d2p

  ; Calculate bxk in cylindrical coordinates

  bxkr = bp*kz - bz*kp
  bxkz = br*kp - bp*kr
  bxkp = bz*kr - br*kz
  
  ; Calculate components in (psi, theta, phi) toroidal coordinates
  ; (bxk) dot grad psi
  bxkpsi   = bxkr * dpsidR + bxkz * dpsidZ
  ; (bxk) dot grad theta
  bxktheta = ( bxkr * dpsidZ - bxkz * dpsidR ) / (Rmaj*Bpol*hthe)

  ;Finally into field-aligned coordinates
  bxcvx1d = bxkpsi
  bxcvy1d = bxktheta
  bxcvz1d = bxkpsi - sinty*bxkpsi - pitch*bxktheta

  ;;;;;;;;;;;;;; PARALLEL CURRENT ;;;;;;;;;;;;;;;;;

  Jpar = B*dfdpsi/MU0 + Rmaj*MU0*dpdpsi

  PRINT, "Parallel current density [Am^-2]: "+STR(MIN(Jpar))+" -> "+STR(MAX(Jpar))

  ;;;;;;;;;;;;;;;; RADIAL MESH ;;;;;;;;;;;;;;;;;;
  
  REPEAT BEGIN
    nrad = get_integer("Number of radial grid points (>5): ")
  ENDREP UNTIL nrad GT 5
  
  PRINT, "Specify radial size of domain in..."
  PRINT, "   1) cm at the outboard midplane"
  PRINT, "   2) normalised psi"
  PRINT, "   3) psi"
  REPEAT BEGIN
    opt = get_integer("Enter option [1-3]: ")
  ENDREP UNTIL (opt GE 1) AND (opt LE 3)
  
  mri = MAX(Rmaj, rmaxi) ; Get index at outboard midplane
  CASE opt OF
    1: BEGIN
      dr = 0.01D*get_float("Enter radial size [cm]:")
      dpsi = Bpol[rmaxi]*Rmaj[rmaxi]*dr
      dpsin = dpsi / (psi_sep - psi_axis)
    END
    2: BEGIN
      dpsin = get_float("Enter radial size [normalised psi]:")
      dpsi = dpsin  * (psi_sep - psi_axis)
      dr = dpsi / (Bpol[rmaxi]*Rmaj[rmaxi])
    END
    3: BEGIN 
      dpsi = get_float("Enter radial size [psi]:")
      dpsin = dpsi / (psi_sep - psi_axis)
      dr = dpsi / (Bpol[rmaxi]*Rmaj[rmaxi])
    END
  ENDCASE
  
  PRINT, ""
  PRINT, "Radial domain"
  PRINT, "-------------"
  PRINT, "  Width [cm]   "+STR(dr)
  PRINT, "        [psin] "+STR(dpsin)
  PRINT, "        [psi]  "+STR(dpsi)
  PRINT, "  Safety factor variation: +/-"+STR(sloop*dpsi / (2.D*!DPI) / 2.D)
  PRINT, ""
  
  ; Convert into cell width
  dr    = dr / DOUBLE(nrad)
  dpsi  = dpsi / DOUBLE(nrad)
  dpsin = dpsin / DOUBLE(nrad)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Put everything into 2D arrays
  
  rmid = DOUBLE(nrad) / 2.D
  rpsi = (FINDGEN(nrad) - rmid) * dpsi

  ; Pressure profile. Only gradient used as drive term
  pressure = DBLARR(nrad, npar)
  Jpar0 = DBLARR(nrad, npar)
  FOR i=0, nrad-1 DO BEGIN
    pressure[i,*] = pres + rpsi[i]*dpdpsi ; Linear in x
    Jpar0[i,*] = Jpar ; Constant in x
  ENDFOR
  
  ; B field components
  Bpxy = DBLARR(nrad, npar)
  Btxy = DBLARR(nrad, npar)
  Bxy  = DBLARR(nrad, npar)
  FOR i=0, nrad-1 DO BEGIN
    Bpxy[i,*] = Bpol
    Btxy[i,*] = Btor
    Bxy[i,*] = B
  ENDFOR
  
  ; Grid spacing
  dx = DBLARR(nrad, npar) + dpsi
  dy = DBLARR(nrad, npar) + 2.D*!DPI/DOUBLE(npar)

  ; Geometrical quantities
  hxy = DBLARR(nrad, npar)
  Rxy = DBLARR(nrad, npar)
  Zxy = DBLARR(nrad, npar)
  FOR i=0, nrad-1 DO BEGIN
    hxy[i,*] = hthe
    Rxy[i,*] = Rmaj
    Zxy[i,*] = Zpos
  ENDFOR
  
  ; Curvature
  bxcvx = DBLARR(nrad, npar)
  bxcvy = bxcvx
  bxcvz = bxcvx
  FOR i=0, nrad-1 DO BEGIN
    bxcvx[i,*] = bxcvx1d
    bxcvy[i,*] = bxcvy1d
    bxcvz[i,*] = bxcvz1d
  ENDFOR
  
  ; Angle for twist-shift as function of x
  ; Taylor expanded q
  ShiftAngle = qloop + rpsi * sloop
  
  ; Integrated shear
  sinty2 = DBLARR(nrad, npar)
  FOR i=0, nrad-1 DO sinty2[i,*] = sinty

  ; Toroidal shift for shifted radial derivatives (optional)
  ; As with twist-shift, this is a linear expansion
  zShift = DBLARR(nrad, npar)
  FOR i=0, nrad-1 DO BEGIN
    zShift[i,*] = qinty-qinty[rmaxi] + rpsi[i]*sinty
  ENDFOR
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Save to file
  ; NOTE: This needs the physics initialisation code
  ;       to calculate the metric tensor components

  handle = file_open(output, /CREATE)
  
  ; Size of the grid
  s = file_write(handle, "nx", nrad)
  s = file_write(handle, "ny", npar)
  
  ; Grid spacing
  s = file_write(handle, "dx", dx)
  s = file_write(handle, "dy", dy)

  ; topology
  ixseps = nrad+1 ; Separatrix outside domain -> Periodic with twist-shift
  s = file_write(handle, "ixseps1", ixseps)
  s = file_write(handle, "ixseps2", ixseps)
  
  ; Shift angles
  s = file_write(handle, "ShiftAngle", ShiftAngle) ; For twist-shift location
  s = file_write(handle, "zShift", zShift) ; for shifted radial derivatives
  s = file_write(handle, "sinty", sinty2)

  ; Geometric quantities
  s = file_write(handle, "Rxy",  Rxy)
  s = file_write(handle, "Zxy",  Zxy) ; Not needed for simulation, useful for plotting
  s = file_write(handle, "hthe", hxy)
  
  ; Field components
  s = file_write(handle, "Bpxy", Bpxy)
  s = file_write(handle, "Btxy", Btxy)
  s = file_write(handle, "Bxy",  Bxy)

  ; Plasma profiles
  s = file_write(handle, "pressure", pressure)
  s = file_write(handle, "Jpar0", Jpar0)
  
  ; Typical quantities
  rmag = Rmaj[rmaxi]   ; maximum major radius
  bmag = ABS(B[rmaxi]) ; B field at the same location
  s = file_write(handle, "bmag", bmag)
  s = file_write(handle, "rmag", rmag)
  
  ; Curvature
  s = file_write(handle, "bxcvx", bxcvx)
  s = file_write(handle, "bxcvy", bxcvy)
  s = file_write(handle, "bxcvz", bxcvz)
  
  ; Psi normalisation (only for post-processing)
  s = file_write(handle, "psi_axis", psi_axis)
  s = file_write(handle, "psi_bndry", psi_sep)
  
  file_close, handle
  
  STOP
END
