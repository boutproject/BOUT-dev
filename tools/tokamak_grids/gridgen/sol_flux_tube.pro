;
; Generate a flux tube input file in the SOL of a tokamak equilibrium
;
; Inputs
;
;   gfile   [string]             Name of the file to read
;   psinorm [float, optional]    Normalised psi of the flux surface
;
; Keywords
;   output [string]       Name of the output file
;   nx [int]              Number of radial grid points
;   ny [int]              Number of points along field-line
;   psiwidth [float]      Radial width of the box in normalised psi
;   /equ  [true/false]    Force input file to be a .equ file. Normally
;                            goes on file ending.
; 

PRO sol_flux_tube, gfile, psinorm, output=output, nx=nx, ny=ny, psiwidth=psiwidth, equ=equ
  IF N_PARAMS() EQ 0 THEN BEGIN
    PRINT, "Arguments are: gfile [, psinorm]"
    RETURN
  ENDIF ELSE IF N_PARAMS() EQ 1 THEN BEGIN
    ; No psinorm
    psinorm = 1.05
  ENDIF

  IF NOT KEYWORD_SET(output) THEN output="fluxtube"+STR(psinorm)+".grd.nc"
  
  IF NOT KEYWORD_SET(nx) THEN nx = 132
  IF NOT KEYWORD_SET(ny) THEN ny = 128

  IF NOT KEYWORD_SET(psiwidth) THEN psiwidth = 0.05

  IF psinorm LE 1.0 THEN BEGIN
    PRINT, "Error: Normalised psi must be greater than 1"
    RETURN
  ENDIF
  
  ; Get the file extension
  s = STRSPLIT(gfile, '.', /extract)
  IF (STRLOWCASE(s[N_ELEMENTS(s)-1]) EQ 'equ') OR KEYWORD_SET(equ) THEN BEGIN
    ; Either file ends in '.equ' or the keyword was set
    
    g = read_equ(gfile)
    
    ; Check the return
    IF SIZE(g, /TYPE) NE 8 THEN BEGIN
      PRINT, "ERROR: Couldn't read EQU file '"+gfile+"'"
      RETURN
    ENDIF
    
    rzgrid = {nr:g.nx, nz:g.ny, $ ; Number of grid points
              r:g.r, z:g.z, $ ; R and Z as 1D arrays
              psi:g.psi, $ ; Poloidal flux in Weber/rad on grid points
              fpol:[g.fpol], $
              nlim:0} ; Wall boundary
  ENDIF ELSE BEGIN
    ; Assume it's an EFIT G-EQDSK file
    PRINT, "Reading a G-EQDSK file"
    g = read_neqdsk(gfile)

    ; Check the return
    IF SIZE(g, /TYPE) NE 8 THEN BEGIN
      PRINT, "ERROR: Couldn't read G-EQDSK file '"+gfile+"'"
      RETURN
    ENDIF
    
    rzgrid = {nr:g.nx, nz:g.ny, $ ; Number of grid points
              r:REFORM(g.r[*,0]), z:REFORM(g.z[0,*]), $ ; R and Z as 1D arrays
              simagx:g.simagx, sibdry:g.sibdry, $       ; Range of psi
              psi:g.psi, $ ; Poloidal flux in Weber/rad on grid points
              pres:g.pres, $ ; Plasma pressure in nt/m^2 on uniform flux grid
              fpol:g.fpol, $
              qpsi:g.qpsi, $     ; q values on uniform flux grid
              nlim:g.nlim, rlim:g.xlim, zlim:g.ylim} ; Wall boundary
  ENDELSE

  ; Plot psi
  nlev = 100
  minf = MIN(rzgrid.psi)
  maxf = MAX(rzgrid.psi)
  levels = findgen(nlev)*(maxf-minf)/FLOAT(nlev-1) + minf
  
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
  
  ; Add these if they're not already in the structure
  str_check_present, rzgrid, "simagx", psi_axis
  str_check_present, rzgrid, "sibdry", psi_sep
  
  ; See how this compares to grid values
  PRINT, "Psi Axis : ", psi_axis, rzgrid.simagx
  PRINT, "Psi Bndry: ", psi_sep, rzgrid.sibdry
  
  psi_axis = rzgrid.simagx
  psi_sep = rzgrid.sibdry

  psi = psi_axis + psinorm*(psi_sep - psi_axis)
  
  ;;;;;;;;;;;;;;; Calculate DCT ;;;;;;;;;;;;;;
  
  PRINT, "Calculating DCT..."
  dctpsi = DCT2D(rzgrid.psi)
  PRINT, "Finished DCT"
  
  ; Create a structure containing interpolation settings and data
  interp_data = {nx:rzgrid.nr, ny:rzgrid.nz, $
                 method:0, $
                 f: rzgrid.psi, $       ; Always include function
                 dct: dctpsi} ; Pass the DCT coefficients
  
  ;;;;;;;;;;;;;;; Get a contour line
  
  contour_lines, rzgrid.psi, findgen(rzgrid.nr), findgen(rzgrid.nz), $
                 levels=[psi], $
                 path_info=info, path_xy=xy
  
  IF N_ELEMENTS(info) GT 1 THEN BEGIN
     ; Find the surface closest to the outboard midplane
      
     ind = closest_line(info, xy, critical.opt_ri[critical.primary_opt], critical.opt_zi[critical.primary_opt])
     info = info[ind]
  ENDIF ELSE info = info[0]
  
  ; Extract r and z indices
  ri = REFORM(xy[0,info.offset:(info.offset+info.n-1)])
  zi = REFORM(xy[1,info.offset:(info.offset+info.n-1)])
  np = N_ELEMENTS(ri)

  ; Refine flux surface location using DCT
  FOR i=0, np-1 DO BEGIN
    follow_gradient, interp_data, rzgrid.R, rzgrid.Z, ri[i], zi[i], psi, $
                     rinew, zinew, status=status
    IF status THEN BEGIN
      PRINT, "ERROR refining surface location"
      RETURN
    ENDIF
    
    ri[i] = rinew
    zi[i] = zinew
  ENDFOR
  
  rpos = INTERPOLATE(rzgrid.R, ri)
  zpos = INTERPOLATE(rzgrid.Z, zi)

  ; Smooth positions
  rpos = SMOOTH(rpos, 4)
  zpos = SMOOTH(zpos, 4)
  
  IF rzgrid.nlim GT 2 THEN BEGIN
    ; Find intersections with boundary
    oplot, rzgrid.xlim, rzgrid.ylim, thick=2, color=3
    
    
    cpos = line_crossings(rpos, zpos, 0, $
                          rzgrid.xlim, rzgrid.ylim, 1, $
                          ncross=ncross, inds1=inds)
    
    print, "Number of crossings: ", ncross
    
    IF ncross NE 2 THEN BEGIN
      PRINT, "HELP! Don't know what to do..."
      STOP
    ENDIF
    
    rpos = rpos[inds[0]:inds[1]]
    zpos = zpos[inds[0]:inds[1]]
    
    ri = ri[inds[0]:inds[1]]
    zi = zi[inds[0]:inds[1]]
    
  ENDIF ELSE BEGIN
    PRINT, "WARNING: No boundary found"
  ENDELSE
    
  oplot, rpos, zpos, color=4, thick=2
  
  ;;;;;;;;; Construct constant arc lengths
  
  drdi = DERIV(rpos)
  dzdi = DERIV(zpos)
  dldi = SQRT(drdi^2 + dzdi^2) ; Length differential
  
  L = int_func(dldi, /simple) ; Integrate function
  
  lpos = max(L) * FINDGEN(ny)/FLOAT(ny-1)
  dldi = max(L) / FLOAT(ny-1)
  
  ; Find index
  inds = INTERPOL(findgen(N_ELEMENTS(L)), L, lpos)
  
  rpos = INTERPOLATE(rpos, inds)
  zpos = INTERPOLATE(zpos, inds)
  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Now got regularly spaced grid points along flux surface
  ; (rpos, zpos)   =   ( Major radius, Height )

  ;;;;;;;;;;; Toroidal field
  ; f = RBt
  
  f = rzgrid.fpol[n_elements(rzgrid.fpol)-1]
  
  Btor = f / rpos
  
  ;;;;;;;;;;; Poloidal field
  
  Bpol = FLTARR(ny)
  
  FOR i=0, ny-1 DO BEGIN
     ; Get local gradient of psi
     local_gradient, interp_data, ri[i], zi[i], status=status, dfdr=dpsidr, dfdz=dpsidz
     
     ; Convert to correct units
     dpsidr = dpsidr / (rzgrid.r[1] - rzgrid.r[0])
     dpsidz = dpsidz / (rzgrid.z[1] - rzgrid.z[0])
     
     Bpol[i] = SQRT(dpsidr^2 + dpsidz^2) / rpos[i]
     
  ENDFOR
  
  ; Total field
  B = SQRT(Btor^2 + Bpol^2)

  ; Change in toroidal angle, following a field-line
  dtdi = dldi * Btor / (Bpol * rpos)
  
  qinty = int_func(dtdi, /simple)
  
  ; Get parallel distance
  dsdi = dldi * B / Bpol
  
  ; Parallel distance
  s = int_func(dsdi, /simple)
  
  hthe = dldi / (2.*!PI / FLOAT(ny))

  ;;;;;;;;;;;;;;;;; Magnetic shear
  
  ; Tricky... zero for now
  
  pitch = hthe * Btor / (Bpol * Rpos)

  sinty = FLTARR(ny)
  
  ;;;;;;;;;;;;;;;;;;;; CURVATURE ;;;;;;;;;;;;;;;;;;;;;;;
  ; Calculate b x kappa using coordinates of a single field-line
  
  ; First derivatives along the line to get tangent vector
  dr = SMOOTH(DERIV(s, rpos),3)  ; R position
  dz = SMOOTH(DERIV(s, zpos),3)  ; Z position 
  dp = SMOOTH(DERIV(s, qinty),3) ; Toroidal angle
  
  ; Second derivatives
  d2r = SMOOTH(DERIV(s, dr),3)
  d2z = SMOOTH(DERIV(s, dz),3)
  d2p = SMOOTH(DERIV(s, dp),3)
  
  ; Components of b (tangent vector)
  br = dr
  bz = dz
  bp = Rpos*dp

  ; Components of curvature (unit vectors)
  kr = d2r - rpos*dp^2
  kz = d2z
  kp = 2.*dr*dp + rpos*d2p

  ; Calculate bxk in cylindrical coordinates

  bxkr = bp*kz - bz*kp
  bxkz = br*kp - bp*kr
  bxkp = bz*kr - br*kz
  
  ; Calculate components in (psi, theta, phi) toroidal coordinates
  ; (bxk) dot grad psi
  bxkpsi   = bxkr * dpsidR + bxkz * dpsidZ
  ; (bxk) dot grad theta
  bxktheta = ( bxkr * dpsidZ - bxkz * dpsidR ) / (rpos*Bpol*hthe)

  ;Finally into field-aligned coordinates
  bxcvx1d = SMOOTH(bxkpsi,3)
  bxcvy1d = SMOOTH(bxktheta,3)
  bxcvz1d = SMOOTH(bxkpsi - sinty*bxkpsi - pitch*bxktheta,3)

  ;;;;;;;;;;;;;;;; RADIAL MESH ;;;;;;;;;;;;;;;;;;
  
  ; Convert normalised psi to psi
  dpsi = (psiwidth * (psi_sep - psi_axis)) / FLOAT(nx-1)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Put everything into 2D arrays
  
  ; B field components
  Bpxy = FLTARR(nx, ny)
  Btxy = FLTARR(nx, ny)
  Bxy  = FLTARR(nx, ny)
  FOR i=0, nx-1 DO BEGIN
     Bpxy[i,*] = Bpol
     Btxy[i,*] = Btor
     Bxy[i,*] = B
  ENDFOR
  
  ; Grid spacing
  dx = FLTARR(nx, ny) + dpsi
  dy = FLTARR(nx, ny) + 2.*!PI/FLOAT(ny)
  
  ; Geometrical quantities
  hxy = FLTARR(nx, ny)
  Rxy = FLTARR(nx, ny)
  Zxy = FLTARR(nx, ny)
  FOR i=0, nx-1 DO BEGIN
    hxy[i,*] = hthe
    Rxy[i,*] = rpos
    Zxy[i,*] = zpos
  ENDFOR
  
  ; Curvature
  bxcvx = FLTARR(nx, ny)
  bxcvy = bxcvx
  bxcvz = bxcvx
  FOR i=0, nx-1 DO BEGIN
     bxcvx[i,*] = bxcvx1d
     bxcvy[i,*] = bxcvy1d
     bxcvz[i,*] = bxcvz1d
  ENDFOR

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Save to file
  ; NOTE: This needs the physics initialisation code
  ;       to calculate the metric tensor components
  
  PRINT, "Writing grid file to: "+output
  
  handle = file_open(output, /CREATE)
  
  ; Size of the grid
  s = file_write(handle, "nx", nx)
  s = file_write(handle, "ny", ny)
  
  ; Grid spacing
  s = file_write(handle, "dx", dx)
  s = file_write(handle, "dy", dy)

  ; topology
  ixseps = -1 ; Separatrix inside domain -> not periodic
  s = file_write(handle, "ixseps1", ixseps)
  s = file_write(handle, "ixseps2", ixseps)
  
  ; Shift angles
  ;s = file_write(handle, "ShiftAngle", ShiftAngle) ; For twist-shift location
  ;s = file_write(handle, "zShift", zShift) ; for shifted radial derivatives
  ;s = file_write(handle, "sinty", sinty2)

  ; Geometric quantities
  s = file_write(handle, "Rxy",  Rxy)
  s = file_write(handle, "Zxy",  Zxy) ; Not needed for simulation, useful for plotting
  s = file_write(handle, "hthe", hxy)
  
  ; Field components
  s = file_write(handle, "Bpxy", Bpxy)
  s = file_write(handle, "Btxy", Btxy)
  s = file_write(handle, "Bxy",  Bxy)
  
  ; Typical quantities
  rmag = max(abs(rpos))   ; maximum major radius
  bmag = max(abs(B)) ; B field at the same location
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
  
END
