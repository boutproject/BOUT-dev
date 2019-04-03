; flux tube input file in the SOL of a tokamak equilibrium
;
; Inputs
;
;   gfile   [string]             Name of the file to read
;   psinorm [double, optional]   Normalised psi of the flux surface
;                                psinorm = (psi - psi_axis)/(psi_sep - psi_axis)
;
; Keywords
;   output [string]       Name of the output file
;   nx [int]              Number of radial grid points
;   ny [int]              Number of points along field-line
;   psiwidth [double]     Radial width of the box in normalised psi
;   /equ  [true/false]    Force input file to be a .equ file. Normally
;                            goes on file ending.
;
;   wall_file [string]    File containing wall coordinates if not given in gfile
;   flip_Bt [Bool]        Set this to artificially reverse the sign of the toroidal field
;   scaleX [double]       Linearly scale the x domain. Needed for Zshift calculation
;
; Features
;   Uses derivatives along field line curve to calculate curvature and 
;   magnetic shear
;
;   Calculates all quantities prior to interpolation and beyond the target plates to 
;   eliminate errors from reduced smoothing near the target
;
; Warnings
;   The code cannot currently facilitate tokamak equilibria with negative
;   poloidal field. 
;
;   The zshift calculation has not been properly tested



PRO sol_flux_tube, gfile, psinorm, output=output, nx=nx, ny=ny, psiwidth=psiwidth, equ=equ, wall_file=wall_file,flip_Bt=flip_Bt,scaleX=scaleX
	

  IF N_PARAMS() EQ 0 THEN BEGIN
    PRINT, "Arguments are: gfile [, psinorm]"
    RETURN
  ENDIF ELSE IF N_PARAMS() EQ 1 THEN BEGIN
    ; No psinorm
    psinorm = 1.05D
  ENDIF

  IF NOT KEYWORD_SET(output) THEN output="fluxtube"+STR(psinorm)+".grd.nc"
  
  IF NOT KEYWORD_SET(nx) THEN nx = 132
  IF NOT KEYWORD_SET(ny) THEN ny = 128

  IF NOT KEYWORD_SET(psiwidth) THEN psiwidth = 0.05D

  IF psinorm LE 1.0D THEN BEGIN
    PRINT, "Error: Normalised psi must be greater than 1"
    RETURN
  ENDIF
  
  ; Get the file extension
  s = STRSPLIT(gfile, '.', /extract)
  IF (STRLOWCASE(s[N_ELEMENTS(s)-1]) EQ 'equ') OR KEYWORD_SET(equ) THEN BEGIN
    ; Either file ends in '.equ' or the keyword was set
    
    g = read_equ(gfile)
    equ = 1 
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
  
  ; Add these if they're not already in the structure
  str_check_present, rzgrid, "simagx", psi_axis
  str_check_present, rzgrid, "sibdry", psi_sep
  
  ; See how this compares to grid values
  PRINT, "Psi Axis : ", psi_axis, rzgrid.simagx
  PRINT, "Psi Bndry: ", psi_sep, rzgrid.sibdry
  
  psi_axis = rzgrid.simagx
  psi_sep = rzgrid.sibdry
  psi = psi_axis + psinorm*(psi_sep - psi_axis)
  print, "Psi :", psi


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

  rpos = INTERPOLATE(rzgrid.R, ri, /DOUBLE)
  zpos = INTERPOLATE(rzgrid.Z, zi, /DOUBLE)


  ;Check that indexing is assending in poloidal angle
  ;Poloidal angle should go clockwise ie from +ve Z to -ve Z
  if zpos[0] LT zpos[N_ELEMENTS(zpos)-1] THEN BEGIN
	;Need to reverse indices
	print,"Reversing indices"
	dummy = DBLARR(n_elements(ri))
	for i=0,n_elements(ri)-1 do begin
		dummy[i] = ri[n_elements(ri)-1-i]
	ENDFOR
	ri = dummy

        dummy = DBLARR(n_elements(zi))
        for i=0,n_elements(zi)-1 do begin
                dummy[i] = zi[n_elements(zi)-1-i]
        ENDFOR
        zi = dummy
  rpos = INTERPOLATE(rzgrid.R, ri, /DOUBLE)
  zpos = INTERPOLATE(rzgrid.Z, zi, /DOUBLE)
 ENDIF

  ; Smooth positions
  rpos = SMOOTH(rpos, 4)
  zpos = SMOOTH(zpos, 4)

  
  IF rzgrid.nlim GT 2 THEN BEGIN
    ; Find intersections with boundary
    oplot, rzgrid.rlim, rzgrid.zlim, thick=2, color=3
    
    
    cpos = line_crossings(rpos, zpos, 0, $
                          rzgrid.rlim, rzgrid.zlim, 1, $
                          ncross=ncross, inds1=inds)
    
    print, "Number of crossings: ", ncross
    
    IF ncross NE 2 THEN BEGIN
      PRINT, "HELP! Don't know what to do..."
      STOP
    ENDIF 
    
  ENDIF ELSE BEGIN
    PRINT, "WARNING: No boundary found, please enter boundary indices: "
    Print, "Total no points: ",n_elements(rpos)
    inds = DBLARR(2)
    inds_ok = 'N'
    oplot, rpos, zpos, color=4, thick=2
    IF KEYWORD_SET(wall_file) THEN BEGIN ;Plot wall
    W = file_import("SXDWall.nc")
    oplot,W.nlimr,W.nlimz_lower,thick=2
    oplot,W.nlimr,W.nlimz_upper,thick=2
    ENDIF

    while inds_ok ne 'Y' DO BEGIN
    	READ, inds0, PROMPT="Lower index: "
    	READ, inds1, PROMPT="Upper index: "
	oplot,[rpos[inds0],rpos[inds1]],[zpos[inds0],zpos[inds1]],psym=1,thick=4
	READ,inds_ok, PROMPT="Is this ok (Y/N)?"
 ENDWHILE 
 inds[0]=inds0
 inds[1]=inds1
 ENDELSE
    

  oplot,rpos[inds[0]:inds[1]],zpos[inds[0]:inds[1]],thick=2,color=4  
  ;;;;;;;;; Construct arc lengths
  
  drdi = DERIV(rpos)
  drdi = SMOOTH(drdi,4)
  dzdi = DERIV(zpos)
  dzdi = smooth(dzdi,4)
  dldi = SQRT(drdi^2 + dzdi^2) ; Length differential
  dldi = smooth(dldi,4) 
  
  L = int_func(dldi, /simple) ; Field line length in poloidal plane



  ;;;;;;;;;;; Toroidal field
  ; f = RBt and f'
  ngrid = n_elements(rzgrid.fpol)
  psigrid = psi_axis + (psi_sep - psi_axis)*FINDGEN(ngrid)/DOUBLE(ngrid)

  f = INTERPOLATE(rzgrid.fpol, psinorm*ngrid, /DOUBLE)
  ;Term not included in .equ file, needs checking
  IF NOT KEYWORD_SET(equ) THEN dfdpsi = INTERPOLATE(DERIV(psigrid,rzgrid.fpol), psinorm*ngrid, /DOUBLE)
  
  ;;;;;;;;;;; Poloidal field

  npoints = N_ELEMENTS(ri)
  Bpol = DBLARR(npoints)
  drposdpsi = DBLARR(npoints)
  dzdpsi = DBLARR(npoints)
  dBpoldpsi = DBLARR(npoints)
  dBpdz = DBLARR(npoints)
  dBpdr = DBLARR(npoints)

  IF NOT KEYWORD_SET(equ) THEN dfdpsi = dfdpsi
  FOR i=0, npoints-1 DO BEGIN
     ; Get local gradient of psi
     grad = EvalCosP(dctpsi,x0 = ri[i],y0 = zi[i])

     ;Set length elements for normalization  
     dr = rzgrid.r[1]-rzgrid.r[0]
     dz = rzgrid.z[1]-rzgrid.z[0]      

     ;Get differentials 
     
     dpsidr = grad[1]
     dpsidr = dpsidr /dr
     dpsidz = grad[2]
     dpsidz = dpsidz / dz
     d2psidr2 = grad[3]
     d2psidr2 = d2psidr2 / dr^2
     d2psidz2 = grad[4]
     d2psidz2 = d2psidz2 / dz^2
     d2psidrdz = grad[5]
     d2psidrdz = d2psidrdz / (dr*dz)

     ;Calculate Poloidal B terms
     Bpol[i] = SQRT(dpsidr^2 + dpsidz^2) / (rpos[i])   
 
     ;dr/dpsi
     drposdpsi[i] = dpsidr/(rpos[i]^2*Bpol[i]^2)

     ;dz/dpsi
     dzdpsi[i] = dpsidz/(rpos[i]^2*Bpol[i]^2)

     ;dBp/dz
     dBpdz[i] = (1/(Bpol[i]*(rpos[i]^2)))*(dpsidz*d2psidz2 + dpsidr* d2psidrdz)

     ;dBp/dr
     dBpdr[i] = (1/(Bpol[i]*(rpos[i]^2)))*(dpsidr*d2psidr2 + dpsidz* d2psidrdz) - Bpol[i]/rpos[i]

     
  ENDFOR
  
  drposdpsi = SMOOTH(drposdpsi,4)
  dzdpsi = SMOOTH(dzdpsi,4)
  
  dBpoldpsi = dBpdr*drposdpsi + dBpdz*dzdpsi
  dBpoldpsi = SMOOTH(dBpoldpsi,4)

  ;Toroidal field
  Btor =  f / rpos
  
  IF KEYWORD_SET(equ) THEN BEGIN
	dBtordpsi =   -(f/rpos^2)*drposdpsi
  ENDIF ELSE BEGIN dBtordpsi =  - (f/rpos^2)*drposdpsi + (1/rpos)*dfdpsi
  ENDELSE  ;<--- include this term if not using equ file
  dBtordpsi = SMOOTH(dBtordpsi,4)
  B = SQRT(Btor^2 + Bpol^2)
 
  IF KEYWORD_SET(flip_Bt) THEN BEGIN
	Btor = -Btor
	dBtordpsi = -dBtordpsi
  ENDIF

  ;Recalculate Grad psi
	
  dpsidR = drposdpsi*(rpos^2*Bpol^2) 
  dpsidZ = dzdpsi*rpos^2*Bpol^2 
 
  
  dBdpsi = (Bpol*dBpoldpsi + Btor*dBtordpsi)/B
  dBdR = dBdpsi*dpsidR
  dBdZ = dBdpsi*dpsidZ
 
  ; Change in toroidal angle, following a field-line
  dtdi = dldi * Btor / (Bpol * rpos)  
  qinty = int_func(dtdi, /simple)
	 
  ; Get parallel distance
  dsdi = dldi * B / Bpol
  
  ; Parallel distance
  s = int_func(dsdi, /simple)
  
  ;Calculate hthe, ensuring that theta = 0,2*pi is at the divertor targets
  hthe = (L[inds[1]] - L[inds[0]])/(2*!DPi)  
  print,"hthe = ",hthe

  ;;;;;;;;;;;;;;;;;;;; CURVATURE ;;;;;;;;;;;;;;;;;;;;;;;

  ; Calculate b x kappa using coordinates of a single field-line
  
  ; First derivatives along the line to get tangent vector
  dr = SMOOTH(DERIV(s, rpos),4)  ; R position
  dz = SMOOTH(DERIV(s, zpos),4)  ; Z position 
  dp = SMOOTH(DERIV(s, qinty),4) ; Toroidal angle  

  ; Second derivatives
  d2r = SMOOTH(DERIV(s, dr),4)
  d2z = SMOOTH(DERIV(s, dz),4)
  d2p = SMOOTH(DERIV(s, dp),4)
 
  ; Components of b (tangent vector)
  br = dr
  bz = dz
  bp = rpos*dp
  
  ; Components of curvature in cylindrical coordinates
  kr = SMOOTH(d2r - rpos*dp^2,4)
  kz = SMOOTH(d2z,4)
  kp = SMOOTH(2.D*dr*dp + rpos*d2p,4)
  
  ;Components of curvature in toroidal coordinates
  ;Not needed for calculation but useful for diagnostic purposes
    
  kpsi = kr*dpsidR + kz*dpsidZ 
  ktheta = ( kr * dpsidZ - kz * dpsidR ) / (rpos*Bpol*hthe)
  
  ;;;;;;;;;;;;;;;;; Magnetic shear;;;;;;;;;;;;;;;;;;;;;;;;;

  ;Local pitch  
  nu = hthe * Btor / (Bpol * rpos)
  nu = SMOOTH(nu,4)

  dhthedpsi = -(d2z*dr - d2r*dz)*(hthe*(B^3))/((Bpol^4)*rpos)
  dhthedpsi = SMOOTH(dhthedpsi,4)
  
  
  ;Calculate dnu/dpsi and integrate        
 
  dnudpsi = nu*((1/Btor)*dBtordpsi - drposdpsi*(1/rpos) - (1/Bpol)*dBpoldpsi + (1/hthe)*dhthedpsi);
  dnudpsi = SMOOTH(dnudpsi,4)
 
  sinty = int_func(L/(hthe),dnudpsi,/simple)
  
  ; Now we have sinty, continue cuvature calculation

  ; Calculate bxk in cylindrical coordinates

  bxkr = bp*kz - bz*kp
  bxkz = br*kp - bp*kr
  bxkp = bz*kr - br*kz

  ; Calculate components in (psi, theta, phi) toroidal coordinates

  ; (bxk) dot grad psi
  bxkpsi   = bxkr * dpsidR + bxkz * dpsidZ
  ; (bxk) dot grad theta
  bxktheta = ( bxkr * dpsidZ - bxkz * dpsidR ) / (rpos*Bpol*hthe)
   
  ;Finally into field-aligned (ballooning) coordinates
  bxcvx1d = SMOOTH(bxkpsi,4)
  bxcvy1d = SMOOTH(bxktheta,4)
  bxcvz1d = SMOOTH(bxkp,4); add the rest to bxcvz later: -sinty*bxkpsi - nu*bxktheta

  ;;;;;;;;;;;;;;;; RADIAL MESH ;;;;;;;;;;;;;;;;;;
  
  ; Convert normalised psi to psi
  dpsi = (psiwidth * (psi_sep - psi_axis)) / DOUBLE(nx-1)
 
  ;;;;;;;;;;;;;;;; INTERPOLATE ALL QUANTITIES ONTO FIELD LINE ;;;;;

  ;Pick out quantities between the two targets
  rpos = rpos[inds[0]:inds[1]]
  zpos = zpos[inds[0]:inds[1]]
  B = B[inds[0]:inds[1]]
  Btor = Btor[inds[0]:inds[1]]
  Bpol = Bpol[inds[0]:inds[1]]
  dBdpsi = dBdpsi[inds[0]:inds[1]]
  dBdR = dBdR[inds[0]:inds[1]]
  dBdZ = dBdZ[inds[0]:inds[1]]  
  nu = nu[inds[0]:inds[1]]
  sinty = sinty[inds[0]:inds[1]]
  dnudpsi = dnudpsi[inds[0]:inds[1]]
  dpsidR = dpsidR[inds[0]:inds[1]] 
  dpsidZ = dpsidZ[inds[0]:inds[1]]
  qinty = qinty[inds[0]:inds[1]]
  bxcvx1d = bxcvx1d[inds[0]:inds[1]]
  bxcvy1d = bxcvy1d[inds[0]:inds[1]]
  kpsi = kpsi[inds[0]:inds[1]]
  ktheta = ktheta[inds[0]:inds[1]]
  s = s[inds[0]:inds[1]]
  bxcvz1d = bxcvz1d[inds[0]:inds[1]]
  dldi = dldi[inds[0]:inds[1]]
  kphi = kp[inds[0]:inds[1]]
 
  L = int_func(dldi, /simple) 
  lpos = max(L) * FINDGEN(ny)/DOUBLE(ny-1)

  ;Interpolate onto grid equally spaced in poloidal angle
 
  inds = INTERPOL(findgen(N_ELEMENTS(L)), L, lpos)
  rpos = INTERPOLATE(rpos, inds, /DOUBLE)
  zpos = INTERPOLATE(zpos, inds, /DOUBLE)
  s = INTERPOLATE(s,inds, /DOUBLE)
  B = INTERPOLATE(B, inds, /DOUBLE)
  Btor = INTERPOLATE(Btor, inds, /DOUBLE)
  Bpol = INTERPOLATE(Bpol, inds, /DOUBLE)
  nu = INTERPOLATE(nu, inds, /DOUBLE)
  sinty = INTERPOLATE(sinty, inds, /DOUBLE)
  sinty = sinty - sinty[ny/2] ; take theta_0 at outboard midplane
  dnudpsi =  INTERPOLATE(dnudpsi, inds, /DOUBLE)
  dpsidR =  INTERPOLATE(dpsidR, inds, /DOUBLE)
  dpsidZ = INTERPOLATE(dpsidZ,inds, /DOUBLE)
  hthe = MAX(lpos)/(2*!DPi)
  qinty =  INTERPOLATE(qinty, inds, /DOUBLE)
  qinty = qinty - qinty[FLOOR(ny/2)]
  dBdpsi = INTERPOLATE(dBdpsi,inds, /DOUBLE)
  dBdR = INTERPOLATE(dBdR,inds, /DOUBLE)  
  dBdZ = INTERPOLATE(dBdZ,inds, /DOUBLE)
  

  ;Add in missing terms in curvature
  bxcvx1d =  INTERPOLATE(bxcvx1d, inds, /DOUBLE)
  bxcvy1d =  INTERPOLATE(bxcvy1d, inds, /DOUBLE)
  bxcvz1d =  INTERPOLATE(bxcvz1d, inds, /DOUBLE) - (sinty*bxcvx1d + nu*bxcvy1d)

  kpsi = INTERPOLATE(kpsi,inds, /DOUBLE)
  ktheta = INTERPOLATE(ktheta,inds, /DOUBLE)
  kphi = INTERPOLATE(kphi,inds, /DOUBLE)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Put everything into 2D arrays
  
  ; B field components
  Bpxy = DBLARR(nx, ny)
  Btxy = DBLARR(nx, ny)
  Bxy  = DBLARR(nx, ny)
  FOR i=0, nx-1 DO BEGIN
     Bpxy[i,*] = Bpol
     Btxy[i,*] = Btor
     Bxy[i,*] = B
  ENDFOR
  
  ; Grid spacing
  dx = DBLARR(nx, ny) + dpsi
  dy = DBLARR(nx, ny) + 2.D*!DPI/DOUBLE(ny)
  
  ; Geometrical quantities
  hxy = DBLARR(nx, ny)
  Rxy = DBLARR(nx, ny)
  Zxy = DBLARR(nx, ny)
  
  FOR i=0, nx-1 DO BEGIN
    hxy[i,*] = hthe
    Rxy[i,*] = rpos
    Zxy[i,*] = zpos
  ENDFOR
  
  ; Curvature and other quantities
  bxcvx = DBLARR(nx, ny)
  bxcvy = bxcvx
  bxcvz = bxcvx
  sinty2 = bxcvx
  vxy = bxcvx
  sxy = bxcvx 
  S_parxy = bxcvx
  zshift = bxcvx
  kpsi_1 = bxcvx
  ktheta_1 = bxcvx
IF NOT KEYWORD_SET(scaleX) THEN BEGIN
	print,"WARNING: No scaleX set, setting to 1"
	scaleX = 1
ENDIF
  FOR i=0, nx-1 DO BEGIN
     bxcvx[i,*] = bxcvx1d
     bxcvy[i,*] = bxcvy1d
     bxcvz[i,*] = bxcvz1d
     sinty2[i,*] = sinty
     vxy[i,*] = nu    
     sxy[i,*] = dnudpsi
     zshift[i,*] = qinty + (i*scaleX*dpsi)*sinty ;First order Taylor Expansion for zshift
     S_parxy[i,*] = s
     kpsi_1[i,*] = kpsi
     ktheta_1[i,*] = ktheta
     
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
 
  ;Curvature vector components, not needed for physics file, but useful 
  s = file_write(handle, "kpsi", kpsi_1)
  s = file_write(handle, "ktheta", ktheta_1)

  ; Grid spacing
  s = file_write(handle, "dx", dx)
  s = file_write(handle, "dy", dy)

  ; topology
  ixseps = -1 ; Separatrix inside domain -> not periodic
  s = file_write(handle, "ixseps1", ixseps)
  s = file_write(handle, "ixseps2", ixseps)
  
  ; Shift angles

  s = file_write(handle, "zShift", zshift) ; for shifted radial derivatives
  s = file_write(handle, "sinty", sinty2)  ; Integrates shear
  s = file_write(handle, "vxy", vxy)       ; Local field line pitch
  s = file_write(handle, "sxy", sxy)       ; Local magnetic shear

  ; Geometric quantities
  s = file_write(handle, "Rxy",  Rxy) ; Major radius
  s = file_write(handle, "Zxy",  Zxy) ; Not needed for simulation, useful for plotting
  s = file_write(handle, "hthe", hxy) 
  s = file_write(handle, "Lxy",S_parxy) ;Distance along magnetic fied line
  
  ; Field components
  s = file_write(handle, "Bpxy", Bpxy)
  s = file_write(handle, "Btxy", Btxy)
  s = file_write(handle, "Bxy",  Bxy)
  
  ; Typical quantities, useful for normalization
  rmag = max(abs(rpos))   ; maximum major radius
  bmag = max(abs(B)) ; B field at the same location
  s = file_write(handle, "bmag", bmag)
  s = file_write(handle, "rmag", rmag)
  
  ; b x kappa
  s = file_write(handle, "bxcvx", bxcvx)
  s = file_write(handle, "bxcvy", bxcvy)
  s = file_write(handle, "bxcvz", bxcvz)

   
  ; Psi normalisation (only for post-processing)
  s = file_write(handle, "psi_axis", psi_axis)
  s = file_write(handle, "psi_bndry", psi_sep)
  
  file_close, handle
  
END

