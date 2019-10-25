; Takes the original R-Z data (from G-EQDSK), and the flux mesh
; from create_grid.pro
;
; o Derives additional quantities on the mesh
; o Enforces force-balance
; o Calculates integrated quantities for field-aligned codes
;
; Inputs
; ------
; 
; rz_grid - a structure containing
;    npsigrid - 1D normalised psi grid
;    fpol     - Poloidal current function
;    pres     - Plasma pressure in nt/m^2
;    qpsi     - q values
; 
; mesh - Structure produced by create_grid.pro
;
;
; Keywords
; --------
;
; poorquality (output) - set to 1 if grid is poor quality
; gui         (input switch) - If set, uses dialogs to question users

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Average over flux-surfaces
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION surface_average, var, mesh
  f = var
  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface_hypnotoad(last=last, xi=xi)
    f[xi,yi] = MEAN(var[xi,yi]) ; Average over this surface
  ENDREP UNTIL last
  RETURN, f
END

function calc_angle, x1, x2, x3, y1, y2, y3
; calculates angle associated with point 1 given 3 points (law of cosines)
	; length of segements between points
	P12 = sqrt((x1-x2)^2 + (y1-y2)^2)
	P13 = sqrt((x1-x3)^2 + (y1-y3)^2)
	P23 = sqrt((x2-x3)^2 + (y2-y3)^2)
	;calculate angle
	return, acos((P12^2 + P13^2 - P23^2)/(2.D*P12*P13))
end

function calc_beta_withgrid, r, z, x
	nx = size(r,/dimensions)
	ny = nx[1]
	nx = nx[0]

	beta = dblarr(ny)
	; loop to calc eta for non-edges
	if (x eq 0) then begin
		for j=1,ny-2 do begin
			beta[j] = calc_angle(r[x,j],r[x+1,j],r[x,j+1],z[x,j],z[x+1,j],z[x,j+1])
			beta[j] += !DPI - calc_angle(r[x,j],r[x+1,j],r[x,j-1],z[x,j],z[x+1,j],z[x,j-1])
			beta[j] /= 2.0D
		endfor
		beta[0] = calc_angle(r[x,0],r[x+1,0],r[x,1],z[x,0],z[x+1,0],z[x,1])
		beta[ny-1] = !DPI - calc_angle(r[x,ny-1],r[x+1,ny-1],r[x,ny-2],z[x,ny-1],z[x+1,ny-1],z[x,ny-2])
	endif else if (x eq nx-1) then begin
		for j=1,ny-2 do begin
			; average angle across the grid point
			beta[j] = calc_angle(r[x,j],r[x-1,j],r[x,j-1],z[x,j],z[x-1,j],z[x,j-1])
			beta[j] += !DPI - calc_angle(r[x,j],r[x-1,j],r[x,j+1],z[x,j],z[x-1,j],z[x,j+1])
			beta[j] /= 2.0D
		endfor
		beta[0] = !DPI - calc_angle(r[x,0],r[x-1,0],r[x,1],z[x,0],z[x-1,0],z[x,1])
		beta[ny-1] = calc_angle(r[x,ny-1],r[x-1,ny-1],r[x,ny-2],z[x,ny-1],z[x-1,ny-1],z[x,ny-2])
	endif else begin
		for j=1,ny-2 do begin
			; average angle across the grid point
			beta[j] = calc_angle(r[x,j],r[x-1,j],r[x,j-1],z[x,j],z[x-1,j],z[x,j-1])
			beta[j] += calc_angle(r[x,j],r[x+1,j],r[x,j+1],z[x,j],z[x+1,j],z[x,j+1])
			beta[j] += !DPI - calc_angle(r[x,j],r[x-1,j],r[x,j+1],z[x,j],z[x-1,j],z[x,j+1])
			beta[j] += !DPI - calc_angle(r[x,j],r[x+1,j],r[x,j-1],z[x,j],z[x+1,j],z[x,j-1])
			beta[j] /= 4.0D
		endfor
		beta[0] = 0.5D*calc_angle(r[x,0],r[x+1,0],r[x,1],z[x,0],z[x+1,0],z[x,1])
		beta[0] = beta[0] + 0.5D*(!DPI - calc_angle(r[x,0],r[x-1,0],r[x,1],z[x,0],z[x-1,0],z[x,1]))
		beta[ny-1] = 0.5D*calc_angle(r[x,ny-1],r[x-1,ny-1],r[x,ny-2],z[x,ny-1],z[x-1,ny-1],z[x,ny-2])
		beta[ny-1] = beta[ny-1] + 0.5D*(!DPI - calc_angle(r[x,ny-1],r[x+1,ny-1],r[x,ny-2],z[x,ny-1],z[x+1,ny-1],z[x,ny-2]))
	endelse
	
	return, beta
end

function calc_beta, Rxy, Zxy, mesh, rz_grid, method
	; calculate beta using local field line gradient

	s = size(Rxy, /dim)
	nx = s[0]
	ny = s[1]
	beta = dblarr(nx,ny)

	if(method EQ 0) then begin
		; interpolate from cartesian grid with psi values on it (ie. from efit)
		interp_data = {nx: rz_grid.nr, ny:rz_grid.nz, method:0, f:rz_grid.psi}
		
		for j=0,ny-1 do begin
			dRdr = DERIV(Rxy[*,j])
			dZdr = DERIV(Zxy[*,j])
			for i=0,nx-1 do begin
				local_gradient, interp_data, mesh.Rixy[i,j], mesh.Zixy[i,j], status=status, dfdr=dfdr, dfdz=dfdz
				dPsidR = dfdr/INTERPOLATE(DERIV(rz_grid.r),i, /DOUBLE)
				dPsidZ = dfdz/INTERPOLATE(DERIV(rz_grid.z),j, /DOUBLE)
	
				angle1 = atan(dPsidR,dPsidZ)
				angle2 = atan(dZdr[i],-dRdr[i])
	
				beta[i,j] = angle1 - angle2 - !DPI/2.D
			endfor
		endfor	

	endif else if(method EQ 1) then begin
		npol_withguards = mesh.npol + mesh.n_y_boundary_guards
		npol = round(total(npol_withguards,/cumulative))
		Nnpol = n_elements(npol)
		status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
		REPEAT BEGIN
			yi = gen_surface_hypnotoad(last=last, xi=xi, period=period)
			; find beta using one field line at xi, with y range yi
			; for better angle calculation, need to split yi into sections based on gridding
			if (xi GE mesh.nrad[0]) then begin  ; if outside the separatrix
				for i=1,Nnpol-1 do begin
					loc = where((yi GE npol[i-1]) AND (yi LT npol[i]))
					if(total(loc) NE -1) then begin
					yi_curr = yi[loc]
					beta[xi,yi_curr] = !DPI/2.D - calc_beta_withgrid(Rxy[*,yi_curr], Zxy[*,yi_curr], xi)
					endif
				endfor
				loc = where(yi LT npol_withguards[0])
				if(total(loc) NE -1) then begin
					yi_curr = yi[loc]
					beta[xi,yi_curr] = !DPI/2.D - calc_beta_withgrid(Rxy[*,yi_curr], Zxy[*,yi_curr], xi)
				endif
			endif else begin
				beta[xi,yi] = !DPI/2.D - calc_beta_withgrid(Rxy[*,yi], Zxy[*,yi], xi)
			endelse
		ENDREP UNTIL last
		;beta = smooth(beta,5) ; smooth beta, it's ugly
	endif else begin
		print,"*** ERROR: UNKNOWN METHOD FOR BETA CALCULATION ***"
		beta = 0.0D
	endelse

	return, beta

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate f = R * Bt
; Using LSODE method
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Function for use by LSODE to integrate Bt
FUNCTION bt_differential, x, Bt
  COMMON bt_com, psi, ax, bx
  a = INTERPOL(ax, psi, x)
  b = INTERPOL(bx, psi, x)
  
  RETURN, a*Bt + b/Bt
END

FUNCTION solve_f, Rxy, psixy, pxy, Bpxy, hthe
  
  MU = 4.d-7*!DPI
  
  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]
  
  a = -DDX(psixy, Rxy) / Rxy
  b = -MU*DDX(psixy, pxy) - Bpxy*DDX(Bpxy*hthe)/hthe
  
  CATCH, theError
  IF theError EQ 0 THEN BEGIN
    ; Call LSODE to follow gradient

  ENDIF ELSE BEGIN
     
  ENDELSE
END

FUNCTION force_balance, psixy, Rxy, Bpxy, Btxy, hthe, pxy
  MU =4.d-7*!DPI
  
  a = DDX(psixy, Rxy) / Rxy
  b = MU*DDX(psixy, pxy) - Bpxy*DDX(psixy, Bpxy*hthe)/hthe
  
  RETURN, DDX(psixy, Btxy) + a*Btxy + b/Btxy
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate toroidal field
; Using NEWTON method
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION Bt_func, Bt
  COMMON fnewt_com, psi, a, b
  
  RETURN, DERIV(psi, Bt) + a*Bt + b / Bt
END

FUNCTION newton_Bt, psixy, Rxy, Btxy, Bpxy, pxy, hthe, mesh
  COMMON fnewt_com, psi, a, b
  MU = 4.d-7*!DPI
  
  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]
  
  axy = DDX(psixy, Rxy) / Rxy
  bxy = MU*DDX(psixy, pxy) - Bpxy*DDX(psixy, Bpxy*hthe)/hthe
  
  Btxy2 = DBLARR(nx, ny)
  FOR i=0, ny-1 DO BEGIN
    psi = psixy[*,i]
    a = axy[*,i]
    b = bxy[*,i]
    PRINT, "Solving f for y=", i
    Btxy2[*,i] = NEWTON(Btxy[*,i], "Bt_func")
  ENDFOR
  
  ; Average f over flux surfaces
  fxy = surface_average(Btxy2*Rxy, mesh)

  RETURN, fxy / Rxy
END

function intx, Rxy, data, simple=simple

  nx = size(data,/dimensions)
  ny = nx[1]
  nx = nx[0]

  result = dblarr(nx,ny)
  result[*,*] = 0.0D
  if keyword_set(simple) then begin
    for i=0, ny-1 do begin
      for j=1, nx-1 do begin
        result[j, i] = result[j-1, i] + 0.5D*(Rxy[j, i] - Rxy[j-1, i])*(data[j, i] + data[j-1, i])
      endfor
    endfor
  endif else begin
    for i=0, ny-1 do begin
      for j=1, nx-1 do begin
        result[j,i] = int_tabulated(Rxy[0:j,i],data[0:j,i])
      endfor
    endfor
  endelse

  return, result
end

function inty, Zxy, data, simple=simple

  nx = size(data,/dimensions)
  ny = nx[1]
  nx = nx[0]

  result = dblarr(nx,ny)
  result[*,*] = 0.0D
  if keyword_set(simple) then begin
    for i=1, ny-1 do begin
      for j=0, nx-1 do begin
        result[j, i] = result[j, i-1] + 0.5D*(Zxy[j, i] - Zxy[j, i-1])*(data[j, i] + data[j, i-1])
      endfor
    endfor
  endif else begin
    for i=1, ny-1 do begin
      for j=0, nx-1 do begin
        result[j,i] = int_tabulated(Zxy[j,0:i],data[j,0:i])
      endfor
    endfor
  endelse

  return, result
end

; Integrate a function over y
FUNCTION my_int_y, var, yaxis, mesh, loop=loop, nosmooth=nosmooth, simple=simple
  f = var
  
  s = SIZE(var, /dim)
  nx = s[0]
  loop = DBLARR(nx)
  loop[*] = !VALUES.F_NAN ; Prevent accidental use of unset values
  
  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface_hypnotoad(last=last, xi=xi, period=period)
    
    f[xi,yi] = inty(yaxis[xi,yi],var[xi,yi], /simple)

    IF ~period THEN BEGIN
      ; reset integral to zero at the first grid point, subtracting intgral
      ; through y-boundary guard cells
      f[xi,yi] = f[xi,yi] - f[xi,yi[mesh.y_boundary_guards]]
    ENDIF

    IF NOT KEYWORD_SET(nosmooth) THEN BEGIN
      f[xi,yi] = SMOOTH(SMOOTH(f[xi,yi], 5, /edge_truncate), 5, /edge_truncate)
    ENDIF
    
    IF period THEN BEGIN
      ;; Only set loop integral in closed (periodic) domains i.e. the
      ;; core. Otherwise it may be overwritten by values in a PF region
       
      loop[xi] = f[xi,yi[N_ELEMENTS(yi)-1]] - f[xi,yi[0]]
    ENDIF
  ENDREP UNTIL last
  
  RETURN, f
END

; derivative function that takes in an axis as well
function dfdy, f, y, mesh

	s = size(f, /dim)
	nx = s[0]
	ny = s[1]
	result = dblarr(nx,ny)

	status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
	REPEAT BEGIN
		yi = gen_surface_hypnotoad(last=last, xi=xi)
		result[xi,yi] = DERIV(y[xi,yi],f[xi,yi])
	ENDREP UNTIL last
	return, result
end

; derivative of function in y, separating into flux regions
function dfdy_seps, f, y, mesh

	s = size(f, /dim)
	nx = s[0]
	ny = s[1]
	result = dblarr(nx,ny)
	N_ints = n_elements(mesh.npol)

	for j=0,nx-1 do begin
		ylow = 0
		for i=0,N_ints-1 do begin
			ylocs = indgen(mesh.npol[i] + mesh.n_y_boundary_guards[i])+ylow
			result[j,ylocs] = DERIV(y[j,ylocs],f[j,ylocs])
			ylow += mesh.npol[i] + mesh.n_y_boundary_guards[i]
		endfor
	endfor
	return, result

end

function dfdx, f, x
	
	nx = size(f,/dimensions)
	ny = nx[1]
	nx = nx[0]

	result = dblarr(nx,ny)
	result[*,*] = 0.0D
	for i=0, ny-1 do begin
		result[*,i] = DERIV(reform(x[*,i]),reform(f[*,i]));
	endfor
	
	return, result
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Solve for pressure and f=R*Bt using force balance
; Using CURVEFIT routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Calculate jxb force at every point, given f on flux-surfaces
; X = X  - dummy, not used
; profiles = [mu0*p(surface), f(surface)]
; force = force(nx, ny), reformed into single vector
PRO jxb_funct, X, profiles, force, pder
  COMMON jxb_com, nx, ny, indxy, psixy, Rxy, hthe, axy

  nsurf = N_ELEMENTS(profiles) / 2
  
  dpdx = profiles[0:(nsurf-1)]
  f = profiles[nsurf:*]

  ; Put profiles into 2D arrays
  fxy = DBLARR(nx, ny)
  dpxy = DBLARR(nx, ny)
  
  FOR x=0, nx-1 DO BEGIN
     FOR y=0, ny-1 DO BEGIN
        i = indxy[x,y]
        fxy[x,y] = f[i]
        dpxy[x,y] = dpdx[i]
     ENDFOR
  ENDFOR

  ; Components of B
  Btxy = fxy / Rxy

  force = Btxy*hthe*DDX(psixy, Btxy) $
    + fxy^2*axy $
    + hthe*dpxy
  
  pder = DBLARR(nx, ny, 2*nsurf)

  ; Diagonal dependencies (dF / dfi)
  dFdfi = (hthe/Rxy)*DDX(psixy, Btxy) $
         + 2.D*fxy*axy
  
  ; Set the elements of pder
  FOR x=0, nx-1 DO BEGIN
    FOR y=0, ny-1 DO BEGIN
      ; Get indices into profiles
      xp = x+1 < nx-1
      xm = x-1 > 0
      i = indxy[x,y]
      ip = indxy[xp,y]
      im = indxy[xm,y]
      
      ; f components
      pder[x,y, nsurf+i]  = dFdfi[x,y]
      dx = psixy[xp,y] - psixy[xm,y]
      pder[x,y, nsurf+ip] = pder[x,y, nsurf+ip] + hthe[x,y]*fxy[x,y]/ (Rxy[x,y]*Rxy[xp,y]*dx)
      pder[x,y, nsurf+im] = pder[x,y, nsurf+im] - hthe[x,y]*fxy[x,y]/ (Rxy[x,y]*Rxy[xm,y]*dx)
      
      ; p component
      pder[x,y, i] = hthe[x,y]
    ENDFOR
  ENDFOR
  
  force = REFORM(force, nx*ny)
  pder = REFORM(pder, nx*ny, 2*nsurf)
END

FUNCTION fit_profiles, mesh, psixy, Rxy, hthe, Bpxy, Btxy, dpdx
  COMMON jxb_com, nx, ny, indxy, psi, R, h, axy
  
  MU = 4.d-7*!DPI
  
  psi = psixy
  r = Rxy
  h = hthe

  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]

  ; Map between location in xy and surface number
  indxy = INTARR(nx, ny)
  
  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  i = 0
  REPEAT BEGIN
    yi = gen_surface_hypnotoad(last=last, xi=xi, period=period)
    indxy[xi,yi] = i
    
    IF i EQ 0 THEN BEGIN
      farr = [MEAN(Btxy[xi,yi]*Rxy[xi,yi])]
      parr = [MEAN(dpdx[xi,yi])]
    ENDIF ELSE BEGIN
      farr = [farr, MEAN(Btxy[xi,yi]*Rxy[xi,yi])]
      parr = [parr, MEAN(dpdx[xi,yi])]
    ENDELSE
    
    i = i + 1
  ENDREP UNTIL last
  nsurf = N_ELEMENTS(farr)
  
  profiles = [MU*parr, farr]

  ; Calculate useful quantities
  axy = hthe*DDX(psixy, Rxy)/(Rxy^3)
  
  fit = CURVEFIT(FINDGEN(nx*ny), $
                 REFORM(-Bpxy*DDX(psixy, Bpxy*hthe), nx*ny), $
                 weights, $
                 profiles, $
                 function_name="jxb_funct", /noder)
  
  Btxy2 = DBLARR(nx, ny)
  dpdx2 = DBLARR(nx, ny)
  
  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  i = 0
  REPEAT BEGIN
    yi = gen_surface_hypnotoad(last=last, xi=xi, period=period)
    Btxy2[xi, yi] = profiles[nsurf+i] / Rxy[xi,yi]
    dpdx2[xi, yi] = profiles[i]
    i = i + 1
  ENDREP UNTIL last

  STOP
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Correct hthe using force balance

FUNCTION new_hfunc, h
  COMMON hvars, psi, fixpos, h0, a, b

  IF fixpos EQ 0 THEN BEGIN
    h2 = [h0, h]
  ENDIF ELSE IF fixpos EQ N_ELEMENTS(psi)-1 THEN BEGIN
    h2 = [h, h0]
  ENDIF ELSE BEGIN
    h2 = [h[0:(fixpos-1)], h0, h[fixpos:*]]
  ENDELSE

  f = a*h2 + b*DERIV(psi, h2)
  
  IF fixpos EQ 0 THEN BEGIN
    f = f[1:*]
  ENDIF ELSE IF fixpos EQ N_ELEMENTS(psi)-1 THEN BEGIN
    f = f[0:(N_ELEMENTS(f)-2)]
  ENDIF ELSE BEGIN
    f = [f[0:(fixpos-1)], f[(fixpos+1):*]]
  ENDELSE

  RETURN, f
END

FUNCTION correct_hthe, Rxy, psixy, Btxy, Bpxy, hthe, pressure, fixhthe=fixhthe
  COMMON hvars, xarr, fixpos, h0, a, b
  
  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]

  MU = 4.d-7*!DPI

  IF NOT KEYWORD_SET(fixhthe) THEN fixhthe = 0
  IF fixhthe LT 0 THEN fixhthe = 0
  IF fixhthe GT nx-1 THEN fixhthe = nx-1

  fixpos = fixhthe
  PRINT, "FIX = ", fixhthe
  
  axy = Btxy*DDX(psixy, Btxy) + Bpxy*DDX(psixy, Bpxy) $
    + Btxy^2*DDX(psixy, Rxy)/Rxy + MU*DDX(psixy, pressure)
  bxy = Bpxy^2

  nh = DBLARR(nx, ny)
  nh[fixhthe,*] = hthe[fixhthe,*]
  FOR i=0, ny-1 DO BEGIN
    PRINT, "Correcting y index ", i
    xarr = psixy[*,i]
    a = axy[*,i]
    b = bxy[*,i]
    h0 = hthe[fixhthe,i]
    
    IF fixhthe EQ 0 THEN BEGIN
      htmp = REFORM(hthe[1:*,i])
    ENDIF ELSE IF fixhthe GE nx-1 THEN BEGIN
      ; fix last point
      htmp = hthe[0:(nx-2),i]
      fixhthe = nx-1
    ENDIF ELSE BEGIN
      ; fix somewhere in the middle  
      htmp = [hthe[0:(fixhthe-1),i], hthe[(fixhthe+1):*,i]]
    ENDELSE
    
    htmp = NEWTON(htmp, "new_hfunc")
    
    IF fixhthe EQ 0 THEN BEGIN
      nh[1:*] = htmp
    ENDIF ELSE IF fixhthe GE nx-1 THEN BEGIN
      nh[0:(nx-2), i] = htmp
    ENDIF ELSE BEGIN
      nh[0:(fixhthe-1), i] = htmp[0:(fixhthe-1)]
      nh[(fixhthe+1):*, i]  = htmp[fixhthe:*]
    ENDELSE
    
    w = WHERE(nh[*,i] LT 0.0D, count)
    IF count GT 0 THEN BEGIN
      PRINT, "Error in hthe solver: Negative solution at y = ", i
      ;STOP
    ENDIF
  ENDFOR

  RETURN, nh
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Refine an equilibrium
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION grid_newt, data

  COMMON grid_newt_com, nx, ny, psixy, gs_f, dpdpsi, R0, Z0, min_f, xfix
  
  n = nx*ny
  dxin = REFORM(data, nx-1, ny)
  dx = DBLARR(nx, ny)
  IF xfix LE 0 THEN BEGIN
     dx[1:*,*] = dxin
  ENDIF ELSE IF xfix GE (nx-1) THEN BEGIN
     dx[0:(nx-2),*] = dxin
  ENDIF ELSE BEGIN
     dx[0:(xfix-1),*] = dxin[0:(nfix-1),*]
     dx[(xfix+1):*,*] = dxin[nfix:*,*]
  ENDELSE

  xpos = dx
  FOR i=0, nx-1 DO xpos[i,*] = xpos[i,*] + i

  Rxy = DBLARR(nx, ny)
  Zxy = Rxy
  
  FOR y=0, ny-1 DO BEGIN
     Rxy[*,y] = INTERPOL(R0[*,y], FINDGEN(nx), xpos[*,y], /spline)
     Zxy[*,y] = INTERPOL(Z0[*,y], FINDGEN(nx), xpos[*,y], /spline)
  ENDFOR

  ; calculate Bpxy, Btxy and hthe

  Btxy = DBLARR(nx, ny)
  FOR x=0, nx-1 DO Btxy[x,*] = gs_f[x] / Rxy[x,*]
  hthe = calc_hthe(Rxy, Zxy)
  Bpxy = calc_bp(psixy, Rxy, Zxy)
  
  F = -1.0D*calc_force(psixy, Bpxy, Btxy, hthe, Rxy, dpdpsi)

  fm = MAX(ABS(F))

  IF (fm LT min_f) OR (min_f LT 0.0D) THEN BEGIN
      min_f = fm
      PRINT, MAX(ABS(Rxy - R0)), MAX(ABS(Zxy - Z0)), MAX(ABS(F))
  ENDIF

  ;!P.multi=[0,0,2,0,0]
  ;surface, Bpxy, chars=3
  ;surface, F, chars=3

  IF yind LT 0 THEN val = REFORM(F, n) ELSE val = F[*,yind]
  
  RETURN, val
END

FUNCTION refine_equlibrium, mesh
  
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Main grid processing routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO process_grid, rz_grid, mesh, output=output, poorquality=poorquality, $
                  gui=gui, parent=parent, reverse_bt=reverse_bt, $
                  curv=curv, smoothpressure=smoothpressure, $
                  smoothhthe=smoothhthe, smoothcurv=smoothcurv, $
                  settings=settings
  
  IF NOT KEYWORD_SET(settings) THEN BEGIN
    ; Create an empty structure
    settings = {dummy:0}
  ENDIF
  ; Check settings
  str_check_present, settings, 'calcp', -1
  str_check_present, settings, 'calcbt', -1
  str_check_present, settings, 'calchthe', -1
  str_check_present, settings, 'calcjpar', -1
  str_check_present, settings, 'orthogonal_coordinates_output', -1
  ; settings.y_boundary_guards is required, so don't set default value
  
  ;CATCH, err
  ;IF err NE 0 THEN BEGIN
  ;  PRINT, "PROCESS_GRID failed"
  ;  PRINT, "   Error message: "+!ERROR_STATE.MSG
  ;  CATCH, /cancel
  ;  RETURN
  ;ENDIF

  MU = 4.d-7*!DPI

  poorquality = 0

  IF NOT KEYWORD_SET(output) THEN output="bout.grd.nc"
  
  ; Size of the mesh
  nx = FIX(TOTAL(mesh.nrad))
  ny = FIX(TOTAL(mesh.npol))
  ny_total = FIX(TOTAL(mesh.npol+mesh.n_y_boundary_guards)) ; including y-boundary cells

  ; Find the midplane
  ymid = 0
  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface_hypnotoad(period=period, last=last, xi=xi)
    IF period THEN BEGIN
      rm = MAX(mesh.Rxy[xi,yi], ymid)
      ymid = yi[ymid]
      BREAK
    ENDIF
  ENDREP UNTIL last
  

  Rxy = mesh.Rxy
  Zxy = mesh.Zxy

  ; Note: The mesh psi normalisation uses the separatrix as psi_n = 1
  ; but the input EQDSK file may have a different normalisation
  psixy = mesh.psixy*mesh.fnorm + mesh.faxis                              ; Non-normalised psi
  psixy_eq = (psixy - rz_grid.simagx) / (rz_grid.sibdry - rz_grid.simagx) ; Normalised using EQDSK file conventions
  
  pressure = DBLARR(nx, ny_total)
  
  ; Use splines to interpolate pressure profile
  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  REPEAT BEGIN
    ; Get the next domain
    yi = gen_surface_hypnotoad(period=period, last=last, xi=xi)
    IF period AND (psixy_eq[xi,yi[0]] GE 0) AND (psixy_eq[xi,yi[0]] LE 1) THEN BEGIN
      ; Pressure only given on core surfaces
      ; Since psi normalised differently, it might go out of range 
      pressure[xi,yi] = SPLINE(rz_grid.npsigrid, rz_grid.pres, psixy_eq[xi,yi[0]], /double)
    ENDIF ELSE BEGIN
      pressure[xi,yi] = rz_grid.pres[N_ELEMENTS(rz_grid.pres)-1]
    ENDELSE
  ENDREP UNTIL last
  
  ; Add a minimum amount
  IF MIN(pressure) LT 1.0d-2*MAX(pressure) THEN BEGIN
    PRINT, "****Minimum pressure is very small:", MIN(pressure)
    PRINT, "****Setting minimum pressure to 1% of maximum"
    pressure = pressure + 1e-2*MAX(pressure)
  ENDIF
  
  IF KEYWORD_SET(smoothpressure) THEN BEGIN
    p0 = pressure[*,ymid] ; Keep initial pressure for comparison
    REPEAT BEGIN
      !P.multi=[0,0,2,0,0]
      PLOT, p0, xtitle="X index", ytitle="pressure at y="+STRTRIM(STRING(ymid),2)+" dashed=original", color=1, lines=1
      OPLOT, pressure[*,ymid], color=1
      PLOT, DERIV(p0), xtitle="X index", ytitle="DERIV(pressure)", color=1, lines=1
      OPLOT, DERIV(pressure[*,ymid]), color=1
      sm = get_yesno("Smooth pressure profile?", gui=gui, dialog_parent=parent)
      IF sm THEN BEGIN
        ; Smooth the pressure profile
        
        p2 = pressure
        FOR i=0, 5 DO BEGIN
          status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
          REPEAT BEGIN
            ; Get the next domain
            yi = gen_surface_hypnotoad(period=period, last=last, xi=xi)
            
            IF (xi GT 0) AND (xi LT (nx-1)) THEN BEGIN
              FOR j=0,N_ELEMENTS(yi)-1 DO BEGIN
                p2[xi,yi[j]] = 0.5D*pressure[xi,yi[j]] + $
                  0.25D*(pressure[xi-1,yi[j]] + pressure[xi+1,yi[j]])
              ENDFOR
            ENDIF
            
            ; Make sure it's still constant on flux surfaces
            p2[xi,yi] = MEAN(p2[xi,yi])
          ENDREP UNTIL last
          pressure = p2
        ENDFOR
      ENDIF
    ENDREP UNTIL sm EQ 0
  ENDIF

  IF MIN(pressure) LT 0.0D THEN BEGIN
    PRINT, ""
    PRINT, "============= WARNING =============="
    PRINT, "Poor quality equilibrium: Pressure is negative"
    PRINT, ""
    poorquality = 1
  ENDIF
  
  dpdpsi = DDX(psixy, pressure)

  ;IF MAX(dpdpsi)*mesh.fnorm GT 0.0D THEN BEGIN
  ;  PRINT, ""
  ;  PRINT, "============= WARNING =============="
  ;  PRINT, "Poor quality equilibrium: Pressure is increasing radially"
  ;  PRINT, ""
  ;  poorquality = 1
  ;ENDIF

  ; Grid spacing
  dx = DBLARR(nx, ny_total)
  FOR y=0, ny_total-1 DO BEGIN
    dx[0:(nx-2),y] = psixy[1:*,y] - psixy[0:(nx-2),y]
    dx[nx-1,y] = dx[nx-2,y]
  ENDFOR
  
  ; Sign
  bpsign = 1.D
  xcoord = psixy
  IF MIN(dx) LT 0.D THEN BEGIN
    bpsign = -1.D
    dx = -dx ; dx always positive
    xcoord = -xcoord
  ENDIF

  dtheta = 2.D*!DPI / DOUBLE(ny)
  dy = DBLARR(nx, ny_total) + dtheta
  
  ; B field components
  ; Following signs mean that psi increasing outwards from
  ; core to edge results in Bp clockwise in the poloidal plane
  ; i.e. in the positive Grad Theta direction.
  
  Brxy = mesh.dpsidZ / Rxy
  Bzxy = -mesh.dpsidR / Rxy
  Bpxy = SQRT(Brxy^2 + Bzxy^2)
  ; Determine direction (dot B with grad y vector)
  
  dot = Brxy[0,ymid]*(Rxy[0,ymid+1] - Rxy[0,ymid-1]) + $
    Bzxy[0,ymid]*(Zxy[0,ymid+1] - Zxy[0,ymid-1])
  
  
  IF dot LT 0.D THEN BEGIN
    PRINT, "**** Poloidal field is in opposite direction to Grad Theta -> Bp negative"
    Bpxy = -Bpxy
    IF bpsign GT 0 THEN STOP ; Should be negative
    bpsign = -1.0D
  ENDIF ELSE BEGIN
    IF bpsign LT 0 THEN STOP ; Should be positive
    bpsign = 1.D
  ENDELSE

  ; Get toroidal field from poloidal current function fpol
  Btxy = DBLARR(nx, ny_total)
  fprime = Btxy
  fp = DERIV(rz_grid.npsigrid*(rz_grid.sibdry - rz_grid.simagx), rz_grid.fpol)
  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  REPEAT BEGIN
    ; Get the next domain
    yi = gen_surface_hypnotoad(period=period, last=last, xi=xi)

    IF period AND (psixy_eq[xi,yi[0]] GE 0) AND (psixy_eq[xi,yi[0]] LE 1) THEN BEGIN
      ; In the core
      ;fpol = INTERPOL(rz_grid.fpol, rz_grid.npsigrid, mesh.psixy[xi,yi], /spline)
      fpol = SPLINE(rz_grid.npsigrid, rz_grid.fpol, psixy_eq[xi,yi[0]], /double)
      fprime[xi,yi] = SPLINE(rz_grid.npsigrid, fp, psixy_eq[xi,yi[0]], /double)
    ENDIF ELSE BEGIN
      ; Outside core. Could be PF or SOL
      fpol = rz_grid.fpol[N_ELEMENTS(rz_grid.fpol)-1]
      fprime[xi,yi] = 0.D
    ENDELSE
    Btxy[xi,yi] = fpol / Rxy[xi,yi]
  ENDREP UNTIL last
  
  ; Total B field
  Bxy = SQRT(Btxy^2 + Bpxy^2)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Go through the domains to get a starting estimate
  ; of hthe
  hthe = DBLARR(nx, ny_total)

  ; Pick a midplane index
  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  REPEAT BEGIN
    ; Get the next domain
    yi = gen_surface_hypnotoad(period=period, last=last, xi=xi)
    
    IF period THEN BEGIN
      ; In the core
      rmax = MAX(Rxy[xi,yi], ymid)
      ymidplane = yi[ymid]
      BREAK
    ENDIF
  ENDREP UNTIL last

  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  REPEAT BEGIN
    ; Get the next domain
    yi = gen_surface_hypnotoad(period=period, last=last, xi=xi)
    
    n = N_ELEMENTS(yi)
    
    ; Get distance along this line
    
    IF period THEN BEGIN
      ; Periodic, so can use FFT
      ;drdi = REAL_PART(fft_deriv(Rxy[xi, yi]))
      ;dzdi = REAL_PART(fft_deriv(Zxy[xi, yi]))
      drdi = (DERIV([REFORM(Rxy[xi,yi[n-2:*]]), $
                     REFORM(Rxy[xi, yi]), $
                     REFORM(Rxy[xi,yi[0:1]])]))[2:(n+1)]
      dzdi = (DERIV([REFORM(Zxy[xi,yi[n-2:*]]), $
                     REFORM(Zxy[xi, yi]), $
                     REFORM(Zxy[xi,yi[0:1]])]))[2:(n+1)]
    ENDIF ELSE BEGIN
      ; Non-periodic
      drdi = DERIV(Rxy[xi, yi])
      dzdi = DERIV(Zxy[xi, yi])
    ENDELSE
    
    dldi = REFORM(SQRT(drdi^2 + dzdi^2))
    
    IF 0 THEN BEGIN

      ; Need to smooth to get sensible results
      IF period THEN BEGIN
        n = N_ELEMENTS(dldi)
        dldi = (SMOOTH([ dldi[(n-2):*], dldi, dldi[0:1] ], 5))[2:(n+1)]
        dldi = (SMOOTH([ dldi[(n-2):*], dldi, dldi[0:1] ], 5))[2:(n+1)]
        dldi = (SMOOTH([ dldi[(n-2):*], dldi, dldi[0:1] ], 5))[2:(n+1)]
      ENDIF ELSE BEGIN
        dldi = SMOOTH(dldi, 5)
        dldi = SMOOTH(dldi, 5)
        dldi = SMOOTH(dldi, 5)
      ENDELSE
    ENDIF
    
    hthe[xi, yi] = dldi / dtheta ; First estimate of hthe
    
    ; Get outboard midplane
    IF period AND xi EQ 0 THEN BEGIN
      m = MAX(Rxy[0,yi], ymidplane)
      ymidplane = yi[ymidplane]
    ENDIF
  ENDREP UNTIL last

  PRINT, "Midplane index ", ymidplane

  fb0 = force_balance(psixy, Rxy, Bpxy, Btxy, hthe, pressure)
  PRINT, "Force imbalance: ", MEAN(ABS(fb0)), MAX(ABS(fb0))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Correct pressure using hthe
  
  PRINT, "Calculating pressure profile from force balance"

  CATCH, err
  IF err THEN BEGIN
    CATCH, /cancel
    PRINT, "WARNING: Pressure profile calculation failed: ", !ERROR_STATE.MSG 
  ENDIF ELSE BEGIN

    ; Calculate force balance
    dpdx = ( -Bpxy*DDX(xcoord, Bpxy * hthe) - Btxy*hthe*DDX(xcoord, Btxy) - (Btxy*Btxy*hthe/Rxy)*DDX(xcoord, Rxy) ) / (MU*hthe)
    
    ; Surface average
    dpdx2 = surface_average(dpdx, mesh)
    
    pres = DBLARR(nx, ny_total)
    ; Integrate to get pressure
    FOR i=0, ny_total-1 DO BEGIN
      pres[*,i] = int_func(psixy[*,i], dpdx2[*,i], /simple)
      pres[*,i] = pres[*,i] - pres[nx-1,i]
    ENDFOR
    
    status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
    REPEAT BEGIN
      ; Get the next domain
      yi = gen_surface_hypnotoad(period=period, last=last, xi=xi)
      
      ma = MAX(pres[xi,yi])
      FOR i=0, N_ELEMENTS(yi)-1 DO BEGIN
        pres[*,yi[i]] = pres[*,yi[i]] - pres[xi,yi[i]] + ma
      ENDFOR
    ENDREP UNTIL last
    
    pres = pres - MIN(pres)
  
    ; Some sort of smoothing here?
  
    fb0 = force_balance(psixy, Rxy, Bpxy, Btxy, hthe, pres)
    PRINT, "Force imbalance: ", MEAN(ABS(fb0)), MAX(ABS(fb0))
  
    !P.MULTI=[0,0,2,0,0]
    SURFACE, pressure, xtitle="X", ytitle="Y", title="Input pressure", chars=2, color=1
    SURFACE, pres, xtitle="X", ytitle="Y", title="New pressure", chars=2,color=1
  
    calcp = settings.calcp
    
    IF calcp EQ -1 THEN BEGIN
      calcp = get_yesno("Keep new pressure?", gui=gui, dialog_parent=parent)
    ENDIF ELSE WAIT, 2
    IF calcp EQ 1 THEN BEGIN
      pressure = pres
      dpdpsi = dpdx2
    ENDIF
  ENDELSE
  CATCH, /cancel
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Correct f = RBt using force balance

  calcbt = settings.calcbt
  IF calcbt EQ -1 THEN calcbt = get_yesno("Correct f=RBt using force balance?", gui=gui, dialog_parent=parent)
  IF calcbt EQ 1 THEN BEGIN

    new_Btxy = newton_bt(psixy, Rxy, Btxy, Bpxy, pres, hthe, mesh)
    
    fb0 = force_balance(psixy, Rxy, Bpxy, new_Btxy, hthe, pressure)
    PRINT, "force imbalance: ", MEAN(ABS(fb0)), MAX(ABS(fb0))
    
    !P.MULTI=[0,0,2,0,0]
    SURFACE, Btxy, xtitle="X", ytitle="Y", title="Input Bt", chars=2,color=1
    SURFACE, new_Btxy, xtitle="X", ytitle="Y", title="New Bt", chars=2,color=1
    
    calcbt = settings.calcbt
    IF calcbt EQ -1 THEN calcbt = get_yesno("Keep new Bt?", gui=gui, dialog_parent=parent)
    IF calcbt EQ 1 THEN BEGIN
      Btxy = new_Btxy
      Bxy = SQRT(Btxy^2 + Bpxy^2)
    ENDIF
  ENDIF
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; CALCULATE HTHE
  ; Modify hthe to fit force balance using initial guess
  ; Does not depend on signs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  calchthe = settings.calchthe
  IF calchthe EQ -1 THEN calchthe = get_yesno("Adjust hthe using force balance?", gui=gui, dialog_parent=parent) 
  IF calchthe EQ 1 THEN BEGIN
    ; This doesn't behave well close to the x-points
    fixhthe = FIX(nx / 2)
    nh = correct_hthe(Rxy, psixy, Btxy, Bpxy, hthe, pressure, fixhthe=fixhthe)
    
    fb0 = force_balance(psixy, Rxy, Bpxy, Btxy, nh, pressure)
    PRINT, "Force imbalance: ", MEAN(ABS(fb0)), MAX(ABS(fb0))
    
    PRINT, "Maximum difference in hthe: ", MAX(ABS(hthe - nh))
    PRINT, "Maximum percentage difference: ", 100.D*MAX(ABS((hthe - nh)/hthe))

    !P.multi=[0,0,1,0,0]
    PLOT, hthe[*,0], title="Poloidal arc length at midplane. line is initial estimate", color=1
    OPLOT, nh[*,0], psym=1, color=2
    OPLOT, nh[*,0], color=2

    IF get_yesno("Keep new hthe?", gui=gui, dialog_parent=parent) THEN BEGIN
      hthe = nh
    ENDIF
  ENDIF
  
  IF KEYWORD_SET(smoothhthe) THEN BEGIN
    ; Smooth hthe to prevent large jumps in X or Y. This
    ; should be done by creating a better mesh in the first place
    
    ; Need to smooth in Y and X otherwise smoothing in X
    ; produces discontinuities in Y
    hold = hthe
    
    IF 1 THEN BEGIN
      ; Nonlinear smoothing. Tries to smooth only regions with large
      ; changes in gradient
      
      hthe = smooth_nl(hthe, mesh);
      
    ENDIF ELSE BEGIN
      ; Just use smooth in both directions
      
      FOR i=0, ny_total-1 DO BEGIN
        hthe[*,i] = SMOOTH(SMOOTH(hthe[*,i],10),10)
      ENDFOR
      
      status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
      REPEAT BEGIN
        ; Get the next domain
        yi = gen_surface_hypnotoad(period=period, last=last, xi=xi)
        
        n = N_ELEMENTS(yi)
        
        IF period THEN BEGIN
          hthe[xi,yi] = (SMOOTH([reform(hthe[xi,yi[(n-4):(n-1)]]), reform(hthe[xi,yi]), reform(hthe[xi,yi[0:3]])], 4))[4:(n+3)]
        ENDIF ELSE BEGIN
          hthe[xi,yi] = SMOOTH(hthe[xi,yi], 4)
        ENDELSE
      ENDREP UNTIL last
    ENDELSE
  ENDIF

  ; Need yxy values at all points for nonorthogonal calculations
  yxy = DBLARR(nx, ny_total)
  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  REPEAT BEGIN
    ; Get the next domain
    yi = gen_surface_hypnotoad(period=period, last=last, xi=xi)
    IF period EQ 0 THEN BEGIN
      ; not periodic, set yxy to zero after y-boundary guard cells
      yxy[xi,yi] = (DINDGEN(N_ELEMENTS(yi)) - settings.y_boundary_guards)*dtheta
    ENDIF ELSE BEGIN
      ; periodic, no y-boundary guard cells
      yxy[xi,yi] = (DINDGEN(N_ELEMENTS(yi)))*dtheta
    ENDELSE
  ENDREP UNTIL last

  ; Calculate hrad and dhrad for thetaxy calculation
  hrad = dblarr(nx,ny_total) 
  dhrad = dblarr(nx,ny_total) 
  for j=0,ny_total-1 do begin
    for i=0, nx-1 do begin
      if(i eq 0) then begin
        r2 = Rxy[i,j]
        z2 = Zxy[i,j]
        r3 = Rxy[i+1,j]
        z3 = Zxy[i+1,j]
	hrad[i,j] = sqrt((r2-r3)^2+(z2-z3)^2)
	dhrad[i,j] = sqrt((r2-r3)^2+(z2-z3)^2)
      endif else if(i eq nx-1) then begin
        r1 = Rxy[i-1,j]
        z1 = Zxy[i-1,j]
        r2 = Rxy[i,j]
        z2 = Zxy[i,j]
	hrad[i,j] = hrad[i-1,j] + sqrt((r2-r1)^2+(z2-z1)^2)
	dhrad[i,j] = sqrt((r2-r1)^2+(z2-z1)^2)
      endif else begin
        r1 = Rxy[i-1,j]
        z1 = Zxy[i-1,j]
        r2 = Rxy[i,j]
        z2 = Zxy[i,j]
        r3 = Rxy[i+1,j]
        z3 = Zxy[i+1,j]
	hrad[i,j] = hrad[i-1,j] + sqrt(((r2-r1)/2.D + (r3-r2)/2.D)^2+((z2-z1)/2.D + (z3-z2)/2.D)^2)
	dhrad[i,j] = sqrt(((r2-r1)/2.D + (r3-r2)/2.D)^2+((z2-z1)/2.D + (z3-z2)/2.D)^2)
      endelse
    endfor
  endfor

  ; Calculate beta (angle between x and y coord)
  ; beta_method - 0 to use local grad psi for calculation
  ; 		- 1 to use gridpoint locations to calculate angle
  beta_method = 0
retrybetacalc:
  beta = calc_beta(Rxy,Zxy,mesh,rz_grid,beta_method) ; use local gradient

  ; Calculate eta (poloidal non-orthogonality parameter)
  eta = sin(beta) ; from geometry
  yshift = intx(hrad, eta, /simple) ; b/c angle was calculated real space, integrate in real space as well (hrad instead of psixy)
  thetaxy = yxy + yshift

  G = 1.D - dfdy_seps(yshift,thetaxy,mesh)
  dyshiftdy = dfdy_seps(yshift,yxy,mesh)
;   G = 1.D - dfdy(yshift,thetaxy,mesh)
;   dyshiftdy = dfdy(yshift,yxy,mesh)

  ; Calculate field-line pitch
  pitch = hthe * Btxy / (Bpxy * Rxy)
  
  ; derivative with psi
  dqdpsi = DDX(psixy, pitch)

  ; Calculate zshift (qinty), sinty = d(zshift)/dpsi, and H = d(zshift)/dtheta
  qinty = my_int_y(pitch*(1.D + dyshiftdy), yxy, mesh, /nosmooth, loop=qloop)
  sinty = DDX(psixy,qinty)
  H = dfdy_seps(qinty,thetaxy,mesh)
;   H = dfdy(qinty,thetaxy,mesh)

  ; NOTE: This is only valid in the core
  pol_angle = DBLARR(nx,ny_total)
  FOR i=0, nx-1 DO pol_angle[i, *] = 2.0D*!DPI * qinty[i,*] / qloop[i]
  
  ;;;;;;;;;;;;;;;;;;;; THETA_ZERO ;;;;;;;;;;;;;;;;;;;;;;
  ; re-set zshift to be zero at the outboard midplane
  
  PRINT, "MIDPLANE INDEX = ", ymidplane
  
  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  REPEAT BEGIN
    ; Get the next domain
    yi = gen_surface_hypnotoad(period=period, last=last, xi=xi)
    
    w = WHERE(yi EQ ymidplane, count)
    IF count GT 0 THEN BEGIN
      ; Crosses the midplane
      qinty[xi, yi] = qinty[xi, yi] - qinty[xi, ymidplane]
      sinty[xi, yi] = sinty[xi, yi] - sinty[xi, ymidplane]
;       H[xi, yi] = H[xi, yi] - H[xi, ymidplane]
    ENDIF ELSE BEGIN
      ; Doesn't include a point at the midplane
      ; Set the value at the first grid-point to zero
      qinty[xi, yi] = qinty[xi, yi] - qinty[xi,yi[settings.y_boundary_guards]]
      sinty[xi, yi] = sinty[xi, yi] - sinty[xi,yi[settings.y_boundary_guards]]
;       H[xi, yi] = H[xi, yi] - H[xi,yi[0]]
    ENDELSE
  ENDREP UNTIL last

  ; Calculate metrics - check jacobian

  orthogonal_coordinates_output = settings.orthogonal_coordinates_output
  IF orthogonal_coordinates_output EQ -1 THEN orthogonal_coordinates_output = get_yesno("Output for simulations in orthogonal coordinates using ShiftedMetric?", gui=gui, dialog_parent=parent)
  IF orthogonal_coordinates_output EQ 1 THEN BEGIN
    print,""
    print,"*******************WARNING****************************************"
    print,"Calculating metrics for ShiftedMetric style orthogonal coordinates"
    print,"******************************************************************"
    print,""
    ; for orthogonal coordinates
    I = 0.D
  ENDIF ELSE BEGIN
    ; for field-aligned coordinates
    I = sinty
  ENDELSE

  g11 = (Rxy*Bpxy)^2;
  g22 = G^2/hthe^2 + eta^2*g11;
  g33 = I^2*g11 + H^2/hthe^2 + 1.0D/Rxy^2;
  g12 = -eta*g11;
  g13 = -I*g11;
  g23 = I*eta*g11 - G*H/hthe^2;

  J = hthe / Bpxy / G

  g_11 = 1.0D/g11 + (hthe*eta/G)^2 + (Rxy*H*eta/G + I*Rxy)^2;
  g_22 = hthe^2/G^2 + Rxy^2*H^2/G^2;
  g_33 = Rxy^2;
  g_12 = hthe^2*eta/G^2 + Rxy^2*H/G*(H*eta/G + I);
  g_13 = Rxy^2*(H*eta/G+I);
  g_23 = H*Rxy^2/G;

  ; check to make sure jacobian is good
  Jcheck = 1.D / sqrt(g11*g22*g33 + 2.0D*g12*g13*g23 - g11*g23*g23 - g22*g13*g13 - g33*g12*g12);
  whr = where(abs(J-Jcheck) gt 0.01D,count)
  if(count gt 0) then begin
    if(beta_method EQ 0) then begin
	print,""
	print,"*****************************************************************"
	print,"WARNING: Jacobians not consistent - trying other beta_calc method"
	print,"*****************************************************************"
	print,""
	beta_method = 1
        goto, retrybetacalc
    endif else begin
	print,""
	print,"********************************************************************"
	print,"WARNING: Jacobians not consistent - both beta_calc methods attempted"
	print,"********************************************************************"
	print,""
    endelse
  endif
  PRINT, ""
  PRINT, "==== Calculating curvature ===="
  
  ;;;;;;;;;;;;;;;;;;;; CURVATURE ;;;;;;;;;;;;;;;;;;;;;;;
  ; Calculating b x kappa
  
  IF NOT KEYWORD_SET(curv) THEN BEGIN
    
    PRINT, "*** Calculating curvature in toroidal coordinates"
    
    
    curvature, nx, ny_total, DOUBLE(Rxy), DOUBLE(Zxy), DOUBLE(brxy), DOUBLE(bzxy), DOUBLE(btxy), $
      DOUBLE(psixy), DOUBLE(thetaxy), hthe, $
      bxcv=bxcv, mesh=mesh

    bxcvx = bpsign*bxcv.psi 
    bxcvy = bxcv.theta
    bxcvz = bpsign*(bxcv.phi - I*bxcv.psi - pitch*bxcv.theta)

    ; x borders
    bxcvx[0,*] = bxcvx[1,*]
    bxcvx[nx-1,*] = bxcvx[nx-2,*]
    
    bxcvy[0,*] = bxcvy[1,*]
    bxcvy[nx-1,*] = bxcvy[nx-2,*]
    
    bxcvz[0,*] = bxcvz[1,*]
    bxcvz[nx-1,*] = bxcvz[nx-2,*]

  ENDIF ELSE IF curv EQ 1 THEN BEGIN
    ; Calculate on R-Z mesh and then interpolate onto grid
    ; ( cylindrical coordinates)

    PRINT, "*** Calculating curvature in cylindrical coordinates"
    
    bxcv = rz_curvature(rz_grid)
    
    ; DCT methods cause spurious oscillations
    ; Linear interpolation seems to be more robust
    bxcv_psi = INTERPOLATE(bxcv.psi, mesh.Rixy, mesh.Zixy, /DOUBLE)
    bxcv_theta = INTERPOLATE(bxcv.theta, mesh.Rixy, mesh.Zixy, /DOUBLE) / hthe
    bxcv_phi = INTERPOLATE(bxcv.phi, mesh.Rixy, mesh.Zixy, /DOUBLE)
    
    ; If Bp is reversed, then Grad x = - Grad psi
    bxcvx = bpsign*bxcv_psi
    bxcvy = bxcv_theta
    bxcvz = bpsign*(bxcv_phi - I*bxcv_psi - pitch*bxcv_theta)
  ENDIF ELSE IF curv EQ 2 THEN BEGIN
    ; Curvature from Curl(b/B)
    
    bxcvx = bpsign*(Bpxy * Btxy*Rxy * DDY(1.D / Bxy, mesh) / hthe)
    bxcvy = -bpsign*Bxy*Bpxy * DDX(xcoord, Btxy*Rxy/Bxy^2) / (2.D*hthe)
    bxcvz = Bpxy^3 * DDX(xcoord, hthe/Bpxy) / (2.D*hthe*Bxy) - Btxy*Rxy*DDX(xcoord, Btxy/Rxy) / (2.D*Bxy) - I*bxcvx
    
  ENDIF ELSE BEGIN
    ; calculate in flux coordinates.
    
    PRINT, "*** Calculating curvature in flux coordinates"
    
    dpb = DBLARR(nx, ny_total)      ; quantity used for y and z components
    
    FOR i=0, ny_total-1 DO BEGIN
      dpb[*,i] = MU*dpdpsi/Bxy[*,i]
    ENDFOR
    dpb = dpb + DDX(xcoord, Bxy)

    bxcvx = bpsign*(Bpxy * Btxy*Rxy * DDY(1.D / Bxy, mesh) / hthe)
    bxcvy = bpsign*(Bpxy*Btxy*Rxy*dpb / (hthe*Bxy^2))
    bxcvz = -dpb - I*bxcvx
  ENDELSE
  

  IF KEYWORD_SET(smoothcurv) THEN BEGIN
    ; Smooth curvature to prevent large jumps
    
    ; Nonlinear smoothing. Tries to smooth only regions with large
    ; changes in gradient
    
    bz = bxcvz + I * bxcvx
    
    PRINT, "Smoothing bxcvx..."
    bxcvx = smooth_nl(bxcvx, mesh)
    PRINT, "Smoothing bxcvy..."
    bxcvy = smooth_nl(bxcvy, mesh)
    PRINT, "Smoothing bxcvz..."
    bz = smooth_nl(bz, mesh)
    
    bxcvz = bz - I * bxcvx
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; CALCULATE PARALLEL CURRENT
  ; 
  ; Three ways to calculate Jpar0:
  ; 1. From fprime and pprime
  ; 2. From Curl(B) in field-aligned coords
  ; 3. From the curvature
  ; 
  ; Provides a way to check if Btor should be reversed
  ;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  PRINT, ""
  PRINT, "==== Calculating parallel current ===="
  
  jpar0 = - Bxy * fprime / MU - Rxy*Btxy * dpdpsi / Bxy
  
  ; Set to zero in PF and SOL
  status = gen_surface_hypnotoad(mesh=mesh) ; Start generator
  REPEAT BEGIN
    ; Get the next domain
    yi = gen_surface_hypnotoad(period=period, last=last, xi=xi)
    
    IF NOT period THEN jpar0[xi,yi] = 0.0D
  ENDREP UNTIL last
  
  ; Curl(B) expression for Jpar0 (very noisy usually)
  j0 = bpsign*((Bpxy*Btxy*Rxy/(Bxy*hthe))*( DDX(xcoord, Bxy^2*hthe/Bpxy) - bpsign*Btxy*Rxy*DDX(xcoord,Btxy*hthe/(Rxy*Bpxy)) ) $
        - Bxy*DDX(xcoord, Btxy*Rxy)) / MU
  
  ; Create a temporary mesh structure to send to adjust_jpar
  tmp_mesh = CREATE_STRUCT(mesh, $
                           "bxcvx", bxcvx, "bxcvy", bxcvy, "bxcvz", bxcvz, $
                           "Bpxy", Bpxy, "Btxy", Btxy, "Bxy", Bxy, $
                           "dx", dx, "dy", dy, $
                           "hthe", hthe, "jpar0", jpar0, "pressure", pressure)
  tmp_mesh.psixy = psixy
  
  adjust_jpar, tmp_mesh, jpar=jpar, /noplot
  
  !P.multi=[0,2,2,0,0]
  surface, jpar0, xtitle="X", ytitle="Y", title="Jpar from F' and P'", chars=2, color=1
  surface, jpar, xtitle="X", ytitle="Y", title="Jpar from curvature", chars=2, color=1
  
  PLOT, jpar0[0,*], tit="jpar at x=0.D Solid from f' and p'", yr=[MIN([jpar0[0,*],jpar[0,*]]), $
                                                                 MAX([jpar0[0,*],jpar[0,*]])]
  OPLOT, jpar[0,*], psym=1
  
  PLOT, jpar0[*,ymidplane], tit="Jpar at y="+STR(ymidplane)+" Solid from f' and p'", $
    yr=[MIN([jpar0[*,ymidplane],jpar[*,ymidplane]]), $
        MAX([jpar0[*,ymidplane],jpar[*,ymidplane]])]
  OPLOT, jpar[*,ymidplane], psym=1
  
  !P.multi=0
  
  calcjpar = settings.calcjpar
  IF calcjpar EQ -1 THEN calcjpar = get_yesno("Use Jpar from curvature?", gui=gui, dialog_parent=parent)
  IF calcjpar EQ 1 THEN BEGIN
    Jpar0 = Jpar
  ENDIF
  
  IF 0 THEN BEGIN
    
    ; Try smoothing jpar0 in psi, preserving zero points and maxima
    jps = jpar0
    FOR y=0,ny_total-1 DO BEGIN
      j = jpar0[*,y]
      js = j
      ma = MAX(ABS(j), ip)
      IF (ma LT 1.d-4) OR (ip EQ 0) THEN BEGIN
        jps[*,y] = j
        CONTINUE
      ENDIF
      
      level = 1.D
      ;i0 = MAX(WHERE(ABS(j[0:ip]) LT level))
      i1 = MIN(WHERE(ABS(j[ip:*]) LT level))
      
      ;IF i0 LE 0 THEN i0 = 1
      i0 = 1
      
      IF i1 EQ -1 THEN i1 = nx-2 ELSE i1 = i1 + ip
      
      IF (ip LE i0) OR (ip GE i1) THEN CONTINUE
      
      ; Now preserve starting and end points, and peak value
      div = FIX((i1-i0)/10)+1 ; reduce number of points by this factor
      
      inds = [i0] ; first point
      FOR i=i0+div, ip-div, div DO inds = [inds, i]
      inds = [inds, ip] ; Put in the peak point
      
      ; Calculate spline interpolation of inner part
      js[0:ip] = spline_mono(inds, j[inds], INDGEN(ip+1), $
                             yp0=(j[i0] - j[i0-1]), ypn_1=0.0D)
      
      inds = [ip] ; peak point
      FOR i=ip+div, i1-div, div DO BEGIN
        inds = [inds, i]
      ENDFOR
      
      inds = [inds, i1] ; Last point
      js[ip:i1] = spline_mono(inds, j[inds], ip+INDGEN(i1-ip+1), $
                              yp0=0.0D, ypn_1=(j[i1+1]-j[i1]))
      
      jps[*,y] = js
    ENDFOR
  ENDIF
  
  ;;;;;;;;;;;;;;;;;;;; TOPOLOGY ;;;;;;;;;;;;;;;;;;;;;;;
  ; Calculate indices for backwards-compatibility
  
  nr = N_ELEMENTS(mesh.nrad)
  np = N_ELEMENTS(mesh.npol)
  IF (nr EQ 2) AND (np EQ 3) THEN BEGIN
    PRINT, "Single null equilibrium"
    
    ixseps1 = mesh.nrad[0]
    ixseps2 = nx
    
    jyseps1_1 = mesh.npol[0]-1
    jyseps1_2 = mesh.npol[0] + FIX(mesh.npol[1]/2)
    ny_inner = jyseps1_2
    jyseps2_1 = jyseps1_2
    jyseps2_2 = ny - mesh.npol[2]-1

  ENDIF ELSE IF (nr EQ 3) AND (np EQ 6) THEN BEGIN
    PRINT, "Double null equilibrium"
    
    ixseps1 = mesh.nrad[0]
    ixseps2 = ixseps1 + mesh.nrad[1]
    
    jyseps1_1 = mesh.npol[0]-1
    jyseps2_1 = jyseps1_1 + mesh.npol[1]
    
    ny_inner = jyseps2_1 + mesh.npol[2] + 1
    
    jyseps1_2 = ny_inner + mesh.npol[3] - 1
    jyseps2_2 = jyseps1_2 + mesh.npol[4]
    
  ENDIF ELSE IF (nr EQ 1) AND (np EQ 1) THEN BEGIN
    
    PRINT, "Single domain"
    
    ixseps1 = nx
    ixseps2 = nx
    
    jyseps1_1 = -1
    jyseps1_2 = FIX(ny/2)
    jyseps2_1 = FIX(ny/2)
    ny_inner = FIX(ny/2)
    jyseps2_2 = ny - 1
    
  ENDIF ELSE BEGIN
    PRINT, "***************************************" 
    PRINT, "* WARNING: Equilibrium not recognised *"
    PRINT, "*                                     *"
    PRINT, "*  Check mesh carefully!              *"
    PRINT, "*                                     *"
    PRINT, "*  Contact Ben Dudson                 *"
    PRINT, "*      benjamin.dudson@york.ac.uk     *"
    PRINT, "***************************************" 
    ixseps1 = -1
    ixseps2 = -1
    
    jyseps1_1 = -1
    jyseps1_2 = FIX(ny/2)
    jyseps2_1 = FIX(ny/2)
    ny_inner = FIX(ny/2)
    jyseps2_2 = ny - 1
  ENDELSE

  PRINT, "Generating plasma profiles:"
          
  PRINT, "  1. Flat temperature profile"
  PRINT, "  2. Flat density profile"
  PRINT, "  3. Te proportional to density"
  REPEAT BEGIN
    opt = get_integer("Profile option:")
  ENDREP UNTIL (opt GE 1) AND (opt LE 3)
  
  IF opt EQ 1 THEN BEGIN
    ; flat temperature profile
    
    PRINT, "Setting flat temperature profile"
    REPEAT BEGIN
      Te_x = get_float("Temperature (eV):")
      
      
      ; get density
      Ni = pressure / (2.D*Te_x* 1.602d-19*1.0d20)
      
      PRINT, "Maximum density (10^20 m^-3):", MAX(Ni)
      
      done = get_yesno("Is this ok?")
    ENDREP UNTIL done EQ 1
    
    Te = DBLARR(nx, ny_total)+Te_x
    Ti = Te
    Ni_x = MAX(Ni)
    Ti_x = Te_x
  ENDIF ELSE IF opt EQ 2 THEN BEGIN
    PRINT, "Setting flat density profile"
    
    REPEAT BEGIN
      ni_x = get_float("Density [10^20 m^-3]:")
      
      ; get temperature
      Te = pressure / (2.D*ni_x* 1.602d-19*1.0d20)
      
      PRINT, "Maximum temperature (eV):", MAX(Te)
    ENDREP UNTIL get_yesno("Is this ok?") EQ 1
    
    Ti = Te
    Ni = DBLARR(nx, ny_total) + ni_x
    Te_x = MAX(Te)
    Ti_x = Te_x
  ENDIF ELSE BEGIN
    PRINT, "Setting te proportional to density"
    
    REPEAT BEGIN
      te_x = get_float("Maximum temperature [eV]:")
      ni_x = max(pressure) / (2.D*Te_x* 1.602d-19*1.0d20)
      
      PRINT, "Maximum density [10^20 m^-3]:", ni_x

      shape = sqrt(pressure / max(pressure))
      Te = te_x * shape
      Ni = ni_x * shape
    ENDREP UNTIL get_yesno("Is this ok?") EQ 1
    Ti = Te
    Ti_x = Te_x
  ENDELSE
  
  ; excluding y-boundary guard cells
  rmag = MAX(ABS([[Rxy[*,settings.y_boundary_guards:ny_inner+settings.y_boundary_guards-1]], [Rxy[*,ny_inner+3*settings.y_boundary_guards:-settings.y_boundary_guards-1]]]))
  PRINT, "Setting rmag = ", rmag
  
  ; excluding y-boundary guard cells
  bmag = MAX(ABS([[Bxy[*,settings.y_boundary_guards:ny_inner+settings.y_boundary_guards-1]], [Bxy[*,ny_inner+3*settings.y_boundary_guards:-settings.y_boundary_guards-1]]]))
  PRINT, "Setting bmag = ", bmag

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; save to file

  PRINT, "Writing grid to file "+output

  handle = file_open(output, /CREATE)

  ; Size of the grid

  s = file_write(handle, "nx", nx)
  s = file_write(handle, "ny", ny)
  s = file_write(handle, "y_boundary_guards", settings.y_boundary_guards)
  IF settings.y_boundary_guards GT 0 THEN BEGIN
    ; set number of y-guard cells, instead of allowing this to be set in
    ; BOUT.inp or by default. Ensures only compatibile grids with
    ; y_boundary_guards=0 or y_boundary_guards=MYG can be loaded by versions of
    ; BOUT++ older than v4.3. Note double-null grids with y_boundary_guards>0
    ; cannot be read by versions of BOUT++ earlier than v4.3.
    s = file_write(handle, "MYG", settings.y_boundary_guards)
  ENDIF

  ; Topology for original scheme
  s = file_write(handle, "ixseps1", ixseps1)
  s = file_write(handle, "ixseps2", ixseps2)
  s = file_write(handle, "jyseps1_1", jyseps1_1)
  s = file_write(handle, "jyseps1_2", jyseps1_2)
  s = file_write(handle, "jyseps2_1", jyseps2_1)
  s = file_write(handle, "jyseps2_2", jyseps2_2)
  s = file_write(handle, "ny_inner", ny_inner);
  
  ; Grid spacing
  
  s = file_write(handle, "dx", dx)
  s = file_write(handle, "dy", dy)
  
  s = file_write(handle, "ShiftAngle", qloop)
  s = file_write(handle, "zShift", qinty)
  s = file_write(handle, "pol_angle", pol_angle)
  s = file_write(handle, "ShiftTorsion", dqdpsi)
  s = file_write(handle, "Gxy", G)
  s = file_write(handle, "Hxy", H)
  s = file_write(handle, "etaxy", eta)
  s = file_write(handle, "beta", beta)
  s = file_write(handle, "yshift", yshift)

  s = file_write(handle, "Rxy",  Rxy)
  s = file_write(handle, "Zxy",  Zxy)
  s = file_write(handle, "Bpxy", Bpxy)
  s = file_write(handle, "Btxy", Btxy)
  s = file_write(handle, "Bxy",  Bxy)
  s = file_write(handle, "hthe", hthe)
  s = file_write(handle, "hrad", hrad)
  s = file_write(handle, "sinty", sinty)
  s = file_write(handle, "psixy", psixy)
  s = file_write(handle, "thetaxy", thetaxy)
  s = file_write(handle, "yxy", yxy)

  ; Topology for general configurations
  s = file_write(handle, "yup_xsplit", mesh.yup_xsplit)
  s = file_write(handle, "ydown_xsplit", mesh.ydown_xsplit)
  s = file_write(handle, "yup_xin", mesh.yup_xin)
  s = file_write(handle, "yup_xout", mesh.yup_xout)
  s = file_write(handle, "ydown_xin", mesh.ydown_xin)
  s = file_write(handle, "ydown_xout", mesh.ydown_xout)
  s = file_write(handle, "nrad", mesh.nrad)
  s = file_write(handle, "npol", mesh.npol)

  ; plasma profiles

  s = file_write(handle, "pressure", pressure)
  s = file_write(handle, "Jpar0", Jpar0)
  s = file_write(handle, "Ni0", Ni)
  s = file_write(handle, "Te0", Te)
  s = file_write(handle, "Ti0", Ti)
  s = file_write(handle, "Ni_x", Ni_x)
  s = file_write(handle, "Te_x", Te_x)
  s = file_write(handle, "Ti_x", Ti_x)
  s = file_write(handle, "bmag", bmag)
  s = file_write(handle, "rmag", rmag)

  ; Curvature
  s = file_write(handle, "bxcvx", bxcvx)
  s = file_write(handle, "bxcvy", bxcvy)
  s = file_write(handle, "bxcvz", bxcvz)

  ; type of coordinate system used to calculate metric tensor terms
  IF orthogonal_coordinates_output EQ 0 THEN BEGIN
    parallel_transform = "identity"
  ENDIF ELSE IF orthogonal_coordinates_output EQ 1 THEN BEGIN
    parallel_transform = "shiftedmetric"
  ENDIF ELSE BEGIN
    PRINT, "ERROR: Unrecognized orthogonal_coordinates_output value", $
           orthogonal_coordinates_output
  ENDELSE
  s = file_write_string(handle, "parallel_transform", parallel_transform)

  ; Metric tensor terms
  s = file_write(handle, "g11", g11)
  s = file_write(handle, "g22", g22)
  s = file_write(handle, "g33", g33)
  s = file_write(handle, "g12", g12)
  s = file_write(handle, "g13", g13)
  s = file_write(handle, "g23", g23)
  s = file_write(handle, "g_11", g_11)
  s = file_write(handle, "g_22", g_22)
  s = file_write(handle, "g_33", g_33)
  s = file_write(handle, "g_12", g_12)
  s = file_write(handle, "g_13", g_13)
  s = file_write(handle, "g_23", g_23)
  s = file_write(handle, "J", J)

  ; Psi range
  s = file_write(handle, "psi_axis", mesh.faxis)
  psi_bndry = mesh.faxis + mesh.fnorm
  s = file_write(handle, "psi_bndry", psi_bndry)

  ; save some version information
  ;
  ; BOUT++ version information: this is set when BOUT++ is configured.
  ; Hypnotoad doesn't require BOUT++ to have been configured, and IDL code may
  ; have changed since 'configure' was run, so this is not 100% reliable, but
  ; still worth saving as a sanity check
  hypnotoad_info = ROUTINE_INFO('hypnotoad', /SOURCE)
  hypnotoad_path = FILE_DIRNAME(hypnotoad_info.path)

  ; BOUT++ git hash
  SPAWN, STRJOIN(['cd ',hypnotoad_path, '&& git describe --always --abbrev=0 --dirty --match "NOT A TAG"']), bout_git_hash, EXIT_STATUS=status
  IF status THEN BEGIN
    ; hyponotoad_path is not in a git repository.
    ; BOUT++ may have been downloaded as a .tar: try to get git hash from
    ; bout-config at location relative to Hypnotoad.
    SPAWN, STRJOIN([hypnotoad_path, PATH_SEP(), '..', PATH_SEP(), '..', PATH_SEP(), '..', PATH_SEP(), 'bin/bout-config --git']), bout_git_hash, EXIT_STATUS=status
    IF status THEN BEGIN
      ; bout-config not found at relative path.
      ; BOUT++ may have been installed as a library, then bout-config should be in the $PATH
      SPAWN, 'bout-config --git', bout_git_hash, EXIT_STATUS=status
      IF status THEN BEGIN
        PRINT, '---------------------------------------------------------------'
        PRINT, 'WARNING: could not find git hash of BOUT++, not saving.'
        PRINT, 'If your BOUT++ is a git repository, something has gone wrong'
        PRINT, 'with "git describe".'
        PRINT, 'If your BOUT++ was extracted from a .tar, the bin/bout-config'
        PRINT, 'executable has failed.'
        PRINT, 'If you have installed BOUT++ as a library, "bout-config" should'
        PRINT, 'be in your $PATH, but cannot be found.'
        PRINT, '---------------------------------------------------------------'
        bout_git_hash = 0
      ENDIF
    ENDIF
  ENDIF
  IF bout_git_hash THEN BEGIN
    ; bout_git_hash as returned from SPAWN seems to have some funny character
    ; in, maybe a trailing newline. This character causes an error when trying
    ; to write as a NetCDF attribute. STRJOIN seems to fix this.
    bout_git_hash = STRJOIN(bout_git_hash, '')

    s = file_write_attribute(handle, "git_hash", bout_git_hash)
  ENDIF

  ; BOUT++ version number
  SPAWN, STRJOIN([hypnotoad_path, PATH_SEP(), '..', PATH_SEP(), '..', PATH_SEP(), '..', PATH_SEP(), 'bin/bout-config --version']), bout_version, EXIT_STATUS=status
  IF status THEN BEGIN
    ; bout-config not found at relative path.
    ; BOUT++ may have been installed as a library, then bout-config should be in the $PATH
    SPAWN, 'bout-config --version', bout_version, EXIT_STATUS=status
    IF status THEN BEGIN
      PRINT, '---------------------------------------------------------------'
      PRINT, 'WARNING: could not find version number of BOUT++, not saving.'
      PRINT, 'Could not find bin/bout-config in BOUT++ directory containing'
      PRINT, 'Hypnotoad.'
      PRINT, 'If you have installed BOUT++ as a library, "bout-config" should'
      PRINT, 'be in your $PATH, but cannot be found.'
      PRINT, '---------------------------------------------------------------'
      bout_version = 0
    ENDIF
  ENDIF
  IF bout_version THEN BEGIN
    bout_version_array = LONG(STRSPLIT(bout_version, '.', /EXTRACT))
    s = file_write_attribute(handle, "BOUT_version", bout_version_array)
  ENDIF

  ; Hypnotoad version number
  s = file_write_attribute(handle, "Hypnotoad_version", hypnotoad_version())

  file_close, handle
  PRINT, "DONE"
  
  !P.multi=[0,0,1,0,0]

  ;STOP
END
