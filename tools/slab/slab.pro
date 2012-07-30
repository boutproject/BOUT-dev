; Slab geometry grid generator
; 
; Optional keywords:
; =================
; 
; output   =   Set output file name
; thin     =   Use thin radial box approximation
;              so Bpxy = constant, but gradient is non-zero
; 
; ni  = Ion density in 10^20 m^-3 
; Ti  = temperature in eV (Te = Ti)
; 
; Rmaj    = Major radius [meters]
; rminor  = Minor radius [m]
; dr      = Radial width of box [m]
; r_wid   = Radial extent, normalised to gyro-radius r_wid = dr / rho_i
; 
; q       = Safety factor q = r*Bt/(R*Bp) at middle of box
; dq      = Change in q. Will go from q-dq/2 to q+dq/2
; 
; L_T     = Temperature length scale [m]
; eta_i   = Ratio of density to temp. length scales eta = L_n / L_T

PRO slab, output=output, thin=thin, $
          ni=ni, Ti=Ti, $
          Rmaj=Rmaj, rminor=rminor, dr=dr, $
          r_wid=r_wid, $
          q=q, dq=dq, $
          L_T=L_T, eta_i=eta_i, $
          nx=nx, ny=ny

  IF NOT KEYWORD_SET(nx) THEN nx = 68 ; Radial grid points
  IF NOT KEYWORD_SET(ny) THEN ny = 32 ; Poloidal (parallel) grid points
  
  IF NOT KEYWORD_SET(output) THEN output = "slab_"+STR(nx)+"x"+STR(ny)+".nc"

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ; Ion density in 10^20 m^-3 
  IF NOT KEYWORD_SET(ni) THEN ni = 1.
  
  ; Temperature in eV (Te = Ti)
  IF NOT KEYWORD_SET(Ti) THEN Ti = 1000
  
  ; Major radius [meters], used for setting box Z size
  IF NOT KEYWORD_SET(Rmaj) THEN Rmaj = 5

  ; Minor radius [m], Y size of domain
  IF NOT KEYWORD_SET(rminor) THEN rminor = 1
  
  ; Radial width of the domain [m]
  IF NOT KEYWORD_SET(dr) THEN dr = 0.1
  
  ; Radial extent, normalised to gyro-radius r_wid = dr / rho_i
  ; Determines magnetic field strength
  IF NOT KEYWORD_SET(r_wid) THEN r_wid = 100
  
  ; Safety factor q = r*Bt/(R*Bp) at middle of box
  IF NOT KEYWORD_SET(q) THEN q = 3

  ; Change in q. Will go from q-dq/2 to q+dq/2
  IF NOT KEYWORD_SET(dq) THEN dq = 1.
  
  ; Temperature length scale [m]
  IF NOT KEYWORD_SET(L_T) THEN L_T = 0.1
  
  ; Ratio of density to temp. length scales eta = L_n / L_T
  IF NOT KEYWORD_SET(eta_i) THEN eta_i = 1
  
  Mi = 2.*1.67262158e-27   ; Ion mass [kg]. Deuterium
  
  MU0 = 4.e-7*!PI
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  L_n = eta_i * L_T      ; Density length scale [m]
  rho_i = dr / r_wid     ; Ion Larmor radius [m]
  Bt0 = SQRT(2.*Ti*Mi / 1.602e-19) / rho_i ; Toroidal field from rho_i [T]
  Bp = rminor * Bt0 / (q * Rmaj) ; Poloidal field [T]
  
  theta = 2.*!PI * FINDGEN(ny) / FLOAT(ny)
  
  Rxy = FLTARR(nx, ny) + Rmaj
  Zxy = FLTARR(nx, ny)
  
  dy = FLTARR(nx, ny) +  2.*!PI / FLOAT(ny)
  hthe = FLTARR(nx, ny) + rminor

  Btxy = FLTARR(nx, ny) + Bt0
  
  ; Minor radius offset
  drprof = dr*((FINDGEN(nx) / FLOAT(nx-1)) - 0.5)
  ; q profile
  qprof = q + dq*((FINDGEN(nx) / FLOAT(nx-1)) - 0.5)
  ShiftAngle = qprof * 2.*!PI
  
  IF KEYWORD_SET(thin) THEN BEGIN
     ; Constant Bp, but shift angle varies
     Bpxy = FLTARR(nx, ny) + Bp
  ENDIF ELSE BEGIN
     ; Vary Bpxy to get shear
     Bpxy = FLTARR(nx, ny)
     FOR i=0, ny-1 DO Bpxy[*,i] = Bp * q / qprof
     PRINT, "Poloidal field varies from "+STR(MIN(Bpxy))+" to "+STR(MAX(Bpxy))
  ENDELSE
  
  dx = Bp * (dr / FLOAT(nx-1)) * Rxy

  Bxy = SQRT(Btxy^2 + Bpxy^2)
  
  zShift = FLTARR(nx, ny)
  FOR i=1, ny-1 DO zShift[*,i] = ShiftAngle * theta[i]/(2.*!PI)
  
  y0 = FIX(ny/2)
  zs0 = zshift[*,y0]
  FOR i=0, ny-1 DO zshift[*,i] = zshift[*,i] - zs0

  Ni0 = FLTARR(nx, ny)
  FOR i=0, ny-1 DO Ni0[*,i] = Ni * EXP(-drprof / L_n)
  Ti0 = FLTARR(nx, ny)
  FOR i=0, ny-1 DO Ti0[*,i] = Ti * EXP(-drprof / L_T)
  Te0 = Ti0
  
  pressure = Ni0 * (Ti0 + Te0) * 1.602e-19*1.0e20 ; In Pascals

  Jpar0 = FLTARR(nx, ny)
  
  IF KEYWORD_SET(thin) THEN BEGIN
    ; Bp = const across box, but dBp/dr is non-zero
    
    jphi = -(dq/dr)*Bp/q/MU0
    
    Jpar0 = (Bt0/Bxy)*jphi
  ENDIF ELSE BEGIN
    r = rminor ;+ drprof
    ;jphi = DERIV(drprof, r*Bpxy[*,0])/(r*MU0)
    jphi = -(dq/dr)*Bp*q/(qprof^2 * MU0)
    FOR y=0, ny-1 DO Jpar0[*,y] = (Bt0/Bxy[*,y]) * jphi
  ENDELSE

  ; Shape       : Rxy, Zxy
  ; Differencing: hthe, dx, dy
  ; Profiles    : Ni0, Ti0, Te0, pressure, Jpar0
  ; B field     : Btxy, Bpxy, Bxy
  ; q profile   : qprof
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Curvature

  ; Bxy is constant in x, so need to supply logB too
  
  logB = FLTARR(nx, ny)
  
  FOR y=0, ny-1 DO BEGIN
    FOR x=0,nx-1 DO BEGIN
      rpos = (FLOAT(x)/FLOAT(nx-1) - 0.5) * dr
      R = Rmaj - (rminor + rpos)*COS(theta[y])
      
      Bt = Bt0 * Rmaj / R
      
      logB[x,y] = ALOG(SQRT(Bt^2 + Bp^2))
    ENDFOR
  ENDFOR
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Topology: Just in the core
  
  ixseps1 = nx
  ixseps2 = nx
  jyseps1_1 = -1
  jyseps1_2 = FIX(ny/2)
  jyseps2_1 = jyseps1_2
  jyseps2_2 = ny-1
  ny_inner = jyseps1_2
  
  ; Only one region
  yup_xsplit   = [nx]
  ydown_xsplit = [nx]
  yup_xin = [0]
  yup_xout = [-1]
  ydown_xin = [0]
  ydown_xout = [-1]
  nrad = [nx]
  npol = [ny]
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  PRINT, "Writing grid to file "+output

  handle = file_open(output, /CREATE)
  ; Size of the grid

  s = file_write(handle, "nx", nx)
  s = file_write(handle, "ny", ny)

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
  s = file_write(handle, "dpsi", dx)
  s = file_write(handle, "dy", dy)
  
  s = file_write(handle, "ShiftAngle", ShiftAngle)
  s = file_write(handle, "zShift", zShift)
  s = file_write(handle, "qinty", zShift)
;  s = file_write(handle, "pol_angle", pol_angle)
;  s = file_write(handle, "ShiftTorsion", dqdpsi)

  s = file_write(handle, "Rxy",  Rxy)
  s = file_write(handle, "Zxy",  Zxy)
  s = file_write(handle, "Bpxy", Bpxy)
  s = file_write(handle, "Btxy", Btxy)
  s = file_write(handle, "Bxy",  Bxy)
  s = file_write(handle, "hthe", hthe)
;  s = file_write(handle, "sinty", sinty)
;  s = file_write(handle, "psixy", psixy)

  ; Topology for general configurations
  s = file_write(handle, "yup_xsplit", yup_xsplit)
  s = file_write(handle, "ydown_xsplit", ydown_xsplit)
  s = file_write(handle, "yup_xin", yup_xin)
  s = file_write(handle, "yup_xout", yup_xout)
  s = file_write(handle, "ydown_xin", ydown_xin)
  s = file_write(handle, "ydown_xout", ydown_xout)
  s = file_write(handle, "nrad", nrad)
  s = file_write(handle, "npol", npol)

  ; plasma profiles

  s = file_write(handle, "pressure", pressure)
  s = file_write(handle, "Jpar0", Jpar0)
  s = file_write(handle, "Ni0", Ni0)
  s = file_write(handle, "Te0", Te0)
  s = file_write(handle, "Ti0", Ti0)
  s = file_write(handle, "Ni_x", Ni)
  s = file_write(handle, "Te_x", Ti)
  s = file_write(handle, "Ti_x", Ti)
  s = file_write(handle, "bmag", Bt0)
  s = file_write(handle, "rmag", Rmaj)

  ; Curvature
  s = file_write(handle, "logB", logB)
;  s = file_write(handle, "bxcvx", bxcvx)
;  s = file_write(handle, "bxcvy", bxcvy)
;  s = file_write(handle, "bxcvz", bxcvz)

  ; Psi range
;  s = file_write(handle, "psi_axis", mesh.faxis)
;  s = file_write(handle, "psi_bndry", psi_bndry)

  file_close, handle
  PRINT, "DONE"
END

