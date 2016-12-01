; Generate cyclone test cases
;
; A.M.Dimits et. al. Phys. Plasmas 7 (2000) p969
;
; Concentric circle "equilibrium" for local flux-surface calculations
; 

; Integrate 1 / (1 - eps*cos(theta))^2  where eps = r / R
FUNCTION eps_func, theta, y
  COMMON eps_com, eps
  RETURN, 1. / ((1. - eps*COS(theta))^2)
END

FUNCTION eps_integral, eps, theta=theta
  COMMON eps_com, ep
  ep = eps
  IF NOT KEYWORD_SET(theta) THEN theta = 2.*!PI
  result = [0.]
  return, (LSODE(result, 0.0, theta, 'eps_func'))[0]
END

PRO cyclone, output=output, varyBp=varyBp

  nx = 68 ; Radial grid points
  ny = 32 ; Poloidal (parallel) grid points
  
  IF NOT KEYWORD_SET(output) THEN output = "cyclone_"+STR(nx)+"x"+STR(ny)+".nc"

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ; Ion density in 10^20 m^-3 
  ni = 1.
  
  ; Temperature in eV (Te = Ti)
  Ti = 1000
  
  ; Major radius [meters]
  Rmaj = 4

  ; Safety factor q = r*Bt/(R*Bp)
  q = 1.4
  
  ; Magnetic shear s = (r/q) dq/dr
  s = 0.776
  
  ; Ratio of density to temp. length scales eta = L_n / L_T
  eta_i = 3.114
  
  ; Inverse aspect ratio epsilon = r / R
  epsilon = 0.18
  
  ; Ratio of major radius to L_T  Rnorm  = R / L_T
  Rnorm = 6.92
  
  ; Normalised ion gyro-radius rho_norm = rho_i / L_T
  rho_norm = 0.01
  
  ; Radial extent, normalised to gyro-radius r_wid = dr / rho_i
  r_wid = 100
  
  Mi = 2.*1.67262158e-27   ; Ion mass [kg]. Deuterium
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  rminor = Rmaj * epsilon  ; Minor radius [m]
  L_T = Rmaj / Rnorm       ; Temp. length scale [m]
  L_n = eta_i * L_T        ; Density length scale [m]
  rho_i = rho_norm * L_T   ; Ion Larmor radius [m]
  Bt0 = SQRT(2.*Ti*Mi / 1.602e-19) / rho_i ; Toroidal field from rho_i [T]
  Bp = rminor * Bt0 * eps_integral(epsilon)/ (2*!PI * q * Rmaj) ; Poloidal field [T]
  
  dr = r_wid * rho_i       ; Width of domain [m]
  
  theta = 2.*!PI * FINDGEN(ny) / FLOAT(ny)
  
  Rxy = FLTARR(nx, ny)
  FOR i=0, ny-1 DO Rxy[*,i] = Rmaj - rminor*COS(theta[i])
  Zxy = FLTARR(nx, ny)
  FOR i=0, ny-1 DO Zxy[*,i] = rminor*SIN(theta[i])
  
  dy = FLTARR(nx, ny) +  2.*!PI / FLOAT(ny)
  hthe = FLTARR(nx, ny) + rminor

  Btxy = Bt0 * Rmaj / Rxy
  PRINT, "Toroidal field varies from "+STR(Bt0*Rmaj/(Rmaj + rminor)) + $
         " to "+STR(Bt0*Rmaj/(Rmaj - rminor))
  ; Minor radius offset
  drprof = dr*((FINDGEN(nx) / FLOAT(nx-1)) - 0.5)
  ; q profile
  qprof = q + (s*q/rminor) * drprof
  PRINT, "q varies from "+STR(MIN(qprof))+" to "+STR(MAX(qprof))
  ShiftAngle = qprof * 2.*!PI
  
  IF KEYWORD_SET(varyBp) THEN BEGIN
     ; Vary Bpxy to get shear
     Bpxy = FLTARR(nx, ny)
     FOR i=0, ny-1 DO Bpxy[*,i] = Bp * q / qprof
     PRINT, "Poloidal field varies from "+STR(MIN(Bpxy))+" to "+STR(MAX(Bpxy))
  ENDIF ELSE BEGIN
     ; Constant Bp, but shift angle varies
     Bpxy = FLTARR(nx, ny) + Bp
  ENDELSE
  
  dx = Bp * (dr / FLOAT(nx-1)) * Rxy

  Bxy = SQRT(Btxy^2 + Bpxy^2)
  
  zShift = FLTARR(nx, ny)
  qint = eps_integral(epsilon)
  FOR i=1, ny-1 DO zShift[*,i] = ShiftAngle * eps_integral(epsilon, theta=theta[i]) / qint
  
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
  s = file_write(handle, "dy", dy)
  
  s = file_write(handle, "ShiftAngle", ShiftAngle)
  s = file_write(handle, "zShift", zShift)
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

