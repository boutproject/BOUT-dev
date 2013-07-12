; Reads an R-Z EFIT equilibrium from IDAM
; Traces are:
;
; EFM_PSI(R,Z)  - psi as function of R and Z
; EFM_P(PSI)_(C)  - Pressure as a function of psi
; EFM_F(PSI)_(C)  - Current function f as a function of psi
; EFM_Q(PSI)_(C)  - Safety factor as function of psi
; EFM_PSI_AXIS    - psi on the magnetic axis
; EFM_PSI_BOUNDARY  - psi on the plasma boundary (separatrix)
; EFM_LIMITER(N)  - Number of points in the limiter boundary
; EFM_LIMITER(R)  - Major radius of limiter
; EFM_LIMITER(Z)  - Height of limiter


FUNCTION read_idam, shot, time
  d = getdata("EFM_PSI(R,Z)", shot)
  
  ; Find the nearest time point
  m = MIN(ABS(d.time - time), tind)
  PRINT, "Time index = "+STR(tind)+" t = " + STR(d.time[tind])
  
  psi = REFORM(d.data[*,*,tind])
  r = d.x
  z = d.y
  nr = N_ELEMENTS(r)
  nz = N_ELEMENTS(z)

  d = getdata("EFM_PSI_AXIS", shot)
  simagx = d.data[tind]
  
  d = getdata("EFM_PSI_BOUNDARY", shot)
  sibdry = d.data[tind]
  
  ; It seems that efm_p(psi) etc. can have different number of
  ; time-points to efm_psi (!)

  d = getdata("EFM_P(PSI)_(C)", shot)
  m = MIN(ABS(d.time - time), tind)
  PRINT, "Time index = "+STR(tind)+" t = " + STR(d.time[tind])
  
  pres = d.data[*,tind]
  npsigrid = d.x
  
  d = getdata("EFM_F(PSI)_(C)", shot)
  fpol = d.data[*,tind]
  
  d = getdata("EFM_Q(PSI)_(C)", shot)
  qpsi = d.data[*,tind]
  
  
  d = getdata("EFM_LIMITER(R)", shot)
  rlim = reform(d.data)
  nlim = N_ELEMENTS(rlim)
  
  d = getdata("EFM_LIMITER(Z)", shot)
  zlim = reform(d.data)
  
  ; Analyse critical points

  critical = analyse_equil(psi, r, z)

  ; Put into structure
  
  rz_grid = {nr:nr, nz:nz, $
             r:r, z:z, $
             simagx:simagx, sibdry:sibdry, $
             psi:psi, $
             npsigrid:npsigrid, $
             fpol:fpol, pres:pres, qpsi:qpsi, $
             nlim:nlim, rlim:rlim, zlim:zlim, $
             critical:critical}

  RETURN, rz_grid
END
