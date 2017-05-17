; SLAB GRID GENERATOR. EVERYTHING IN SI UNITS 
;
; NX - number of points in radial direction. 
;      For radial need nx = m * nproc + 4
;
; NY - Number of points along field-line
;
; NZ - Number of points in symmetry direction.
;      Set in options file, not grid
;
; LX - Height of box [m]
; LY - Length of box along field-line [m]
; LZ - Length of box in symmetry direction [m]
;
; p0 - Maximum pressure [Pa]
; pedge - Minimum pressure [Pa]
; pwidth - Width of pedestal [m]
;
; density - Constant density everywhere [m^-3]
; B0 - Constant B field (must be non-zero) [T]

PRO generate
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; SETTINGS

    nx = 32 ;64
    ny = 32 ;64
    nz = 33
    output = "slab_grid.nc"

    Lx = 0.6      ;0.6 is standard, have increased it to move edges of box away from interesting region
    Ly = 5.4
    Lz = 0.06
    
    xmin=0.3-Lx/2.0
    
    ; Grid spacing
    dx = FLTARR(nx, ny) + Lx / FLOAT(nx-1)
    dy = FLTARR(nx, ny) + Ly / FLOAT(ny-1)

  ;  pedge = p0 / 100.  ; Just need to prevent negatives    
  ;  pwidth = Lx / 10.
    
  ;  density = 1e20  ; NOTE: SHOULD BE SAME AS IN BOUT.inp FILE
  ;  B0 = 0.01

  ;  gravity_sin = 2e4;1e4  ; sin component of gravity [m/s^2]
   ; gravity_cos = 0.0;1e4  ; cos component of gravity [m/s^2]

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    MU0 = 1.0 ;4.e-7*!PI
    
   ; rho = density * 2.0*1.6726e-27 ; Deuterium mass density

    F2D = FLTARR(nx, ny)
    F3D = FLTARR(nx,ny,nz)


    	rho_a = 2.3
	l_rho0 = 2.5
	rho_b = 0.5
	x_rho0 = 0.6
	l_rho1 = 0.4
	p0 = 0.3
	p_b = 0.8
	x_p0 = 0.15
	l_p1 = 0.2
        gravity = 0.25 ;0.21

        s = 0.0	;shear
;        Ls = 1.0/s	;inverse shear
        xs0 = 0.32	;position of straight field line
    
    ; Set pressure profile. Function of x only
    x = Lx * FINDGEN(nx) / FLOAT(nx-1) + xmin ; Distance [m]
    y = Ly * FINDGEN(ny) / FLOAT(ny-1) -Ly/2; 
    z = Lz * FINDGEN(nz) / FLOAT(nz-1) -Lz/2; 
;    z = Lz * FINDGEN(nz) / FLOAT(nz-1);

   B=FLTARR(nx)
   Bz_prof=FLTARR(nx)
   By_prof=FLTARR(nx)
   Bx_prof=FLTARR(nx)

;--------initialise fields-------------------------------------------------------------------

    p = p0*1 + p0*p_b*exp(-(x-x_p0)*(x-x_p0)/(l_p1*l_p1))
    
    rho = rho_a*exp(-x/l_rho0) + rho_b*exp(-(x-x_rho0)*(x-x_rho0)/(l_rho1*l_rho1))


    ; Put into 2D function
    pressure = F2D
    FOR i=0, nx-1 DO pressure[i,*] = p[i]
    
    density = F2D
    for i=0, nx-1 DO density[i,*] = rho[i]


    ;magnetic field

    for i =0,nx-1 do B[i] =sqrt(mu0*(1+2*p[0] - 2*p[i] -2*gravity*(-rho_a*l_rho0*(exp(-x[i]/l_rho0) -1) + rho_b*l_rho1*(sqrt(!PI))*0.5*(erf((x[i]-x_rho0)/l_rho1) - erf(-x_rho0/l_rho1)))))
    

    modB=F2D
    for i = 0, nx-1 DO modB[i,*] = B[i]


    zmag = fltarr(nx)
    norm = fltarr(nx)

    FOR i =0, nx-1 DO BEGIN 

      zmag[i] = (x[i]-xs0)*s
      norm[i] = sqrt(1+zmag[i]*zmag[i])

      Bx_prof[i] = 0.0
      By_prof[i] = B[i]/norm[i]
      Bz_prof[i] = B[i]*zmag[i]/norm[i]
    ENDFOR


    Bx=F2D
    By=F2D
    Bz=F2D

    for i=0, nx-1 DO By[i,*] = By_prof[i]
    for i=0, nx-1 DO Bx[i,*] = Bx_prof[i]
    for i=0, nx-1 DO Bz[i,*] = Bz_prof[i]



   Jpar0=F2D
   jpar=fltarr(nx)

   for i=0, nx-1 DO jpar[i] = s*B[i]/(MU0*(1+s*s*x[i]*x[i]))

   for i=0, nx-1 DO Jpar0[i,*] = jpar[i]
;--------------------------------------------------------------------------------------------------



;------------initialise velocities------------------------------------------------------------------
   Vx0 = 1.0E-3
   Vy0 = 1.0E-3
   Vz0 = 1.0E-3
   Vpar0= 1.0E-3
   x0 = .3112
   ax1 = 138.66885708
   kz = 104.7197551196597
   ky = 0.7265193360796418
   beta = 0.59089827
   delta = 1.861762821403754



   Vx=F3D
   Vy=F3D
   Vz=F3D

   V0x=0
   V0y=0
   V0z=0


   for i=0,nx-1 do begin
      for j=0,ny-1 do begin
         for k=0,nz-1 do begin
            Vx(i,j,k) = Vx0 * exp(-ax1*(x[i]-x0)*(x[i]-x0)) * cos(kz*(z[k]-s*(x[i]-xs0)*y[j])) * (1-cos(ky*y[j])/cos(ky*ly/2))/(1-1/cos(ky*ly/2)) 
         endfor
      endfor
   endfor




   for i=0,nx-1 do begin
      for j=0,ny-1 do begin
         for k=0,nz-1 do begin
            Vz(i,j,k) = Vz0 * exp(-ax1*(x[i]-x0)*(x[i]-x0))*(  2.0* ax1*(x[i]-x0) * sin(kz*(z[k]-s*(x[i]-xs0)*y[j])) *(1/ky) +s*y[j]*cos(kz*(z[k]-s*(x[i]-xs0)*y[j])))* (1-cos(ky*y(j))/cos(ky*ly/2))/(1-1/cos(ky*ly/2)) 
         endfor
      endfor
   endfor

   for i=0,nx-1 do begin
      for j=0,ny-1 do begin
         for k=0,nz-1 do begin

             Vy(i,j,k) = Vy0*(-(s*(x[i]-xs0)/Vz0)*Vz(i,j,k) + sqrt(1+s*s*(x[i]-xs0)*(x[i]-xs0))* exp(- ax1 * (x[i] - x0) * (x[i] - x0)) * cos(kz*(z[k]-s*(x[i]-xs0)*y[j])) *2* delta / ly * ( tan(ky * (ly/2)) / (ky * (ly/2)) * (y[j] + (ly/2)) - (sin(ky * y[j]) + sin(ky * (ly/2))) / (ky * cos(ky * (ly/2)) )) )/ (1. - 1. / cos(ky * (ly/2)))
;            Vy(i,j,k) = Vy0 * exp(- ax1 * (x[i] - x0) * (x[i] - x0)) * cos(kz * z[k]) *2* delta / ly * ( tan(ky * (ly/2)) / (ky * (ly/2)) * (y[j] + (ly/2)) - (sin(ky * y[j]) + sin(ky * (ly/2))) / (ky * cos(ky * (ly/2)) )) / (1. - 1. / cos(ky * (ly/2)))
         endfor
      endfor
   endfor


  Vpar=F3D

    for i=0,nx-1 do begin
      for j=0,ny-1 do begin
         for k=0,nz-1 do begin

             Vpar(i,j,k) = Vpar0*(exp(- ax1 * (x[i] - x0) * (x[i] - x0)) * cos(kz*(z[k]-s*(x[i]-xs0)*y[j])) *2* delta / ly * ( tan(ky * (ly/2)) / (ky * (ly/2)) * (y[j] + (ly/2)) - (sin(ky * y[j]) + sin(ky * (ly/2))) / (ky * cos(ky * (ly/2)) )) )/ (1. - 1. / cos(ky * (ly/2)))

         endfor
      endfor
   endfor


  phi=F3D

   for i=0,nx-1 do begin
      for j=0,ny-1 do begin
         for k=0,nz-1 do begin
            phi(i,j,k) = -Vx0 * exp(-ax1*(x[i]-x0)*(x[i]-x0)) *(1.0/kz) *sin(kz*(z[k]-s*(x[i]-xs0)*y[j])) * (1-cos(ky*y[j])/cos(ky*ly/2))/(1-1/cos(ky*ly/2)) 
         endfor
      endfor
   endfor

;-----------------------------------------------------------------------------------------------------------------------

;------------Funk around with FFTs so that bout knows z variations properly---------------------------------------------

  nf=fix(nz/2)

  IF nf mod 2 NE 1 then nf=nf-1

  vfx=fltarr(nx,ny,nf)

  for i=0, nx-1 do begin
     for j=0,ny-1 do begin
        fx = FFT(reform(Vx[i,j,*]))
      
        Vfx[i,j,0] = real_part(fx[0])
        for k=0,(nf-1)/2 - 1 do begin
           Vfx[i,j,2*k+1] = real_part(fx[k+1])
           Vfx[i,j,2*k+2] = imaginary(fx[k+1])
        endfor
     endfor
  endfor


  vfy=fltarr(nx,ny,nf)

  for i=0, nx-1 do begin
     for j=0,ny-1 do begin
        fy = FFT(reform(Vy[i,j,*]))
      
        Vfy[i,j,0] = real_part(fy[0])
        for k=0,(nf-1)/2 - 1 do begin
           Vfy[i,j,2*k+1] = real_part(fy[k+1])
           Vfy[i,j,2*k+2] = imaginary(fy[k+1])
        endfor
     endfor
  endfor


  vfz=fltarr(nx,ny,nf)

  for i=0, nx-1 do begin
     for j=0,ny-1 do begin
        fz = FFT(reform(Vz[i,j,*]))
      
        Vfz[i,j,0] = real_part(fz[0])
        for k=0,(nf-1)/2 - 1 do begin
           Vfz[i,j,2*k+1] = real_part(fz[k+1])
           Vfz[i,j,2*k+2] = imaginary(fz[k+1])
        endfor
     endfor
  endfor



  vfpar=fltarr(nx,ny,nf)

  for i=0, nx-1 do begin
     for j=0,ny-1 do begin
        for k=0,nf-1 do begin
        vpar[i,j,k]= 0

        endfor
     endfor
  endfor

  for i=0, nx-1 do begin
     for j=0,ny-1 do begin
        fpar = FFT(reform(Vpar[i,j,*]))
      
        Vfpar[i,j,0] = real_part(fpar[0])
           k=0
           ;for k=0,(nf-1)/2 - 1 do begin
           Vfpar[i,j,2*k+1] = real_part(fpar[k+1])
           Vfpar[i,j,2*k+2] = imaginary(fpar[k+1])
        ;endfor
     endfor
  endfor

  phif=fltarr(nx,ny,nf)

  for i=0, nx-1 do begin
     for j=0,ny-1 do begin
        for k=0,nf-1 do begin
        phif[i,j,k]= 0

        endfor
     endfor
  endfor

  for i=0, nx-1 do begin
     for j=0,ny-1 do begin
        fphi = FFT(reform(phi[i,j,*]))
      
        phif[i,j,0] = real_part(fphi[0])
         k=0
        ;for k=0,(nf-1)/2 - 1 do begin
           phif[i,j,2*k+1] = real_part(fphi[k+1])
           phif[i,j,2*k+2] = imaginary(fphi[k+1])
        ;endfor
     endfor
  endfor

;------------------------------------------------------------------------------------------------------------------------
;-----------gravity----------------------------------------------------------------------------------------------------
    G = fltarr(nx)
    G = -gravity*x

    Garr=F2D
    FOR i=0, nx-1 DO Garr[i,*] = G[i]


;--------------------------initialise vorticity--------------------------------------------------------

  U0=F3D
  fst=F3D
  sec=F3D
  thd=F3D
  fth=F3D
  fvt=F3D


   for i=0,nx-1 do begin
      for j=0,ny-1 do begin
         for k=0,nz-1 do begin
            fst[i,j,k] = 2*ax1*x[i]*x[i]*s*Vy[i,j,k] - 2*s*ax1*x[i]*Vz[i,j,k] + kz*Vx[i,j,k]*tan(kz*(z[k]-s*(x[i]-xs0)*y[j]))
            sec[i,j,k] = (-(1/kz)*4*ax1*x[i]*x[i]*s*s + kz*s*s*s*s*x[i]*x[i]*y[j]*y[j] + s*s*kz*(x[i]*x[i] + y[j]*y[j]) * 2*ax1/kz)*sin(kz*(z[k]-s*(x[i]-xs0)*y[j]))
            thd[i,j,k] = (-2*ax1*s*x[i]*y[j]*(1+s*s*x[i]*x[i]) +s*s*s*x[i]*y[j] ) * cos(kz*(z[k]-s*(x[i]-xs0)*y[j]))
            fth[i,j,k] = sqrt(1+s*s*x[i]*x[i])*s*s*x[i]*y[j]*delta*(2*kz/(ky*ly))*sin(kz*(z[k]-s*(x[i]-xs0)*y[j])) + s*s*s*x[i]*x[i]*(2*delta/(ky*ly*sqrt(1+s*s*x[i]*x[i])))*cos(kz*(z[k]-s*(x[i]-xs0)*y[j]))
            fvt[i,j,k] = s*x[i]*cos(kz*(z[k]-s*(x[i]-xs0)*y[j]))
         endfor
      endfor
   endfor



   for i=0, nx-1 do begin
     for j=0, ny-1 do begin
       for k=0, nz-1 do begin

          U0[i,j,k]= (  fst[i,j,k] + exp(-ax1*(x[i]-x0)*(x[i]-x0))*((sec[i,j,k]+thd[i,j,k])*(1 - cos(ky*y[j])/cos(ky*ly/2)) + fth[i,j,k]*((2*y[j]/ly)*tan(ky*ly/2) - sin(ky*y[j])/cos(ky*ly/2)) + fvt[i,j,k]*(ky*sin(ky*y[j])/cos(ky*ly/2)))*(1/(1-1/(cos(ky*ly/2))))   ) / (sqrt(1+s*s*x[i]*x[i]))

       endfor
     endfor
   endfor


  Uf=fltarr(nx,ny,nf)

  for i=0, nx-1 do begin
     for j=0,ny-1 do begin
        fu = FFT(reform(U0[i,j,*]))
      
        Uf[i,j,0] = real_part(fu[0])
        for k=0,(nf-1)/2 - 1 do begin
           Uf[i,j,2*k+1] = real_part(fu[k+1])
           Uf[i,j,2*k+2] = imaginary(fu[k+1])
        endfor
     endfor
  endfor


;----------------------------------------------------------------------------------------------------------
;----------------------------------------------------------------------------------------------------------
;------------------------------------ Analytical initial perturbations for density, pressure and Psi------


;-----Pressure

  P1 = F3D

  for i =0, nx-1 do begin
    for j = 0, ny-1 do begin
      for k = 0, nz-1 do begin
        P1[i,j,k] = ((1./kz)*(-p0)*p_b*2.0*(1/(l_p1*l_p1))*(x[i]- x_p0)*exp(-(x[i]-x_p0)*(x[i]-x_p0)/(l_p1*l_p1)) * exp(-ax1*x[i]*x[i])*kz*sin(kz*z[k])*(1-cos(ky*y[j])/cos(ky*ly/2)) - (5./3.)*(p[i]/(1+(5./3.)*p[i]*mu0/(B[i]*B[i])))*( (-(1./kz)*rho[i]*mu0*gravity*(1/(B[i]*B[i]))*exp(-ax1*x[i]*x[i])*sin(kz*z[k])*(1-cos(ky*y[j])/cos(ky*ly/2))) +exp(-ax1*x[i]*x[i])*cos(kz*z[k])*delta*(2./(ky*ly))*(tan(ky*ly/2.)-ky*cos(ky*y[j])/(cos(ky*ly/2.0)))) )/ (1. - 1. / cos(ky * (ly/2)))
      endfor
    endfor
  endfor


  Pf=fltarr(nx,ny,nf)

  for i=0, nx-1 do begin
     for j=0,ny-1 do begin
        fp = FFT(reform(P1[i,j,*]))
      
        Pf[i,j,0] = real_part(fp[0])
        for k=0,(nf-1)/2 - 1 do begin
           Pf[i,j,2*k+1] = real_part(fp[k+1])
           Pf[i,j,2*k+2] = imaginary(fp[k+1])
        endfor
     endfor
  endfor



;-----Density

  D1 = F3D

  for i =0, nx-1 do begin
    for j = 0, ny-1 do begin
      for k = 0, nz-1 do begin
        D1[i,j,k] = ((1./kz)*((-rho_a/l_rho0)*exp(-x[i]/l_rho0) - 2.0*rho_b*(1/(l_rho1*l_rho1)))*(x[i]-x_rho0)*exp(-(x[i]-x_rho0)*(x[i]-x_rho0)/(l_rho1*l_rho1)) * exp(-ax1*x[i]*x[i])* kz*sin(kz*z[k])*(1-cos(ky*y[j])/cos(ky*ly/2))   - (rho[i]/(1+(5./3.)*p[i]*mu0/(B[i]*B[i])))*( (-(1./kz)*rho[i]*mu0*gravity*(1/(B[i]*B[i]))*exp(-ax1*x[i]*x[i])*sin(kz*z[k])*(1-cos(ky*y[j])/cos(ky*ly/2))) +exp(-ax1*x[i]*x[i])*cos(kz*z[k])*delta*(2./(ky*ly))*(tan(ky*ly/2.)-ky*cos(ky*y[j])/(cos(ky*ly/2.0)))) )/ (1. - 1. / cos(ky * (ly/2)))
      endfor
    endfor
  endfor


  Df=fltarr(nx,ny,nf)

  for i=0, nx-1 do begin
     for j=0,ny-1 do begin
        fd = FFT(reform(D1[i,j,*]))
      
        Df[i,j,0] = real_part(fd[0])
        for k=0,(nf-1)/2 - 1 do begin
           Df[i,j,2*k+1] = real_part(fd[k+1])
           Df[i,j,2*k+2] = imaginary(fd[k+1])
        endfor
     endfor
  endfor




;-----Psi

  psi1 = F3D

  for i =0, nx-1 do begin
    for j = 0, ny-1 do begin
      for k = 0, nz-1 do begin
        psi1[i,j,k] = -(1./kz)* exp(-ax1*x[i]*x[i])* sin(kz*z[k])*(sin(ky*y[j])/cos(ky*ly/2))/ (1. - 1. / cos(ky * (ly/2)))
      endfor
    endfor
  endfor


  Psif=fltarr(nx,ny,nf)

  for i=0, nx-1 do begin
     for j=0,ny-1 do begin
        fpsi = FFT(reform(psi1[i,j,*]))
      
        psif[i,j,0] = real_part(fpsi[0])
        for k=0,(nf-1)/2 - 1 do begin
           psif[i,j,2*k+1] = real_part(fpsi[k+1])
           psif[i,j,2*k+2] = imaginary(fpsi[k+1])
        endfor
     endfor
  endfor
;-------------------------------------------------------------------------------------------------------------------------

;---- TOPOLOGY

;   ixseps= -1  ; Domain outside core -> NOT periodic

    ;Ben says define xmesh such that B = grad_z cross grad_x - this basically involves multiplying our xmesh by B.


    xmeshB = fltarr(nx);

    for i = 0,nx-1 do begin
      xmeshB[i] = x[i]*B[i];
    endfor

;--------------------------- Write to file---------------------------------------------------------------------------
    fp = file_open(output, /create)
    
    status = file_write(fp, "nx", nx)
    status = file_write(fp, "ny", ny)
;    status = file_write(fp, "nz", nz) ; Don't write nz if using FFT in Z
    
    status = file_write(fp, "dx", dx)
    status = file_write(fp, "dy", dy)
;     status = file_write(fp, "ixseps1", ixseps)
;     status = file_write(fp, "ixseps2", ixseps)
    
    status = file_write(fp, "p0", pressure)
    status = file_write(fp, "rho0", density)
    status = file_write(fp, "x", xmeshB)
    status = file_write(fp, "y", y)
    status = file_write(fp, "z", z)
    status = file_write(fp, "B0_vec_x", Bx)
    status = file_write(fp, "B0_vec_y", By)
    status = file_write(fp, "B0_vec_z", Bz)
    status = file_write(fp, "v0xtest", Vx)
    status = file_write(fp, "v0ytest", Vy)
    status = file_write(fp, "v0ztest", Vz)
;     status = file_write(fp, "gravity_x", gravityx)
    status = file_write(fp, "v0_x", Vfx)
    status = file_write(fp, "v0_y", Vfy)               ;if want covariant
    status = file_write(fp, "v0_z", Vfz)
    status = file_write(fp, "G", Garr)
    status = file_write(fp, "Jpar0", Jpar0)
    status = file_write(fp, "U0_test", U0)
    status = file_write(fp, "U0", Uf)
    status = file_write(fp, "B0", modB)
    status = file_write(fp, "Vpar0", Vfpar)
    status = file_write(fp, "psi1", psif)
    status = file_write(fp, "rho1", Df)               ;if want covariant
    status = file_write(fp, "p1", Pf)
    status = file_write(fp, "p1test", P1)
    status = file_write(fp, "Vpar0test", Vpar)
    status = file_write(fp, "psi1test", psi1)
    status = file_write(fp, "rho1test", D1)               ;if want covariant
    status = file_write(fp, "U0test", U0)
    status = file_write(fp, "phi0", phif)
;     status = file_write(fp, "v0x", Vfx)
;     status = file_write(fp, "v0y", Vfy)                ;if want contravariant
;     status = file_write(fp, "v0z", Vfz)
;    status = file_write(fp, "gravity_x", gfx)
;     status = file_write(fp, "gravity_x", gravity) 

    ; I is integrated shear, and zero here
    
    file_close, fp
;-------------------------------------------------------------------------------------------------------------------------
END
