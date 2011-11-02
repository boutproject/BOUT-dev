; Create an input file for 2D advection-diffusion equation example
;=================================================================;

@pdb2idl

PRO generate_grid, $
nx=nx, ny=ny, rxy=rxy, zxy=zxy, metric=metric, core=core,$
file=file, save=save, plot=plot, debug=debug

;
;
;

  IF NOT KEYWORD_SET(nx) THEN nx = 30
  IF NOT KEYWORD_SET(ny) THEN ny = 32
  IF NOT KEYWORD_SET(file) THEN file="./slab.grd.nc"




  ;;-domain boundaries 
  Rt0=0.5     ;;-reference location [m]
  Rmin=Rt0-0.5 ;;-left boundary  [m]
  Rmax=Rt0+0.5 ;;-right boundary [m]
  Zmin=0.0     ;;-bottom boundary [m]
  Zmax=1.0     ;;-top boundary [m]

  dR=(Rmax-Rmin)/(nx-1)
  dZ=(Zmax-Zmin)/(ny-1)

  hthe0=(Zmax-Zmin)/(2*!PI)
  hthe = FLTARR(nx, ny) + hthe0



  ;;-actual coordinates are needed only for plotting
  Rxy=fltarr(nx,ny)
  Zxy=fltarr(nx,ny)

  dx=fltarr(nx,ny)
  dy=fltarr(nx,ny)





  ;;-set geometry and magnetic field on the grid
  for ix=0,nx-1 do begin
      for jy=0,ny-1 do begin

          Rxy[ix,jy]=Rmin+dR*ix
          Zxy[ix,jy]=Zmin+dZ*jy

          dx[ix,jy] = dR
          dy[ix,jy] = dZ
          ;;dy[ix,jy] = 2.0*!PI/ny

      endfor
  endfor



  ;;-set background 2D profiles
  V0 = FLTARR(nx, ny)



  Ln=1e10 ;;-density decay length [m]

  for ix=0,nx-1 do begin
      for jy=0,ny-1 do begin

          V0[ix,jy] = 1e0

      endfor
  endfor


  ;;-set elements of metric tensor
  g11=fltarr(nx,ny)
  g22=fltarr(nx,ny)
  g33=fltarr(nx,ny)
  g12=fltarr(nx,ny)
  g13=fltarr(nx,ny)
  g23=fltarr(nx,ny)

  g11[*,*] = dR^2
  g22[*,*] = dZ^2
  g33[*,*] = 1.0
  g12[*,*] = 0.0
  g13[*,*] = 0.0
  g23[*,*] = 0.0
  

  if keyword_set(CORE) then begin
      ;; entire domain inside 'core' - periodic
      ixseps1 = nx
      ixseps2 = nx
      jyseps1_1 = -1
      jyseps2_2 = ny-1
  endif else begin
      ;; Topology: Set all points outside separatrix
      ;; so NOT periodic in Y
      ixseps1 = 0
      ixseps2 = 0
      jyseps1_1 = -1
      jyseps2_2 = ny-1
  endelse
  jyseps1_2=ny/2
  jyseps2_1=ny/2


  if keyword_set(SAVE) then begin

      print, "Writing data to ", file

      f = file_open(file, /create)

      status = file_write(f, "nx", nx)
      status = file_write(f, "ny", ny)


      status = file_write(f, "Rxy", Rxy)
      status = file_write(f, "Zxy", Zxy)


      ;;-supply metric tensor if not calculating internally
      status = file_write(f, "g11", g11)
      status = file_write(f, "g22", g22)
      status = file_write(f, "g33", g33)
      status = file_write(f, "g12", g12)
      status = file_write(f, "g13", g13)
      status = file_write(f, "g23", g23)


      ;;-flags for setting topology (may be needed)
      status = file_write(f, "ixseps1", ixseps1)
      status = file_write(f, "ixseps2", ixseps2)
      status = file_write(f, "jyseps1_1", jyseps1_1)
      status = file_write(f, "jyseps2_2", jyseps2_2)
      status = file_write(f, "jyseps2_1", jyseps2_1)
      status = file_write(f, "jyseps1_2", jyseps1_2)


      ;;-initial profile
      status = file_write(f, "V0", V0)



      ;;-some auxiliary quantities (may be needed)
      status = file_write(f, "hthe", hthe)
      status = file_write(f, "hthe0", hthe0)

      status = file_write(f, "dx", dx)
      status = file_write(f, "dy", dy)


      file_close, f

  endif


;
;
;
if keyword_set(DEBUG) then STOP
END
