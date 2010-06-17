; generates an input file for Orszag-Tang vortex problem

PRO generate, ny=ny, mxg=mxg, file=file
  IF NOT KEYWORD_SET(ny) THEN ny = 128
  IF NOT KEYWORD_SET(mxg) THEN mxg = 2
  IF NOT KEYWORD_SET(file) THEN file="data/otv.grd.pdb"

  nx = 1 + 2*FIX(mxg)

  dy = 1.0 / FLOAT(ny)

  dxarr = FLTARR(nx, ny) + 1.0;
  dyarr = FLTARR(nx, ny) + dy;
  
  yarr = FINDGEN(ny) * dy

  rho = FLTARR(nx, ny) + 25./(36.*!PI)
  P = FLTARR(nx, ny) + 5./(12.*!PI)
  
  v_z = FLTARR(nx, ny)
  Bz = FLTARR(nx, ny)

  FOR y=0, ny-1 DO BEGIN
      v_z[*,y] = SIN(2.*!PI*yarr[y])
      Bz[*,y] = SIN(4.*!PI*yarr[y])
  ENDFOR
  Bz = Bz / SQRT(4.*!PI)

  ; domain inside core (periodic)

  ixseps1 = nx;
  ixseps2 = nx;

  PD_write, file, "nx", nx
  PD_write, file, "ny", ny
  PD_write, file, "dx", dxarr
  PD_write, file, "dy", dyarr

  PD_write, file, "ixseps1", ixseps1
  PD_write, file, "ixseps2", ixseps2

  PD_write, file, "rho0", rho
  PD_write, file, "p0", P
  PD_write, file, "v0_z", v_z
  PD_write, file, "B0z", Bz
END
