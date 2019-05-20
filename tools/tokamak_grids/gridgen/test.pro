WINDOW, xsize=600, ysize=800

;g = read_neqdsk("efit/neqdsk")
g = read_neqdsk("efit/g014220.00200")

rz_grid = {nr:g.nx, nz:g.ny, $  ; Number of grid points
                   r:REFORM(g.r[*,0]), z:REFORM(g.z[0,*]), $  ; R and Z as 1D arrays
                   simagx:g.simagx, sibdry:g.sibdry, $ ; Range of psi
                   psi:g.psi, $  ; Poloidal flux in Weber/rad on grid points
                   npsigrid:(FINDGEN(N_ELEMENTS(g.pres))/(N_ELEMENTS(g.pres)-1)), $ ; Normalised psi grid for fpol, pres and qpsi
                   fpol:g.fpol, $ ; Poloidal current function on uniform flux grid
                   pres:g.pres, $ ; Plasma pressure in nt/m^2 on uniform flux grid
                   qpsi:g.qpsi, $ ; q values on uniform flux grid
                   nlim:g.nlim, rlim:g.xlim, zlim:g.ylim} ; Wall boundary

settings = {nrad:64, npol:64, psi_inner:0.8D, psi_outer:1.1D}
boundary = TRANSPOSE([[rz_grid.rlim], [rz_grid.zlim]])
mesh = create_grid(rz_grid.psi, rz_grid.r, rz_grid.z, settings, boundary=boundary, /strict, /simple)

;process_grid, rz_grid, mesh

