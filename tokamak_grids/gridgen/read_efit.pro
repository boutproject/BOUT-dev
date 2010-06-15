
FUNCTION read_efit, dir=dir, afile=afile, gfile=gfile, $
                    adata=adata, gdata=gdata
  IF NOT KEYWORD_SET(dir) THEN dir="."
  IF NOT KEYWORD_SET(afile) THEN afile=dir+"/aeqdsk"
  IF NOT KEYWORD_SET(gfile) THEN gfile=dir+"/neqdsk"
  
  gdata = read_neqdsk(gfile)
  adata = read_aeqdsk(afile)
  
  ; Put main results into structure
  
  efd = {nx:gdata.nx, ny:gdata.ny, $ ; Number of grid points
         r:r, z:z, $ ; R, Z locations of grid points
         psi:gdata.psi, $ 
         xdim:gdata.xdim, zdim:gdata.zdim, $ ; Size of the domain in meters
         rcentr:gdata.rcentr, bcentr:gdata.bcentr, $ ; R, Bt reference vacuum field
         rgrid1:gdata.rgrid1, $
         zmid:gdata.zmid, $
         rmagx:gdata.rmagx, zmagx:gdata.zmagx, $
         simagx:gdata.simagx, sibdry:gdata.sibdry, $
         fpol:gdata.fpol, $  ; Poloidal current function on uniform flux grid
         pres:gdata.pres, $  ; Plasma pressure in nt/m^2 on uniform flux grid
         qpsi:gdata.qpsi, $  ; q values on uniform flux grid
        }
         
  RETURN, efd
END
