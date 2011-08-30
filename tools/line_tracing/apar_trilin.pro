; Trilinear interpolation for perturbed parallel vector potential
; Input is x,y,z coordinates
; Output is structure ={x:dA/dx,y:dA/dy,z:dA,dz}
FUNCTION apar_trilin, x, y, z, debug=debug, shift=shift

  COMMON BDATA, bd ;-bd={apar:apar, nx:nx, ny:ny, nz:nz}
  COMMON griddata, g, deltaZtor, Ntor
  
  yMin=0.0
  yMax=2*!PI
  xMin=MIN(g.psixy)
  xMax=MAX(g.psixy)
  zmin=0.0
  zmax=deltaZtor                ;;2*!PI/period
  nz=bd.nz

;Find the lower left corner of the cell
  xind = INTERPOL(FINDGEN(g.nx),g.psixy[*,0],x)
  ; must be 0,...,g.nx-2; if ix_int>=g.nx-1 then point is outside domain
  ; should return with an error message if x is larger or smaller than that
  IF ((xind<0) OR (xind GE g.nx-1)) THEN STOP
  ix_int = FIX(xind)

  ; assumes that y is between 0 and 2*!PI
  yind = DOUBLE(g.ny)*(y-yMin)/(yMax-yMin) ;gives floats between 0 and g.ny
  iy_int = FIX(yind)                       ;NOTE that if y=2*!PI exactly, then we 
  IF iy_int EQ g.ny THEN STOP         ;need to take twist-shift into account!

  ; assumes that z is between 0 and deltaZtor
  zind = DOUBLE(nz)*(z-zMin)/(zMax-zMin)   ;gives floats between 0 and nz
  iz_int = FIX(zind) MOD nz                ;if z>=deltaZtor exactly, use 0->delta_z

;Calculate grid spacing
  dx=g.psixy[ix_int+1,0]-g.psixy[ix_int,0]
  dy=(yMax-yMin)/DOUBLE(g.ny)
  dz=(zMax-zMIN)/DOUBLE(nz) 
    
  ;-calculate partial derivatives within cubic cell
  cellVal=bd.Apar[ix_int:(ix_int+1),iy_int:(iy_int+1),iz_int:(iz_int+1)]

;;;;Using Trilinear interpolation
  x_0=g.psixy[ix_int,0]
  xd=(x-x_0)/dx    ;normalized x 
  y_0=iy_int*dy
  yd=(y-y_0)/dy    ;normalized y
  z_0=iz_int*dz
  zd=(z-z_0)/dz    ;normalized z

  F_xyz=           ( cellval[0,0,0]*(1.0-xd)*(1.0-yd)*(1.0-zd)                 $
                    +cellval[1,0,0]*     xd *(1.0-yd)*(1.0-zd)                 $
                    +cellval[0,1,0]*(1.0-xd)*     yd *(1.0-zd)                 $
                    +cellval[1,1,0]*     xd *     yd *(1.0-zd)                 $
                    +cellval[0,0,1]*(1.0-xd)*(1.0-yd)*     zd                  $
                    +cellval[1,0,1]*     xd *(1.0-yd)*     zd                  $
                    +cellval[0,1,1]*(1.0-xd)*     yd *     zd                  $
                    +cellval[1,1,1]*     xd *     yd *     zd  )
  dAdx_tri=1.0d/dx*(-cellval[0,0,0]*(1.0-yd)*(1.0-zd)                          $
                    +cellval[1,0,0]*(1.0-yd)*(1.0-zd)                          $
                    -cellval[0,1,0]*     yd *(1.0-zd)                          $
                    +cellval[1,1,0]*     yd *(1.0-zd)                          $
                    -cellval[0,0,1]*(1.0-yd)*     zd                           $
                    +cellval[1,0,1]*(1.0-yd)*     zd                           $
                    -cellval[0,1,1]*     yd *     zd                           $
                    +cellval[1,1,1]*     yd *     zd  )
  dAdy_tri=1.0d/dy*(-cellval[0,0,0]*(1.0-xd)*(1.0-zd)                          $
                    -cellval[1,0,0]*     xd *(1.0-zd)                          $
                    +cellval[0,1,0]*(1.0-xd)*(1.0-zd)                          $
                    +cellval[1,1,0]*     xd *(1.0-zd)                          $
                    -cellval[0,0,1]*(1.0-xd)*     zd                           $
                    -cellval[1,0,1]*     xd *     zd                           $
                    +cellval[0,1,1]*(1.0-xd)*     zd                           $
                    +cellval[1,1,1]*     xd *     zd  )
  dAdz_tri=1.0d/dz*(-cellval[0,0,0]*(1.0-xd)*(1.0-yd)                          $
                    -cellval[1,0,0]*     xd *(1.0-yd)                          $
                    -cellval[0,1,0]*(1.0-xd)*     yd                           $
                    -cellval[1,1,0]*     xd *     yd                           $
                    +cellval[0,0,1]*(1.0-xd)*(1.0-yd)                          $
                    +cellval[1,0,1]*     xd *(1.0-yd)                          $
                    +cellval[0,1,1]*(1.0-xd)*     yd                           $
                    +cellval[1,1,1]*     xd *     yd  )

  res={f:F_xyz, x:dAdx_tri, y:dAdy_tri, z:dAdz_tri}

  IF KEYWORD_SET(debug) THEN STOP
  RETURN, res
END
