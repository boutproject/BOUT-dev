; Bilinear interpolation in x,y, cubic in z for Apar 
; Input is x,y,z coordinates
; Output is structure ={x:dA/dx,y:dA/dy,z:dA,dz}
FUNCTION apar_bilin_cube, x, y, z, debug=debug, shift=shift

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
  IF iy_int EQ g.ny THEN STOP              ;need to take twist-shift into account!

  ; assumes that z is between 0 and deltaZtor
  zind = DOUBLE(nz)*(z-zMin)/(zMax-zMin)+1   ;gives floats between 1 and nz+1
  iz_int = FIX(zind) MOD (nz+1)              ;if z>=deltaZtor exactly, use 0->delta_z

;Calculate grid spacing
  dx=g.psixy[ix_int+1,0]-g.psixy[ix_int,0]
  dy=(yMax-yMin)/DOUBLE(g.ny)
  dz=(zMax-zMIN)/DOUBLE(nz) 
    
;normalized cell coordinates
  x_0=g.psixy[ix_int,0]
  xd=(x-x_0)/dx    ;normalized x 
  y_0=iy_int*dy
  yd=(y-y_0)/dy    ;normalized y
  z_0=(iz_int-1)*dz
  zd=(z-z_0)/dz    ;normalized z

  ;calculate coefficients of each corner's cubic polynomial in z
  coeffs_00=interp_cube(bd.apar[ix_int  ,iy_int  ,(iz_int-1):(iz_int+2)])
  coeffs_10=interp_cube(bd.apar[ix_int+1,iy_int  ,(iz_int-1):(iz_int+2)])
  coeffs_01=interp_cube(bd.apar[ix_int  ,iy_int+1,(iz_int-1):(iz_int+2)])
  coeffs_11=interp_cube(bd.apar[ix_int+1,iy_int+1,(iz_int-1):(iz_int+2)])

  z_vec=[[zd^0],[zd^1],[zd^2],[zd^3]]
  zp_vec=[[0],[1],[2.*zd],[3.*zd^2]]

  Aval=  (1.0-xd)*(1.0-yd)*z_vec#coeffs_00                                     $
        +(    xd)*(1.0-yd)*z_vec#coeffs_10                                     $
        +(1.0-xd)*(    yd)*z_vec#coeffs_01                                     $
        +(    xd)*(    yd)*z_vec#coeffs_11

  dAdx=1.0d/dx*( -(1.0-yd)*z_vec#coeffs_00                                     $
                 +(1.0-yd)*z_vec#coeffs_10                                     $
                 -(    yd)*z_vec#coeffs_01                                     $
                 +(    yd)*z_vec#coeffs_11 ) 

  dAdy=1.0d/dy*( -(1.0-xd)*z_vec#coeffs_00                                     $
                 -(    xd)*z_vec#coeffs_10                                     $
                 +(1.0-xd)*z_vec#coeffs_01                                     $
                 +(    xd)*z_vec#coeffs_11 ) 

  dAdz=1.0d/dz*(  (1.0-xd)*(1.0-yd)*zp_vec#coeffs_00                           $
                 +(    xd)*(1.0-yd)*zp_vec#coeffs_10                           $
                 +(1.0-xd)*(    yd)*zp_vec#coeffs_01                           $
                 +(    xd)*(    yd)*zp_vec#coeffs_11 )


  res={f:Aval,x:dAdx, y:dAdy, z:dAdz}

  IF KEYWORD_SET(debug) THEN STOP
  RETURN, res
END
