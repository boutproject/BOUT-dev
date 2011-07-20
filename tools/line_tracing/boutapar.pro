function boutApar, x, y, z, debug=debug, shift=shift
;
; Input: x,y,z coordinates
; Output: {dApar/d(x), dApar/d(y), dApar/d(z)}
;====================================================;

  COMMON BDATA, bd ;-bd={apar:apar, nx:nx, ny:ny, nz:nz}
  COMMON griddata, g, deltaZtor, Ntor
  
  yMin=0.0
  yMax=2*!PI
  xMin=g.psixy[0,0]
  xMax=g.psixy[g.nx-1,g.ny-1]
  zmin=0.0
  zmax=deltaZtor                ;;2*!PI/period
  nz=bd.nz
  
  ;-use closest grid point for now
  xind = INTERPOL(FINDGEN(g.nx), g.psixy[*,0], x)
  ix_int = FIX( xind )
  ix_int = (ix_int > 1) < (g.nx-2) ; Keep within range
  
  yind = (g.ny-1)*(y-yMin)/(yMax-yMin)
  iy_int=FIX(yind)

  dz = (zMax-zMin) / FLOAT(nz)
  zind = (z-zMin)/dz
  
  ;-grid spacing
  dx=(g.psixy[ix_int+1,0]-g.psixy[ix_int-1,0])/2.0
  dy=2*!PI/g.ny

  IF KEYWORD_SET(shift) THEN BEGIN
    ; Use derivatives in psi rather than x

    xm = ix_int
    xp = ix_int+1
    
    ym = iy_int
    yp = (iy_int + 1) MOD g.ny
    
    ; Calculate derivatives at ym
    
    ; Get the Z shifts
    xm_zs = g.zshift[xm,ym]
    xp_zs = g.zshift[xp,ym]
    zs = INTERPOL(g.zshift[*,ym], g.psixy[*,0], x)
    
    ; Get Z indices
    xm_zi = zind + (zs - xm_zs)/dz
    xp_zi = zind + (zs - xp_zs)/dz
    
    A = FLTARR(2,2) ; [xm,xp][ym,yp]
    dAdZ = A
    
    ; Get Apar, dA/dZ at (xm, ym, xm_zi) and (xp, ym, xp_zi)
    
    A[0,0]    = interpZApar(xm,ym,xm_zi)
    dAdZ[0,0] = derivZApar(xm,ym,xm_zi)
    A[1,0]    = interpZApar(xp,ym,xp_zi)
    dAdZ[1,0] = derivZApar(xp,ym,xp_zi)
    
    ; Calculate derivatives at yp
    
    xm_zs = g.zshift[xm,yp]
    xp_zs = g.zshift[xp,yp]
    zs = INTERPOL(g.zshift[*,yp], g.psixy[*,0], x)
    
    xm_zi = zind + (zs - xm_zs)/dz
    xp_zi = zind + (zs - xp_zs)/dz
    
    IF yp EQ 0 THEN BEGIN
      ; Crossing the twist-shift location
      shift = INTERPOL(g.Shiftangle, g.psixy[*,0], x)
      xm_zi = zind + shift / dz
      xp_zi = zind + shift / dz
    ENDIF
    
    ; Get Apar, dA/dZ at (xm, yp, xm_zi) and (xp, yp, xp_zi)
    
    A[0,1]    = interpZApar(xm,yp,xm_zi)
    dAdZ[0,1] = derivZApar(xm,yp,xm_zi)
    A[1,1]    = interpZApar(xp,yp,xp_zi)
    dAdZ[1,1] = derivZApar(xp,yp,xp_zi)
  
  
    ; Calculate A, dA/dpsi and dA/dy using 4 corners
    IF (ABS(xind-xm) GT 1.) OR (ABS(yind-ym) GT 1.) THEN STOP
    pd = pdiff_xy(A, xind-xm, yind-ym)
    
    dAdx = pd.dfdx/dx
    dAdy = pd.dfdy/dy
    
    ; Interpolate dAdz
    pd = pdiff_xy(dAdZ, xind-xm, yind-ym)
    dAdz = pd.f
    
  ENDIF ELSE BEGIN
    iz_int=FIX(zind) MOD nz;;-need closed-peridoic domain!!!
    
    ;-calculate partial derivatives within cubic cell
    cellVal=bd.Apar[ix_int:(ix_int+1),iy_int:(iy_int+1),iz_int:(iz_int+1)]

    ;-derivatives with respect to grid index
    dAparInd=CubeDeriv(cellVal)
    
    ;-derivatives with respect to x,y,z
    dAdx=dAparInd.x/dx
    dAdy=dAparInd.y/dy
    dAdz=dAparInd.z/dz
  ENDELSE
  
  res={x:dAdx, y:dAdy, z:dAdz}
;
;
;
if keyword_set(DEBUG) then STOP
return, res
end

