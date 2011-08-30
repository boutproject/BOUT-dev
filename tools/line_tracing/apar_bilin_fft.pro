; Bilinear interpolation with FFT in z for Apar
; Input is x,y,z coordinates
; Output is structure ={x:dA/dx,y:dA/dy,z:dA,dz}
FUNCTION apar_bilin_fft, x, y, z, debug=debug, shift=shift

  COMMON BDATA, bd ;-bd={apar:apar, nx:nx, ny:ny, nz:nz}
  COMMON griddata, g, deltaZtor, Ntor
  
  yMin=0.0
  yMax=2*!PI
  xMin=MIN(g.psixy)
  xMax=MAX(g.psixy)

;Find the lower left corner of the cell
  xind = INTERPOL(FINDGEN(g.nx),g.psixy[*,0],x)
  ; must be 0,...,g.nx-2; if ix_int>=g.nx-1 then point is outside domain
  ; should return with an error message if x is larger or smaller than that
  IF ((xind<0) OR (xind GE g.nx-1)) THEN STOP 
  ix_int = FIX(xind)

  ; assumes that y is between 0 and 2*!PI
  yind = DOUBLE(g.ny)*(y-yMin)/(yMax-yMin) ;gives floats between 0 and g.ny
  iy_int = FIX(yind)                       ;NOTE if y=2*!PI exactly, then we 
  IF iy_int EQ g.ny THEN STOP              ;need to use twist-shift!

;Calculate grid spacing
  dx=g.psixy[ix_int+1,0]-g.psixy[ix_int,0]
  dy=(yMax-yMin)/DOUBLE(g.ny)

;normalized cell coordinates
  x_0=g.psixy[ix_int,0]
  xd=(x-x_0)/dx    ;normalized x 
  y_0=iy_int*dy
  yd=(y-y_0)/dy    ;normalized y

;; Now for each z mode, calculate the bilinear interpolation
  Aval=0.0d
  dAdx=0.0d
  dAdy=0.0d
  dAdz=0.0d

  zper=ROUND(2*!DPI/deltaZtor)

;;Include n=0 component
  c_corners=bd.Apar[ix_int:(ix_int+1),iy_int:(iy_int+1),0]

  Aval=Aval+( c_corners[0,0]*(1.0-xd)*(1.0-yd)                                 $
              c_corners[1,0]*(    xd)*(1.0-yd)                                 $
              c_corners[0,1]*(1.0-xd)*(    yd)                                 $
              c_corners[1,1]*(    xd)*(    yd) )

  dAdx=dAdx+1.0d/dx*(-c_corners[0,0]*(1.0-yd)                                  $
                     +c_corners[1,0]*(1.0-yd)                                  $
                     -c_corners[0,1]*     yd                                   $
                     +c_corners[1,1]*     yd  )

  dAdy=dAdy+1.0d/dy*(-c_corners[0,0]*(1.0-xd)                                  $
                     -c_corners[1,0]*     xd                                   $
                     +c_corners[0,1]*(1.0-xd)                                  $
                     +c_corners[1,1]*     xd  )

;;Now the n=1,..n=nz/2-1 components
  FOR kz=1,(bd.nz/2)-1 DO BEGIN
    ;cosine coefficients are stored in 1,nz/2-1 ascending
    ;sine coefficients are stored from nz-1,...,nz/2+1 descending
    c_corners=bd.Apar[ix_int:(ix_int+1),iy_int:(iy_int+1),kz]
    s_corners=bd.Apar[ix_int:(ix_int+1),iy_int:(iy_int+1),bd.nz-kz]
;;;;Using bilinear interpolation
    Aval=Aval+( c_corners[0,0]*(1.0-xd)*(1.0-yd)                               $
                c_corners[1,0]*(    xd)*(1.0-yd)                               $
                c_corners[0,1]*(1.0-xd)*(    yd)                               $
                c_corners[1,1]*(    xd)*(    yd) )*COS(kz*zper*z)              $
              ( s_corners[0,0]*(1.0-xd)*(1.0-yd)                               $
                s_corners[1,0]*(    xd)*(1.0-yd)                               $
                s_corners[0,1]*(1.0-xd)*(    yd)                               $
                s_corners[1,1]*(    xd)*(    yd) )*SIN(kz*zper*z)

    dAdx=dAdx+1.0d/dx*(-c_corners[0,0]*(1.0-yd)                                $
                       +c_corners[1,0]*(1.0-yd)                                $
                       -c_corners[0,1]*     yd                                 $
                       +c_corners[1,1]*     yd  )*COS(kz*zper*z)               $
             +1.0d/dx*(-s_corners[0,0]*(1.0-yd)                                $
                       +s_corners[1,0]*(1.0-yd)                                $
                       -s_corners[0,1]*     yd                                 $
                       +s_corners[1,1]*     yd  )*SIN(kz*zper*z)

    dAdy=dAdy+1.0d/dy*(-c_corners[0,0]*(1.0-xd)                                $
                       -c_corners[1,0]*     xd                                 $
                       +c_corners[0,1]*(1.0-xd)                                $
                       +c_corners[1,1]*     xd  )*COS(kz*zper*z)               $
             +1.0d/dy*(-s_corners[0,0]*(1.0-xd)                                $
                       -s_corners[1,0]*     xd                                 $
                       +s_corners[0,1]*(1.0-xd)                                $
                       +s_corners[1,1]*     xd  )*SIN(kz*zper*z)

    dAdz=dAdz   -1.0*( c_corners[0,0]*(1.0-xd)*(1.0-yd)                        $
                      +c_corners[1,0]*     xd *(1.0-yd)                        $
                      +c_corners[0,1]*(1.0-xd)*     yd                         $
                      +c_corners[1,1]*     xd *     yd )*kz*zper*SIN(kz*zper*z)$
                    +( s_corners[0,0]*(1.0-xd)*(1.0-yd)                        $
                      +s_corners[1,0]*     xd *(1.0-yd)                        $
                      +s_corners[0,1]*(1.0-xd)*     yd                         $
                      +s_corners[1,1]*     xd *     yd )*kz*zper*COS(kz*zper*z)
  ENDFOR

  res={f:Aval, x:dAdx, y:dAdy, z:dAdz}

  IF KEYWORD_SET(debug) THEN STOP
  RETURN, res
END
