; Bicubic interpolation with FFT in z for Apar
; Input is x,y,z coordinates
; Output is structure ={x:dA/dx,y:dA/dy,z:dA,dz}
FUNCTION apar_bicube_fft, x, y, z, debug=debug, shift=shift

  COMMON BDATA, bd
  COMMON griddata, g, deltaZtor, Ntor
  COMMON bicube_data, ad
  
  yMin=0.0
  yMax=2*!PI
  xMin=MIN(g.psixy)
  xMax=MAX(g.psixy)

;;set the maximum Fourier component to include
  Nmax=bd.nz/2.0-1.0
;  Nmax=3


;Find the lower left corner of the cell
  xind = INTERPOL(FINDGEN(g.nx),g.psixy[*,0],x)
  ; must be 1,...,g.nx-3; otherwise cubic interpolation is outside domain
  ; should return with an error message if x is larger or smaller than that
  IF ( (xind LT 1) OR (xind GE (g.nx-2)) ) THEN STOP 
  ix_int = FIX(xind)

  ; assumes that y is between 0 and 2*!PI
  yind = DOUBLE(g.ny)*(y-yMin)/(yMax-yMin)+1 ;gives floats in [1,g.ny+1)
  iy_int = FIX(yind)                         ;NOTE if y=2*!PI exactly, then we
  IF iy_int EQ (g.ny+1) THEN STOP            ;need to use twist-shift!

  ;normalized cell coordinates
  dx=g.psixy[ix_int+1,0]-g.psixy[ix_int,0]
  x_0=g.psixy[ix_int,0]
  xd=(x-x_0)/dx                       ;normalized x 
  dy=(yMax-yMin)/DOUBLE(g.ny)
  y_0=(iy_int-1)*dy
  yd=(y-y_0)/dy                       ;normalized y

  ;form x_y_vec [x^i*y^j] and its derivatives
  x_y_vec= [ [xd^0*yd^0],[xd^1*yd^0],[xd^2*yd^0],[xd^3*yd^0] ,                 $
             [xd^0*yd^1],[xd^1*yd^1],[xd^2*yd^1],[xd^3*yd^1] ,                 $
             [xd^0*yd^2],[xd^1*yd^2],[xd^2*yd^2],[xd^3*yd^2] ,                 $
             [xd^0*yd^3],[xd^1*yd^3],[xd^2*yd^3],[xd^3*yd^3] ]

  dx_y_vec=[ [0],[xd^0*yd^0],[2*xd^1*yd^0],[3*xd^2*yd^0] ,                     $
             [0],[xd^0*yd^1],[2*xd^1*yd^1],[3*xd^2*yd^1] ,                     $
             [0],[xd^0*yd^2],[2*xd^1*yd^2],[3*xd^2*yd^2] ,                     $
             [0],[xd^0*yd^3],[2*xd^1*yd^3],[3*xd^2*yd^3] ]

  x_dy_vec=[ [    0    ],[    0    ],[    0    ],[    0    ]         ,         $
             [xd^0*1*yd^0],[xd^1*1*yd^0],[xd^2*1*yd^0],[xd^3*1*yd^0] ,         $
             [xd^0*2*yd^1],[xd^1*2*yd^1],[xd^2*2*yd^1],[xd^3*2*yd^1] ,         $
             [xd^0*3*yd^2],[xd^1*3*yd^2],[xd^2*3*yd^2],[xd^3*3*yd^2] ]

;; Now for each z mode, calculate the bicubic interpolation
  Aval=0.0d
  dAdx=0.0d
  dAdy=0.0d
  dAdz=0.0d
  zper=ROUND(2*!DPI/deltaZtor)

;;check to see if n=0 coefficients need to be calculated
  IF ad[ix_int,iy_int-1,0].TF EQ 0 THEN BEGIN
    ;calculate cell-normalized coordinates for interpolation scheme
    x0d=0.0                             ;normalized x_0
    x1d=1.0                             ;normalized x_1
    xmd=(g.psixy[ix_int-1,0]-x_0)/dx    ;normalized x_-1
    x2d=(g.psixy[ix_int+2,0]-x_0)/dx    ;normalized x_2
    xd_vec=[xmd,x0d,x1d,x2d]
    y0d=0.0                             ;normalized y_0
    y1d=1.0                             ;normalized y_1
    ymd=-1.0                            ;normalized y_-1
    y2d=2.0                             ;normalized y_2
    yd_vec=[ymd,y0d,y1d,y2d]
    ;calculate coefficients
    ad[ix_int,iy_int-1,0].TF=1
    ad[ix_int,iy_int-1,0].coeffs=interp_bicube(bd.apar[(ix_int-1):(ix_int+2),  $
                                       (iy_int-1):(iy_int+2),0],xd_vec,yd_vec)
  ENDIF

;;Include n=0 component
  Aval=Aval+x_y_vec#ad[ix_int,iy_int-1,0].coeffs
  dAdx=dAdx+(dx_y_vec#ad[ix_int,iy_int-1,0].coeffs)
  dAdy=dAdy+(x_dy_vec#ad[ix_int,iy_int-1,0].coeffs)

;;Now the n=1,..n=nz/2-1 components
  nz=bd.nz
  FOR kz=1,Nmax DO BEGIN
    ;cosine coefficients are stored in 1,nz/2-1 ascending
    IF ad[ix_int,iy_int-1,kz].TF EQ 0 THEN BEGIN
      ad[ix_int,iy_int-1,kz].TF=1
      ad[ix_int,iy_int-1,kz].coeffs=interp_bicube(                             $
         bd.apar[(ix_int-1):(ix_int+2),(iy_int-1):(iy_int+2),kz],xd_vec,yd_vec)
    ENDIF
    ;sine coefficients are stored from nz-1,...,nz/2+1 descending
    IF ad[ix_int,iy_int-1,nz-kz].TF EQ 0 THEN BEGIN
      ad[ix_int,iy_int-1,nz-kz].TF=1
      ad[ix_int,iy_int-1,nz-kz].coeffs=interp_bicube(                          $
      bd.apar[(ix_int-1):(ix_int+2),(iy_int-1):(iy_int+2),nz-kz],xd_vec,yd_vec)
    ENDIF
    coeffs_ak=ad[ix_int,iy_int-1,kz].coeffs
    coeffs_bk=ad[ix_int,iy_int-1,nz-kz].coeffs

    Aval=Aval+(x_y_vec#coeffs_ak)*COS(kz*zper*z)                               $
             +(x_y_vec#coeffs_bk)*SIN(kz*zper*z)
    dAdx=dAdx+( (dx_y_vec#coeffs_ak)*COS(kz*zper*z)                            $
               +(dx_y_vec#coeffs_bk)*SIN(kz*zper*z) )
    dAdy=dAdy+( (x_dy_vec#coeffs_ak)*COS(kz*zper*z)                            $
               +(x_dy_vec#coeffs_bk)*SIN(kz*zper*z) )
    dAdz=dAdz   -1.0d*( ( x_y_vec#coeffs_ak)*kz*zper*SIN(kz*zper*z) )          $
                     +( ( x_y_vec#coeffs_bk)*kz*zper*COS(kz*zper*z) )
  ENDFOR

  dAdx=dAdx/dx
  dAdy=dAdy/dy

  res={f:Aval, x:dAdx, y:dAdy, z:dAdz}
;
;
;
  IF KEYWORD_SET(debug) THEN STOP
  RETURN, res
END
