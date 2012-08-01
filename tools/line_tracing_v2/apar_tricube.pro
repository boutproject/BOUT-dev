; Tricubic interpolation for Apar
; Input is x,y,z coordinates
; Output is structure ={f:A(x,y,z),x:dA/dx,y:dA/dy,z:dA,dz}
FUNCTION apar_tricube, x, y, z, debug=debug

  COMMON BDATA, bd
  COMMON griddata, g, deltaZtor, Ntor
  COMMON tricube_data, A_int

  xMin=MIN(g.psixy)
  xMax=MAX(g.psixy) 
  yMin=0.0
  yMax=2*!PI
  zMin=0.0
  zMax=deltaZtor
  nz=bd.nz

;Find the lower left corner of the cell
  xind = INTERPOL(FINDGEN(g.nx),g.psixy[*,0],x)
  ; must be 1,...,g.nx-3; otherwise cubic interpolation is outside domain
  ; should return with an error message if x is larger or smaller than that
  IF ( (xind LT 1) OR (xind GE (g.nx-2)) ) THEN STOP, 'X out of bounds!' 
  ix_int = FIX(xind)

  ; assumes that y is between 0 and 2*!PI
  yind = DOUBLE(g.ny)*(y-yMin)/(yMax-yMin)   ;gives floats in [0,g.ny)
  iy_int = FIX(yind)
  IF iy_int EQ (g.ny) THEN STOP, 'Y out of bounds!'

  ; Do not assume that z is between 0 and deltaZtor
  zmod = FMODULO(z,zMax)
  zind = DOUBLE(nz)*(zmod-zMin)/(zMax-zMin)
  iz_int = FIX(zind)

  ;find normalized cell coordinates
  dx=g.psixy[ix_int+1,0]-g.psixy[ix_int,0]
  x_0=g.psixy[ix_int,0]
  xd=(x-x_0)/dx                       ;normalized x 
  dy=(yMax-yMin)/DOUBLE(g.ny)
  y_0=iy_int*dy
  yd=(y-y_0)/dy                       ;normalized y
  dz=(zMax-zMin)/DOUBLE(nz)
  z_0=iz_int*dz 
  zd=(zmod-z_0)/dz

  ;form x_y_z_vec [x^i*y^j*z^k] and its derivatives
  x_y_z  = make_tricube_vec(xd,yd,zd,derivative='f')
  dx_y_z = make_tricube_vec(xd,yd,zd,derivative='dx')
  x_dy_z = make_tricube_vec(xd,yd,zd,derivative='dy')
  x_y_dz = make_tricube_vec(xd,yd,zd,derivative='dz')
 
;; Now check to see if the tricubic interpolant has been calculated already
  IF A_int[ix_int,iy_int,iz_int].TF EQ 0 THEN BEGIN
    ;calculate normalized cell coordinates for interpolation scheme
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
    z0d=0.0                             ;normalized z_0
    z1d=1.0                             ;normalized z_1
    zmd=-1.0                            ;normalized z_-1
    z2d=2.0                             ;normalized z_2
    zd_vec=[zmd,z0d,z1d,z2d]

    ;assemble the function values matrix using appropriate +/- indices for y,z
    vals=FLTARR(4,4,4)
    FOR i_=0,3 DO BEGIN
      FOR j_=0,3 DO BEGIN
        FOR k_=0,3 DO BEGIN
          ;i_, j_, k_ are indices for vals matrix
          ;ix, jy, kz are indices for bd.apar matrix
          ;ix_int,iy_int,iz_int are fixed global indices that identify the cell
          ix = ix_int+i_-1
          jy = iy_int+j_-1
          kz = iz_int+k_-1
          ;kz index can be recast into appropriate domain via fmodulo function
          kz=ROUND(FMODULO(DOUBLE(kz),DOUBLE(nz)))
          ;; if the y index is within the domain simply use the value
          IF ( ( jy GE 0) AND ( jy LE g.ny-1) ) THEN BEGIN
            vals[i_,j_,k_]=bd.apar[ix,jy,kz]
          ENDIF ELSE BEGIN
            twist=g.shiftangle[ix]
            zval=kz*dz
            zper=ROUND(2*!DPI/zMax)
            CASE jy OF
            ;;Case 1: jy = -1, take FFT of data at y=2*PI-dy
              -1 : BEGIN
                     fft_=FFT(REFORM(bd.apar[ix,g.ny-1,*]))
                     val=REAL_PART(fft_[0])   ; n=0 comp.
                     FOR fkz=1,(nz/2)-1 DO BEGIN
                       val=val                                                 $
                          +2.*REAL_PART(fft_[fkz])*COS(fkz*zper*(zval+twist))  $
                          -2.*IMAGINARY(fft_[fkz])*SIN(fkz*zper*(zval+twist))
                     ENDFOR
                     vals[i_,j_,k_]=val
                   END
            ;;Case 2: jy = g.ny, take FFT of data at y=0 
            g.ny : BEGIN
                     fft_=FFT(REFORM(bd.apar[ix,0,*]))
                     val=REAL_PART(fft_[0])   ; n=0 comp.
                     FOR fkz=1,(nz/2)-1 DO BEGIN
                       val=val                                                 $
                          +2.*REAL_PART(fft_[fkz])*COS(fkz*zper*(zval-twist))  $
                          -2.*IMAGINARY(fft_[fkz])*SIN(fkz*zper*(zval-twist))
                     ENDFOR
                     vals[i_,j_,k_]=val
                   END
            ;;Case 3: jy = g.ny+1, take FFT of data at y=0+dy 
           g.ny+1: BEGIN
                     fft_n=FFT(REFORM(bd.apar[ix,1,*]))
                     val=REAL_PART(fft_n[0])   ; n=0 comp.
                     FOR fkz=1,(nz/2)-1 DO BEGIN
                       val=val                                                 $
                          +2.*REAL_PART(fft_[fkz])*COS(fkz*zper*(zval-twist))  $
                          -2.*IMAGINARY(fft_[fkz])*SIN(fkz*zper*(zval-twist))
                     ENDFOR
                     vals[i_,j_,k_]=val
                   END
             else: STOP,'Error in determining tricubic interpolant values'
            ENDCASE        
          ENDELSE
        ENDFOR
      ENDFOR
    ENDFOR
    ;calculate coefficients
    A_int[ix_int,iy_int,iz_int].TF=1
    A_int[ix_int,iy_int,iz_int].coeffs=interp_tricube(vals,xd_vec,yd_vec,zd_vec)
  ENDIF

  Aval=         TRANSPOSE(x_y_z )#A_int[ix_int,iy_int,iz_int].coeffs
  dAdx=1.0d/dx*(TRANSPOSE(dx_y_z)#A_int[ix_int,iy_int,iz_int].coeffs)
  dAdy=1.0d/dy*(TRANSPOSE(x_dy_z)#A_int[ix_int,iy_int,iz_int].coeffs)
  dAdz=1.0d/dz*(TRANSPOSE(x_y_dz)#A_int[ix_int,iy_int,iz_int].coeffs)

  res={f:Aval, x:dAdx, y:dAdy, z:dAdz}
;
  IF KEYWORD_SET(debug) THEN STOP
  RETURN, res
END
