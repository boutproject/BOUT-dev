;Given 64 function values f(i,j,k) for i,j,k=-1,0,1,2
;and normalized coordinates x(i),y(j),z(k) calculate
;coefficients of P(x)=SUM(i=0,3)SUM(j=0,3)SUM(k=0,3)a_ijk*x^i*y^j*z^k
FUNCTION interp_tricube, f, x, y, z

  COMMON tricube_data, A_int, A_inv
  ;;A_inv is stored in this common data so it only needs to be calculated once

;;First, estimate the derivatives numerically
  f_   =FLTARR(8)
  fx_  =FLTARR(8)
  fy_  =FLTARR(8)
  fz_  =FLTARR(8)
  fxy_ =FLTARR(8)
  fxz_ =FLTARR(8)
  fyz_ =FLTARR(8)
  fxyz_=FLTARR(8)
  FOR k=1,2 DO BEGIN
    FOR j=1,2 DO BEGIN
      FOR i=1,2 DO BEGIN
        f_[(i-1)+(j-1)*2+(k-1)*4]=f[i,j,k]
  ; d_x( f(i,j,k) )~= ( f(i+1,j,k)-f(i-1,j,k) )/( x(i+1)-x(i-1) )
       fx_[(i-1)+(j-1)*2+(k-1)*4]=( f[i+1,j,k]-f[i-1,j,k] )/( x[i+1]-x[i-1] )
  ; d_y( f(i,j,k) )~= ( f(i,j+1,k)-f(i,j-1,k) )/( y(j+1)-y(j-1) )
        fy_[(i-1)+(j-1)*2+(k-1)*4]=( f[i,j+1,k]-f[i,j-1,k] )/( y[j+1]-y[j-1] )
  ; d_z( f(i,j,k) )~= ( f(i,j,k+1)-f(i,j,k-1) )/( z(k+1)-z(k-1) )
        fz_[(i-1)+(j-1)*2+(k-1)*4]=( f[i,j,k+1]-f[i,j,k-1] )/( z[k+1]-z[k-1] )
  ; d_xy(f(i,j,k) )~= ( f(i+1,j+1,k)-f(i+1,j-1,k)-f(i-1,j+1,k)+f(i-1,j-1,k) )
  ;                    /( ( x(i+1)-x(i-1) )*( y(j+1)-y(j-1) ) )
        fxy_[(i-1)+(j-1)*2+(k-1)*4]=                                           $
          ( f[i+1,j+1,k]-f[i+1,j-1,k]-f[i-1,j+1,k]+f[i-1,j-1,k] )              $
           /( (x[i+1]-x[i-1])*(y[j+1]-y[j-1]) )
  ; d_xz(f(i,j,k) )~= ( f(i+1,j,k+1)-f(i+1,j,k-1)-f(i-1,j,k+1)+f(i-1,j,k-1) )
  ;                    /( ( x(i+1)-x(i-1) )*( z(k+1)-z(k-1) ) )
        fxz_[(i-1)+(j-1)*2+(k-1)*4]=                                           $
          ( f[i+1,j,k+1]-f[i+1,j,k-1]-f[i-1,j,k+1]+f[i-1,j,k-1] )              $
           /( (x[i+1]-x[i-1])*(z[k+1]-z[k-1]) )
  ; d_yz(f(i,j,k) )~= ( f(i,j+1,k+1)-f(i,j+1,k-1)-f(i,j-1,k+1)+f(i,j-1,k-1) )
  ;                    /( ( y(j+1)-y(j-1) )*( z(k+1)-z(k-1) ) )
        fyz_[(i-1)+(j-1)*2+(k-1)*4]=                                           $
          ( f[i,j+1,k+1]-f[i,j+1,k-1]-f[i,j-1,k+1]+f[i,j-1,k-1] )              $
           /( (y[j+1]-y[j-1])*(z[k+1]-z[k-1]) )
  ;d_xyz(f(i,j,k))=(f(i+1,j+1,k+1)-f(i+1,j-1,k+1)-f(i-1,j+1,k+1)+f(i-1,j-1,k+1)
  ;                -f(i+1,j+1,k-1)+f(i+1,j-1,k-1)+f(i-1,j+1,k-1)-f(i-1,j-1,k-1))
  ;                /( (x(i+1)-x(i-1))*(y(j+1)-y(j-1))*(z(k+1)-z(k-1)) )

        fxyz_[(i-1)+(j-1)*2+(k-1)*4]=                                          $
          ( f[i+1,j+1,k+1]-f[i+1,j-1,k+1]-f[i-1,j+1,k+1]+f[i-1,j-1,k+1]        $
           -f[i+1,j+1,k-1]+f[i+1,j-1,k-1]+f[i-1,j+1,k-1]-f[i-1,j-1,k-1] )      $
           /( (x[i+1]-x[i-1])*(y[j+1]-y[j-1])*(z[k+1]-z[k-1]) )
      ENDFOR
    ENDFOR
  ENDFOR
;;Next, assemble the RHS vector of function+derivative values
  f_vec=[f_,fx_,fy_,fz_,fxy_,fxz_,fyz_,fxyz_]

;;Lastly, multiply the inversion matrix by the vector of function values
  coeffs=TRANSPOSE(A_inv)#f_vec
  RETURN,coeffs
END
