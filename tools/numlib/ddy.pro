; Calculate derivative in y (poloidal angle)
; taking into account branch-cuts
; assumes that the 2nd index represents y

; Take derivative in y, taking into account branch-cuts
;
; Created:  2012/03/19  I. Joseph
; Modified: 2012/03/27  I. Joseph
;

FUNCTION yderiv, mesh, f 

  df = 0*f


  status = gen_surface(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=period)
    ny = size(yi,/N_ELEMENTS)
    dtheta = 2.*!PI / FLOAT(ny)
    IF period THEN BEGIN
       df[xi,yi] = fft_deriv(f[xi,yi])/dtheta
    ENDIF ELSE BEGIN
       df[xi,yi] = DERIV(f[xi,yi])/dtheta
    ENDELSE
  ENDREP UNTIL last

  ftype = size(f, /tname)
  if ftype ne 'COMPLEX' OR ftype ne 'DCOMPLEX' then df=real_part(df)
 
  RETURN, df 
END

FUNCTION ddy, mesh, f

  ndim = size(f,/n_dimensions)
  dim = size(f,/dimensions)
  df = 0*f 

  case ndim of
  2: begin
        df=yderiv(mesh,f)
     end
  3: begin
        for iz=0,dim(2)-1 do begin
          df[*,*,iz] = yderiv(mesh,reform(f[*,*,iz]))
        endfor
     end
  4: begin
        for iz=0,dim(2)-1 do begin
        for it=0,dim(3)-1 do begin
          df[*,*,iz,it] = yderiv(mesh,reform(f[*,*,iz,it]))
        endfor
        endfor
     end
  endcase

  return, df
END
