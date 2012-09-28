; Calculate integral in y (poloidal angle)
; taking into account branch-cuts
; assumes that the 2nd index represents y
;
; Created:  2012/03/19  I. Joseph
; Modified: 2012/03/20  I. Joseph
;

FUNCTION diff, f
 n=size(f,/n_elements)
 return, f(1:n-1)-f(0:n-2)
END

FUNCTION int0, x, f
 n=size(f,/n_elements)
 dx = diff(x)
 df = diff(f)
 RETURN, [0, total((f(0:n-1)+df/2.)*dx, /cum)]
END


FUNCTION integratey, mesh, f 

  fint = 0*f
;  print, mesh.npol
  dtheta = 2.*!PI / FLOAT(TOTAL(mesh.npol))

  status = gen_surface(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=period)
    IF period THEN BEGIN
       fint[xi,yi] = fft_integrate(f[xi,yi])
    ENDIF ELSE BEGIN
       fint[xi,yi] = int0(yi,f[xi,yi])
    ENDELSE
  ENDREP UNTIL last

  RETURN, fint*dtheta
END

FUNCTION inty, mesh, f

  ndim = size(f,/n_dimensions)
  dim = size(f,/dimensions)
  fint = 0*f 

  case ndim of
  2: begin
        fint=integratey(mesh,f)
     end
  3: begin
        for iz=0,dim(2)-1 do begin
          fint[*,*,iz] = integratey(mesh,reform(f[*,*,iz]))
        endfor
     end
  4: begin
        for iz=0,dim(2)-1 do begin
        for it=0,dim(3)-1 do begin
          fint[*,*,iz,it] = integratey(mesh,reform(f[*,*,iz,it]))
        endfor
        endfor
     end
  endcase

  return, fint
END
