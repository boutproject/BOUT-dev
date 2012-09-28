; Calculate derivative in theta
;
; Created:  2012/03/19  I. Joseph
; Modified: 2012/03/27  I. Joseph
;

FUNCTION multiply32, a, b
  adim = size(a,/n_dimensions)
  bdim = size(b,/n_dimensions)
  if (adim eq 3) and (bdim eq 2) then begin
        dim = size(a,/dimensions)
        c3d = fltarr(dim)
        FOR k=0, dim(2)-1 DO c3d[*,*,k] = b*a[*,*,k]
  endif else begin
        dim = size(b,/dimensions)
        c3d = fltarr(dim)
        FOR k=0, dim(2)-1 DO c3d[*,*,k] = a*b[*,*,k]
  endelse

  RETURN, c3d
END

FUNCTION ddtheta, mesh, f

  Qlocal = mesh.Btxy*mesh.Hthe/(mesh.Bpxy*mesh.Rxy)
  dfdtheta = ddy(mesh,f) + multiply32(ddz(mesh,f),Qlocal)
  return, dfdtheta
END
