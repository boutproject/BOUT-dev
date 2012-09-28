; Calculate derivative in x (psi)
; uses IDL 3-point DERIV function
; 
; Assumes 1st index of f corresponds to the x coordinate
; mesh can either be a bout mesh structure or 
; a 1d or 2d array of x values
;
; Created:  2012/03/19  I. Joseph
; Modified: 2012/03/27  I. Joseph
;

FUNCTION xderiv, x,f
 RETURN, DERIV(x,f)
END

FUNCTION psinorm, mesh
  RETURN, (mesh.psixy - mesh.psi_axis)/(mesh.psi_bndry-mesh.psi_axis)
END

FUNCTION ddx, mesh, f, normalized=normalized

  if size(mesh,/tname) eq 'STRUCT' then begin
    if not keyword_set(normalized) then begin
      x = mesh.psixy
    endif else begin  
      x = psinorm(mesh)
    endelse
  endif else begin  
     x = mesh
  endelse

  ndim = size(f,/n_dimensions)
  dim = size(f,/dimensions)
  df = 0*f 

  case ndim of
  1: begin
        if size(x,/dimensions) eq 2 then x=x[*,0]
        df=xderiv(x,f)
     end
  2: begin
        for i=0,dim(1)-1 do begin
          df[*,i] = xderiv(reform(x[*,i]),reform(f[*,i]))
        endfor
     end
  3: begin
        for i=0,dim(1)-1 do begin
        for j=0,dim(2)-1 do begin
          df[*,i,j] = xderiv(reform(x[*,i]),reform(f[*,i,j]))
        endfor
        endfor
     end
  4: begin
        for i=0,dim(1)-1 do begin
        for j=0,dim(2)-1 do begin
        for k=0,dim(3)-1 do begin
          df[*,i,j,k] = xderiv(reform(x[*,i]),reform(f[*,i,j,k]))
        endfor
        endfor
        endfor
     end
  endcase

  return, df
END
