; Calculate derivative in x (psi)
; uses IDL 3-point DERIV function
; 
; Assumes 1st index of f corresponds to the x coordinate
; mesh can either be a bout mesh structure or 
; a 1d or 2d array of x values
;
; Created:  2012/03/19  I. Joseph
; Modified: 2012/03/19  I. Joseph
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

FUNCTION intx, mesh, f, normalized=normalized

  if size(mesh,/tname) eq 'STRUCT' then begin
    if not keyword_set(normalized) then begin
      x =  mesh.psixy
    endif else begin  
      x = (mesh.psixy - mesh.psi_axis)/(mesh.psi_bndry-mesh.psi_axis)
    endelse
  endif else begin  
     x = mesh
  endelse

  ndim = size(f,/n_dimensions)
  dim = size(f,/dimensions)
  fint = 0*f 

  case ndim of
  1: begin
        fint=int0(x,f)
     end
  2: begin
        for i=0,dim(1)-1 do begin
          fint[*,i] = int0(reform(x[*,i]),reform(f[*,i]))
        endfor
     end
  3: begin
        for i=0,dim(1)-1 do begin
        for j=0,dim(2)-1 do begin
          fint[*,i,j] = int0(reform(x[*,i]),reform(f[*,i,j]))
        endfor
        endfor
     end
  4: begin
        for i=0,dim(1)-1 do begin
        for j=0,dim(2)-1 do begin
        for k=0,dim(3)-1 do begin
          fint[*,i,j,k] = int0(reform(x[*,i]),reform(f[*,i,j,k]))
        endfor
        endfor
        endfor
     end
  endcase

  return, fint
END
