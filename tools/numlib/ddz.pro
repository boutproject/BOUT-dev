; Calculate derivative in z (toroidal angle)
; using Fourier method
; for dimensions less than 3, assumes that the final index represents z
; for 4 dimensions, assumes that the 3rd dimension is z (and the final is t)
;
; Created:  2012/03/19  I. Joseph
; Modified: 2012/03/27  I. Joseph
;

FUNCTION ddz, mesh, f

  ndim = size(f,/n_dimensions)
  dim = size(f,/dimensions)
  df = 0*f 
  if ndim lt 4 then begin 
        nz = dim(ndim-1)
  endif else begin
        nz = dim(2)
  endelse
  dz = 2.*!PI / FLOAT(nz)

  case ndim of
  1: begin
        df=fft_deriv(f)
     end
  2: begin
        for i=0,dim(0)-1 do begin
          df[i,*] = fft_deriv(reform(f[i,*]))
        endfor
     end
  3: begin
        for i=0,dim(0)-1 do begin
        for j=0,dim(1)-1 do begin
          df[i,j,*] = fft_deriv(reform(f[i,j,*]))
        endfor
        endfor
     end
  4: begin
        for i=0,dim(0)-1 do begin
        for j=0,dim(1)-1 do begin
        for k=0,dim(3)-1 do begin
          df[i,j,*,k] = fft_deriv(reform(f[i,j,*,k]))
        endfor
        endfor
        endfor
     end
  endcase
  
  ftype = size(f, /tname)
  if ftype ne 'COMPLEX' OR ftype ne 'DCOMPLEX' then df=real_part(df)
  return, df/dz 
END
