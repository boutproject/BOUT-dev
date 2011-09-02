function deriv_fft1d, fk, d, chop=chop, norm = norm
;
; Author:       Ilon Joseph
; Began:        2011/06/22
; Modified:     2011/06/24 
; 
; Generalized derivative of order d in fourier space
; Notes:
;       * ans is double if fk is double
; Options: 
;       * chop filters the result

  n = size(fk,/n_elements)
  t = size(fk,/type)
  if t EQ size(double(0),/type) or t EQ size(dcomplex(0,0),/type) then begin
        j = dcomplex(0,1)
        k = dindgen(n)
  endif else begin
        j = complex(0,1)
        k = findgen(n)
  endelse

  k(n/2:n-1) = k(n/2:n-1)-n 
  if (n mod 2) EQ 0 then k(n/2)=0.0
;  print, k

  ans = fk*(j*k)^d
  if keyword_set(chop) then ans = chop(ans)

  return, ans
end  
