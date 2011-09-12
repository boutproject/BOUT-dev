function deriv_fct1d, fk, d, chop=chop
;
; Author:       Ilon Joseph
; Began:        2011/06/22
; Modified:     2011/06/24 
; 
; Generalized derivative of order a in fourier space
;
; Notes:
;       * ans is double if fk is double
;       * returns a complex ans for fct inversion
; Options: 
;       * chop filters the result

  n = size(fk,/n_elements)
  t = size(fk,/type)
  if t EQ size(double(0),/type) or t EQ size(dcomplex(0,0),/type) then begin
        j=dcomplex(0,1)
        k = dindgen(n)/2.
  endif else begin
        j=complex(0,1)
        k = findgen(n)/2.
  endelse
;  print, k

  ans = fk*(-j*k)^d
  if keyword_set(chop) then ans = chop(ans)

  return, ans
end  
