function filter_fft1d, fk, a, filter_type=filter_type, nmax=nmax
;
; Author:       Ilon Joseph
; Began:        2011/06/22
; Modified:     2011/06/24 
; 
; Filters for smoothing Gibbs oscillations in FFT reconstruction
; answer in double precision
;
; Notes:
;       * returns double
;
; Options:
;       * a = parameter for filter
;       * filter_type = sets the type of filter
;                       currently implemented: lanczos, cesaro

  n = size(fk,/n_elements)
  k = dindgen(n) 
  k(n/2:n-1) = k(n/2:n-1)-n 
;  print, k
  if not keyword_set(nmax) then nmax=n/2
  if nmax GT n/2 then nmax=n/2
  k /= double(nmax)
  mask = where(abs(k) GT nmax)
  if mask(0) NE -1 then k(mask) = 1d0

  if not keyword_set(filter_type) then begin
    filter_type = "lanczos"
  endif

  case filter_type of
  "lanczos": begin
;     Lanczos filter fk*sinc(!dpi*k)*sinc(!dpi*k/a)
;     a=1 or 2 usually work best
      if not keyword_set(a) then a=1d0
      ans = fk * sinc(!dpi*k)*sinc(!dpi*k/a)
    end
  "cesaro": begin
;     Cesaro filter fk*(1-k)^a 
;     choose parameter a based on power law
      if not keyword_set(a) then a=1d0
      ans = fk * (1-abs(k))^a      
    end
  endcase 

  return, ans
end

