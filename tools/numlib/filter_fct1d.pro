function filter_fct1d, fk, a, filter_type=filter_type, nmax=nmax
;
; Author:       Ilon Joseph
; Began:        2011/06/22
; Modified:     2011/06/24 
; 
; Filters for smoothing Gibbs oscillations in DCT reconstruction
; Notes:
;       * returns double
;
; Options:
;       * a = parameter for filter
;       * nmax = cutoff
;       * filter_type = sets the type of filter
;                       currently implemented: lanczos, cesaro


  n = size(fk,/n_elements)
  k = dindgen(n)

  if not keyword_set(filter_type) then begin
    filter_type = "lanczos"
  endif

  if not keyword_set(nmax) then nmax=n
  if nmax GT n then nmax=n
  k /= double(nmax)
  mask = where(abs(k) GT nmax)
  if mask(0) NE -1 then k(mask) = 1d0

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
      if not keyword_set(a) then a=2d0
      ans = fk * (1-k)^a      
    end
  endcase 

  return, ans
end
