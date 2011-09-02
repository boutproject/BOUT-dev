function filter_dct2d, fsig, a, nmax=nmax
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

  s = SIZE(fsig, /DIMENSION)

  nx=s[0]
  ny=s[1]

  if not keyword_set(a) then begin
    a=[1d0,1d0]
  endif else begin
    if size(a,/n_elements) EQ 1 then a=[a,a]
  endelse

  if not keyword_set(nmax) then nmax=[nx,ny]
  if size(nmax,/n_elements) EQ 1 then nmax=[nmax,nmax]
  if nmax[0] GT nx then nmax[0]=nx
  if nmax[1] GT ny then nmax[1]=ny

  kx = dindgen(nx)/double(nmax[0])
  mask = where(kx GT 1d0)
  fx = sinc(!dpi*kx)*sinc(!dpi*kx/a[0])
  if mask(0) NE -1 then fx(mask) = 0d0

  ky = dindgen(ny)/double(nmax[1])
  fy = sinc(!dpi*ky)*sinc(!dpi*ky/a[1])
  mask = where(ky GT 1d0)
  if mask(0) NE -1 then fy(mask) = 0d0

  ans = fsig * (fx # fy)

  return, ans
end
