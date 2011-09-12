function splder, x0, y0, y2, period=period
;
; Author:       Ilon Joseph
; Began:        2011/06/01
; Modified:     2011/06/08 
;
; derivative of cubic spline interpolation based on Numerical Recipes in C
; y2 = second derivatives = output of IDL function spl_init or splgen
;
; Options:
;       * period = period length for periodic spline
;         assumes that x0 has branch cut at boundary which jumps by 1 period



n0 = size(x0,/n_elements)

;print, n0 

if not keyword_set(period) then begin
  dx  = [x0(1:n0-1) - x0(0:n0-2), x0(n0-1)-x0(n0-2)]
  dy  = [y0(1:n0-1) - y0(0:n0-2), y0(n0-1)-y0(n0-2)]
  dy2 = [y2(1:n0-1) - y2(0:n0-2), y2(n0-1)-y2(n0-2)]
endif else begin
  dx  = [x0(1:n0-1) - x0(0:n0-2), x0(0)-x0(n0-1)+period]
  dy  = [y0(1:n0-1) - y0(0:n0-2), y0(0)-y0(n0-1)]
  dy2 = [y2(1:n0-1) - y2(0:n0-2), y2(0)-y2(n0-1)]
endelse


;print, dx
;print, dy
;print, dy2

; Compute the answer

 y1 = dy/dx  - (y2+dy2/3.0)*dx/2.0

return, y1

end
