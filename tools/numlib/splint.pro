function splint, x0, y0, y2, xinterp, deriv=deriv, period=period
;
; Author:       Ilon Joseph
; Began:        2011/06/01
; Modified:     2011/06/07 
;
; cubic spline interpolation based on Numerical Recipies in C
; y2 = second derivatives = output of IDL function spl_init
;
; Options
;       * deriv=0 or unsupplied returns y(xinterp)
;       * deriv=1 returns dy/dx
;       * deriv=2 returns d2y/dx2
;       * period = period length for periodic spline
;         assumes that x0 has branch cut at boundary which jumps by 1 period


; Search for intervals that bracket the points to be interpolated

N0 = size(x0,/n_elements)
Nin = size(xinterp,/n_elements)
;print, N0, Nin

if not keyword_set(period) then begin
        t0=x0
        z0=y0
        z2=y2
        xin = xinterp 
endif else begin
        xin = xinterp mod period
        t0=[x0(N0-1)-period,x0 mod period,x0(0)+period]
        z0=[y0(N0-1),y0,y0(0)]
        z2=[y2(N0-1),y2,y2(0)]
        N0+=2
endelse

klo = 0*indgen(Nin)
khi = klo + N0 - 1

for i = 0, Nin-1 do begin
;  print, "i = ",i," x(i) = ",xin(i)
  while( khi(i)-klo(i) GT 1) do begin
        k = (klo(i) + khi(i))/2                ;integer arithmetic
;        print, k, t0(k), xin(i)
        if (t0(k) GT xin(i)) then begin
          khi(i) = k
        endif else begin
          klo(i) = k
        endelse
  endwhile
endfor  
 

 
  xlo = t0(klo)
  xhi = t0(khi)
  ylo = z0(klo)
  yhi = z0(khi)
  y2lo = z2(klo)
  y2hi = z2(khi)
 
  h = xhi - xlo

;print, 'splint'
;print, t0
;print, z0

;print, klo
;print, khi
;print, xlo
;print, xhi
;print, ylo
;print, yhi
;print, yhi-ylo
;print, y2lo
;print, y2hi

; Compute the answer

a = (xhi - xin)/h
b = (xin - xlo)/h

;print, a
;print, b

  case deriv of
    0: yinterp = b*yhi + a*ylo + ((b^2-1.0)*b*y2hi   + (a^2-1.0)*a*y2lo)*h^2/6.0 
    1: yinterp = (yhi-ylo)/h   + ((3.0*b^2-1.0)*y2hi - (3.0*a^2-1.0)*y2lo)*h/6.0
    2: yinterp = b*y2hi + a*y2lo 
  endcase


return, yinterp

end
