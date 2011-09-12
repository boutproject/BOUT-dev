function splgen, x, y, y1a=y1a, y1b=y1b, period=period, tridiag=tridiag, cycdiag=cycdiag, lapack=lapack
;
; Author:       Ilon Joseph
; Began:        2011/06/03
; Modified:     2011/06/07 
;
; generates cubic spline 2nd derivatives based on Numerical Recipes in C
; equivalent to IDL function spl_init with enhanced capability for periodic bc's
; y2 = second derivatives = output of IDL function splgen
;
; Options
;       * default  = natural spline with internal tridiagonal solver 
;       * y1a      = set dydx on left boundary
;       * y1b      = set dydx on right boundary       
;       * tridiag  = natural spline with tridiag.pro solver
;       * cycdiag  = periodic spline with cycdiag.pro solver
;       * period = period length for periodic spline
;         assumes that x0 has branch cut at boundary which jumps by 1 period
;       * lapack   = use lapack versions of tridiag and cycdiag (la_trisol)        


n  = size(x,/n_elements)
ny = size(y,/n_elements)

if n  NE ny then begin
        print, "splgen:  length(x) != length(y)"
        return, 0 
endif 

  dx    = x[1:n-1]-x[0:n-2]             ; n-1
  dy    = y[1:n-1]-y[0:n-2]             ; n-1
  dydx  = dy/dx                         ; n-1
  dx2   = x[2:n-1]-x[0:n-3]             ; n-2 
  d2ydx = dydx[1:n-2] - dydx[0:n-3]     ; n-2 "initial u"

;  print, "splgen"
;  print, dx2
;  print, d2ydx

if keyword_set(tridiag) then begin
  r = 6.0*d2ydx
  b = 2.0*dx2   	; n-2
  a = [0.0, dx[1:n-3]]  ; n-2
  c = [dx[1:n-3], 0.0]  ; n-2

  if not keyword_set(lapack) then begin
     y2 = tridiag(a,b,c,r)
  endif else begin
     b1=b        ; n-2
     a1=a[1:n-3] ; n-3
     c1=c[0:n-4] ; n-3
     la_tridc, a1, b1, c1, b2, index
     y2 = la_trisol(a1, b1, c1, b2, index, r) 
  endelse
   
  return, [0.0,y2,0.0]
endif

if keyword_set(period) or keyword_set(cycdiag) then begin
  if not keyword_set(period) then begin
        print, "splgen: period must be specified"
        return, 0 
  endif

  dx    = [dx,   x[0]-x[n-1]+period]
  dy    = [dy,   y[0]-y[n-1]]
  dydx  = dy/dx
  dx2   = [x[1]-x[n-1]+period, dx2,   x[0]-x[n-2]+period]
  d2ydx = [dydx[0]-dydx[n-1],  d2ydx, dydx[n-1]-dydx[n-2]]

 
;  print, dx2
;  print, d2ydx

  r = 6.0*d2ydx
  b = 2.0*dx2
  a = [0.0, dx[1:n-1]]
  c = [dx[1:n-1], 0.0]
  alpha =  dx[0] ; = c[n-1]
  beta  =  dx[0] ; = a[0]

  if not keyword_set(lapack) then begin
    y2 = cycdiag(a,b,c,r,alpha,beta)
  endif else begin  
    y2 = cycdiag(a,b,c,r,alpha,beta,/lapack)
  endelse

  return, y2
endif

;Internal tridiagonal solver

  y2  = 0.0*y                           ; n
  sx  = dx(0:n-3)/dx2(0:n-3)            ; n-2
  u   = [0.0,d2ydx,0.0]			; "u"

  if keyword_set(y1a) then begin
        y2[0] =-0.5
        u[0]  = 3.0*(dydx[0] - y1a)/dx[0]
  endif

  for j=0,n-3 do begin
    s       = dx[j]/dx2[j];
    p       = 2.0 + s*y2[j];
    u[j+1]  = (6.0*u[j+1]/dx2[j] - s*u[j])/p;
    y2[j+1] = (s-1.0)/p;
  endfor

  if keyword_set(y1b) then begin
        u[n-1] =3.0*(y1b-dydx[n-2])/dx[n-2]
        y2[n-1]=(2.0*u[n-1]-u[n-2])/(2.0+y2[n-2])
  endif

  for k=n-2,0,-1 do begin
        y2[k] = y2[k+1]*y2[k] + u[k]
  endfor

  return, y2

end

