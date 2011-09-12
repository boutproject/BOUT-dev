function triprod, a, b, c, x
;
; Author:       Ilon Joseph
; Began:        2011/06/03
; Modified:     2011/06/03 
; 
; Product of M*x=r for x where the matrix M is of the tridiagonal form 
; M = [ A, B, C, 0 ...]
;
; Note that A(0) and C(n-1) are not used in the algorithm
;

  n  = size(b,/n_elements)
  na = size(a,/n_elements) 
  nc = size(c,/n_elements) 
  nx = size(x,/n_elements) 

  if n NE nx  then begin
    print, 'tridiag: length(x)  != length(B)'
    return, 0
  endif   
  if ( (n NE na) AND (n NE na+1) ) then begin
    print, 'tridiag: length(A)  != length(B)'
    return, 0
  endif 
  if na EQ n-1 then a=[0.0,a]
  if ( (n NE nc) AND (n NE nc+1) )  then begin
    print, 'tridiag: length(C)  != length(B)'
    return, 0
  endif   
  if nc EQ n-1 then c=[c,0.0]
  
  ans = [0,a(1:n-1)*x(0:n-2)] + b*x + [c(0:n-2)*x(1:n-1),0]

  return, ans
end


function tridiag, a, b, c, r, gamma=gamma, beta=beta, tol=tol, lapack=lapack
;
; Author:       Ilon Joseph
; Began:        2011/06/03
; Modified:     2011/06/07 
; 
; Solves M*x=r for x where the matrix M is of the tridiagonal form
;
; M = [B0, C0, 0,        ...] [X0]  = [R0]
;     [A1, B1, C1, 0,    ...] [X1]  = [R1]
;      ...
;     [..., 0, A(n-1), B(n-1)][X(n-1)] = [R(n-1)]
;
; Note that A(0) and C(n-1) are not used in the algorithm
; No pivoting --> can lead to failure
;
; Based on Numerical Recipes in C
; Equivalent to IDL routine trisol
;
; Options:
;   *tol = reports error if |M*x-r|/|r| < tol
;   *beta, gamma = return beta, gamma working vectors
;   *lapack = use lapack routines la_tridc, la_trisol which include pivoting instead

  n  = size(b,/n_elements)
  na = size(a,/n_elements) 
  nc = size(c,/n_elements) 
  nr = size(r,/n_elements) 

  if n NE nr  then begin
    print, 'tridiag: length(r)  != length(B)'
    return, 0
  endif   
  if ( (n NE na) AND (n NE na+1) ) then begin
    print, 'tridiag: length(A)  != length(B)'
    return, 0
  endif 
  if na EQ n-1 then a=[0.0,a]
  if ( (n NE nc) AND (n NE nc+1) )  then begin
    print, 'tridiag: length(C)  != length(B)'
    return, 0
  endif   
  if nc EQ n-1 then c=[c,0.0]

if not keyword_set(lapack) then begin

  x=0.0*r
  ZERO = 0.0*b[0]

  if b[0] EQ ZERO then begin
    print, 'tridiag: matrix degenerate B(0)=0, rewrite system with N-1 eqs'
    return, 0
  endif

; LU decomp
  gamma = b*0.0
  beta  = b
  x[0] = r[0]/b[0] 

  for j=1,n-1 do begin
    gamma[j]  = c[j-1]/beta[j-1]
    beta[j]  -= a[j]*gamma[j]

    if beta[j] EQ ZERO then begin
        print, 'tridiag failed: either matrix is degenerate or pivoting is required'
        return, 0
    endif

    x[j] = (r[j]-a[j]*x[j-1])/beta[j];
  endfor  

; Backsub
  for j=n-1,1,-1 do begin
    x[j-1] -= x[j]*gamma[j]
  endfor

  if keyword_set(tol) then begin
    r1 = triprod(a, b, c, x)
    error = max(abs(r1-r))/max(abs(r)) 
    if error GT TOL then begin
        print, 'tridiag: large relative error =',error
    endif
  endif

endif else begin
  a1=a(1:n-1)
  b1=b
  c1=c(0:n-2)

  la_tridc,a1,b1,c1,b2,index
  x=la_trisol(a1,b1,c1,b2,index,r)
  
endelse

  return, x
end



