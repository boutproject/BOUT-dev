function cycprod, a, b, c, x, alpha, beta
;
; Author:       Ilon Joseph
; Began:        2011/06/03
; Modified:     2011/06/03 
; 
; Product of M*x where the matrix M is of the cyclic tridiagonal form 
; M = [ B, C, 0 ...beta]
;     [ alpha, ... A, B]

  n  = size(b,/n_elements)
  na = size(a,/n_elements) 
  nc = size(c,/n_elements) 
  nx = size(x,/n_elements) 

  if n NE nx  then begin
    print, 'cycdiag: length(x)  != length(B)'
    return, 0
  endif   
  if ( (n NE na) AND (n NE na+1) ) then begin
    print, 'cycdiag: length(A)  != length(B)'
    return, 0
  endif 
  if na EQ n-1 then a=[0.0,a]
  if ( (n NE nc) AND (n NE nc+1) )  then begin
    print, 'cycdiag: length(C)  != length(B)'
    return, 0
  endif   
  if nc EQ n-1 then c=[c,0.0]

  ans = [beta*x[n-1],a[1:n-1]*x[0:n-2]] + b*x + [c[0:n-2]*x[1:n-1],alpha*x[0]]

  return, ans
end


function cycdiag, a, b, c, r, alpha, beta, tol=tol, lapack=lapack
;
; Author:       Ilon Joseph
; Began:        2011/06/03
; Modified:     2011/06/03 
; 
; Solves M*x=r for x where the matrix M is of the cyclic tridiagonal form
;
; M = [B0, C0, 0,               ...beta] [X(0)]   = [R(0)]
;     [A1, B1, C1, 0,           ...   0] [X(1)]   = [R(1)]
;     [0,...     A(n-2), B(n-2), C(n-2)] [X(n-2)] = [R(n-2)]
;     [alpha, 0, ...  0, A(n-1), B(n-1)] [X(n-1)] = [R(n-1)]
;
; Note that A[0[ and C[n-1] are not used in the algorithm
;
; Based on Numerical Recipes in C
;
; Options
;   *tol = reports error if |M*x-r|/|r| < tol
;   *lapack = use lapack-based la_trisol which includes pivoting

  n  = size(b,/n_elements)
  na = size(a,/n_elements) 
  nc = size(c,/n_elements) 
  nr = size(r,/n_elements) 

  if n LT 3  then begin
    print, 'cycdiag: length too small (< 3)'
    return, 0
  endif  
  if n NE nr  then begin
    print, 'tridiag: length(r)  != length(B)'
    return, 0
  endif   
  if ( (n NE na) AND (n NE na+1) ) then begin
    print, 'cycdiag: length(A)  != length(B)'
    return, 0
  endif 
  if na EQ n-1 then a=[0.0,a]
  if ( (n NE nc) AND (n NE nc+1) )  then begin
    print, 'cycdiag: length(C)  != length(B)'
    return, 0
  endif   
  if nc EQ n-1 then c=[c,0.0]

  ZERO = 0.0*b[0]

  if b[0] EQ ZERO then begin
    print, 'cycdiag: matrix degenerate B(0)=0'
    return, 0
  endif

  ; The NR choice gamma=-b[0] throws a tridiag error when b[1]=2*b[0]
  ; instead use an irrational factor close to -1
  gfac = (sqrt(5.0)+1)/3.0
  gamma   = -b[0]*gfac

  bb      = b
  bb[0]  -= gamma
  bb[n-1]-= alpha*beta/gamma
  u      = 0.0*bb
  u[0]   = gamma
  u[n-1] = alpha

  if not keyword_set(lapack) then begin
    x=tridiag(a,bb,c,r)
    z=tridiag(a,bb,c,u) ; this step can cause underflow
    
  endif else begin
    la_a = a[1:n-1]
    la_b = bb
    la_c = c[0:n-2]
    la_tridc,la_a,la_b,la_c,la2,index
    x=la_trisol(la_a,la_b,la_c,la2,index,r)
    z=la_trisol(la_a,la_b,la_c,la2,index,u) ; this step can cause underflow
  endelse
 
  fac = (x[0]+x[n-1]*beta/gamma)/(1.0+z[0]+z[n-1]*beta/gamma)
  x  -= z*fac

  if keyword_set(tol) then begin
    r1 = cycprod(a, b, c, x, alpha, beta)
    error = max(abs(r1-r))/max(abs(r)) 
    if error GT TOL then begin
        print, 'cycdiag: large relative error =',error
    endif
  endif

  return, x
end
