;; Calculates dispersion relation using FZ_ROOTS
;; Adapted from complexRoots

;; Result depends on mu = (c * kperp / wpe)^2  and sparsperp = wci*wce*(kpar*c/wpe)^2

FUNCTION calc_cubic, x, c
  RETURN, c[0] + c[1]*x + c[2]*x*x + c[3]*x*x*x
END

FUNCTION cubic_root, coef
  a = coef[2] / coef[3]
  b = coef[1] / coef[3]
  c = coef[0] / coef[3]

  p = b - a*a/3.
  q = c + (2.*a*a*a - 9.*a*b) / 27.

  u1 = (-q/2. + SQRT(q*q/4. + p*p*p/27.))^(1./3.)

  u2 = COMPLEX(-0.5, SQRT(3.)/2.) * u1
  u3 = COMPLEX(-0.5, -SQRT(3.)/2.) * u1

  x1 = -p/(3.*u1) + u1 - a/3.
  x2 = -p/(3.*u2) + u2 - a/3.
  x3 = -p/(3.*u3) + u3 - a/3.

  return, [x1,x2,x3]
END

FUNCTION dispersion, mu, sparsperp, linear=linear, epar=epar, gradp=gradp

 if keyword_set(LINEAR) then begin
    s=1e-3 + 1e2*findgen(1001)/1000. 
 endif else begin
    s=1e-3*10^(6*findgen(100)/101)
 endelse

 N = N_ELEMENTS(s)
 
 cubic = 0
 IF KEYWORD_SET(epar) OR KEYWORD_SET(gradp) THEN BEGIN
   C = COMPLEXARR(4)
   cubic = 1
 ENDIF ELSE C = COMPLEXARR(3)

 w1 = COMPLEXARR(N)
 w2 = w1
 w3 = w1

 FOR i=0, N-1 DO BEGIN
   C[0] = COMPLEX(0.0, -1.*s[i])
   C[1] = COMPLEX(0.0, s[i])

   C[2] = COMPLEX(1.0)
   IF KEYWORD_SET(gradp) THEN C[2] = C[2] + COMPLEX(0.0, s[i]/sparsperp) 
   
   IF cubic THEN BEGIN
     ; problem becomes cubic
     C[3] = COMPLEX(0.0, -1.*mu*s[i] / sparsperp)
     IF KEYWORD_SET(epar) THEN C[3] = C[3] - COMPLEX(0.0, s[i]/sparsperp)

     result = FZ_ROOTS(C, EPS=1e-6)

     IF imaginary(result[2]) GT imaginary(result[1]) THEN BEGIN
       tmp = result[2]
       result[2] = result[1]
       result[1] = tmp
     ENDIF

     ;result = cubic_root(c)
     
     ;PRINT, s[i], c[0], c[1], c[2], c[3]

     ;PRINT, "-> ", calc_cubic(result[0], c), calc_cubic(result[1], c), calc_cubic(result[2], c)
   ENDIF ELSE result = FZ_ROOTS(C) 

   w1[i] = result[0]
   w2[i] = result[1]
   IF KEYWORD_SET(epar) OR KEYWORD_SET(gradp) THEN w3[i] = result[2]
 ENDFOR
 
 RETURN, {mu:mu, sparsperp:sparsperp, s:s, w1:w1, w2:w2, w3:w3}
END
