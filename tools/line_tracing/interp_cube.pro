FUNCTION interp_cube, f

;Given function values f(x=-1),f(x=0),f(x=1),f(x=2) calculate
;coefficients of P(x)=a_0+a_1*x+a_2*x^2+a_3*x^3
;estimates derivative based on nearest two neighbors
  f0_prime=(f[2]-f[0])/2.
  f1_prime=(f[3]-f[1])/2.
  A_inv=[[1,0,0,0],[0,0,1,0],[-3,3,-2,-1],[2,-2,1,1]]
  coeffs=transpose(A_inv)#[f[1],f[2],f0_prime,f1_prime]
  RETURN,coeffs
END
