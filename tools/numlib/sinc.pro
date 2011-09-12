function sinc,x
;
; Author:       Ilon Joseph
; Began:        2011/06/22
; Modified:     2011/06/22 
;
; mathematical function sinc=sin(x)/x
; returns an answer with the same precision as x 

  ans = sin(x)
  ans(where(x ne 0.0)) /= x(where(x ne 0.0))
  ans(where(x eq 0.0))  = 1.0

  return, ans
end
