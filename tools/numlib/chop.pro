function chop, x, tol=tol
;
; Author:       Ilon Joseph
; Began:        2011/06/22
; Modified:     2011/06/22 
;
; "chop" filter 
; sets elements of x to 0 if abs(x)<tol*max(abs(x))
;       returns answer of type(x)
;       if tol is unspecified, chop sets tol based on type(x)
;               tol(float)  = 1e-7
;               tol(double) = 1e-15
;       these values were set based on accuracy of fft

 ans = x
 xmax=max(abs(x))

 if not keyword_set(tol) then begin
   st = size(x,/type)
   if st EQ size(float(1),/type) OR st EQ size(complex(1,0),/type) then tol=1e-7
   if st EQ size(double(1),/type) OR st EQ size(dcomplex(1,0),/type)  then tol=1e-15
 endif
 xmin = xmax*abs(tol)
; print, tol, xmin
; print,  where(abs(x) LT xmin)

 mask = where(abs(x) LT xmin)
 if mask(0) ne -1 then begin
      ans(mask)  = 0.0
 endif

 return, ans
end
