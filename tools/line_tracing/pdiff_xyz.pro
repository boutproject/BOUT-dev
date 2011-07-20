function pdiff_xyz, x, y, z, f, debug=debug
;
; Calculate partial derivatives df/dx, df/dy, and df/dz
; for function f given on a set of n points (x,y,z)
; using regression
;
; Inputs: vectors x[n],y[n],z[n],f[n] - same size
; Output: structure {x:df/dx, y:df/dy, z:df/dz}
;-------------------------------------------------

; Combine vectors of independent variable into a 3x8 array  
v = [TRANSPOSE(x), TRANSPOSE(y), TRANSPOSE(z)]  
  

; Compute the fit
res = REGRESS(v, f, CONST=const)  


  
if keyword_set(DEBUG) then begin
    print, "f=", const, " + ", res[0],"*x + ", res[1],"*y + ", res[2],"*z"
endif

;
;
;
return, {x:res[0],y:res[1],z:res[2]}
end

