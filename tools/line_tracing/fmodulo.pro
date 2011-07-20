function FMOD, x,y
;
; Cast x within the period [0,y]
;===================================================;

 return, (x mod y)
end


function FMODULO, x,y
;
; Generalization of FMOD to negative x.
;
; Example: 
;       fmodulo(-2*!PI+0.1, 2*!PI) -> 0.1
; while 
;       fmod(-2*!PI+0.1, 2*!PI) -> -6.18319
;===================================================;

  if (x ge 0.) then begin 
    res=FMOD(x,y)
  endif else begin
    res=FMOD(x+(CEIL(ABS(x/y)))*y,y)
  endelse

 return, res
end
