FUNCTION FMODULO, x,y
; Generalization of modulo to negative x.
;
; Example: 
;       fmodulo(-2*!PI+0.1, 2*!PI) -> 0.1
; while 
;       mod(-2*!PI+0.1, 2*!PI) -> -6.18319
;===================================================;
  IF (x GE 0.) THEN BEGIN 
    res=(x MOD y)
  ENDIF ELSE BEGIN
    res=( (x+(CEIL(ABS(x/y)))*y) MOD y )
  ENDELSE

 RETURN, res
END
