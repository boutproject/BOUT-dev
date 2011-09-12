FUNCTION Differential, y, V
;-right hand side of system of two ODEs describing evolution
;of magnetic field line in [ix,iz] plane;
;
; Inputs:
; y is the poloidal location in [0,2*PI)
; V is a 2-vector [x,z]
;
; Output: 2-vector of the right-hand-side
;=====================================================;

;-coords of the location where evaluation is desired
  x=V[0]
  z=V[1]

  C = Fun_C12(x,y,z)

; (JPS) This is what was originally in the code, but testing with the updated Z
; changes the plots significantly, just need to keep in [0,deltaZtor)
;  C[1] = 0.0 ; For now no Z

  RETURN, C
END
