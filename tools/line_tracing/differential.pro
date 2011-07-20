function Differential, y, V
;
;-right hand side of system of two ODEs describing evolution
;of magnetic field line in [ix,iz] plane;
;
; Inputs:
; iy is the poloidal index
; V is a 2-vector [ix,iz]
;
; Output: 2-vector of the right-hand-side
;=====================================================;


;-coords of the location where evaluation is desired
   x=V[0]
   z=V[1]

   C = Fun_C12(x,y,z)
;   C1=Fun_C1(x,y,z)        ;Eq1: d(ix)/d(iy) = C_1
;   C2=Fun_C2(x,y,z)        ;Eq2: d(iz)/d(iy) = C_2

   C[1] = 0.0 ; For now no Z

;
;
;
return, C
end

