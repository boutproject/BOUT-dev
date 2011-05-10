FUNCTION theta_differential, s, pos
;
; Expression of differentials
; dr/ds = Br/B
; dz/ds = Bz/B
;
; Inputs: x->s (distance along theta line)
; y[0:1] -> [r,z] (coords)
;------------------------------------------------;
  COMMON td_com, dctf, lastgoodpos
  a = local_gradient(dctf, pos[0], pos[1], status=status)
  
  IF status EQ 0 THEN BEGIN
    ; No error in localgradient.
    lastgoodpos = pos
    
    ; If status NE 0 then an error occurred.
    ; Allow a.dfdz to cause an error so escape LSODE
  ENDIF
  
  Br = a.dfdz
  Bz = -a.dfdr
  B = SQRT(Br^2 + Bz^2)

  RETURN, [Br/B, Bz/B]
END
