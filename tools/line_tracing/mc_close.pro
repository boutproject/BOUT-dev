; Close neighbor interpolation for equilibrium field values on grid
; This is NOT true nearest neighbor, but is the original scheme used
; Input is x,y coordinates
; Output structure = {rxy, bxy, bpxy, btxy, sinty, hthe, bxcvx, bxcvz}
FUNCTION mc_close, x, y
  COMMON griddata, g, deltaZtor, Ntor
  COMMON metric_data, mc_data

  ;; Find lower left cell indices
  xind=INTERPOL(FINDGEN(g.nx),g.psixy[*,0],x)   ;float in [0,g.nx)
  ix_int=FIX(xind)
  IF ( (ix_int LT 1) OR (ix_int GE (g.nx-1)) ) THEN STOP  ;exits domain...
  yind=DOUBLE(g.ny)*y/(2.*!DPI)        ;float in [0,g.ny)
  iy_int=FIX(yind)
  IF iy_int EQ g.ny THEN STOP   ;shouldn't get to y=2*PI but stop if we do

  rxy   = g.rxy[ix_int,iy_int]
  bxy   = g.bxy[ix_int,iy_int]
  bpxy  = g.bpxy[ix_int,iy_int]
  btxy  = g.btxy[ix_int,iy_int]
  sinty = g.sinty[ix_int,iy_int]
  hthe  = g.hthe[ix_int,iy_int]
  bxcvx = g.bxcvx[ix_int,iy_int]
  bxcvz = g.bxcvz[ix_int,iy_int]

  RETURN, { rxy:rxy, bxy:bxy, bpxy:bpxy, btxy:btxy, sinty:sinty, hthe:hthe,    $
            bxcvx:bxcvx, bxcvz:bxcvz                                       }
END
