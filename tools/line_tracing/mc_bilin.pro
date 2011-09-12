; Bilinear interpolation for equilibrium field values on grid
; Input is x,y coordinates
; Output structure = {rxy, bxy, bpxy, btxy, sinty, hthe, bxcvx, bxcvz}
FUNCTION mc_bilin, x, y
  COMMON griddata, g, deltaZtor, Ntor
  COMMON metric_data, mc_data

  ;; Find lower left cell indices
  xind=INTERPOL(FINDGEN(g.nx),g.psixy[*,0],x)   ;float in [0,g.nx)
  ix_int=FIX(xind)
  IF ( (ix_int LT 1) OR (ix_int GE (g.nx-1)) ) THEN STOP  ;exits domain...
  yind=DOUBLE(g.ny)*y/(2.*!DPI)        ;float in [0,g.ny)
  iy_int=FIX(yind)
  IF iy_int EQ g.ny THEN STOP   ;shouldn't get to y=2*PI but stop if we do

  ;; Calculate normalized cell coordinates
  dx=g.psixy[ix_int+1,0]-g.psixy[ix_int,0]
  x_0=g.psixy[ix_int,0]
  xd=(x-x_0)/dx                       ;normalized x 
  dy=(2.*!DPI)/DOUBLE(g.ny)
  y_0=iy_int*dy
  yd=(y-y_0)/dy                       ;normalized y

  ;form x_y_vec [ (1-xd)*(1-yd), xd*(1-yd), (1-xd)*yd, xd*yd ]^T
  x_y_vec= [ [(1.0-xd)*(1.0-yd)], [xd*(1.0-yd)], [(1.0-xd)*yd], [xd*yd] ]

  IF mc_data[ix_int,iy_int].TF EQ 0 THEN BEGIN
    mc_data[ix_int,iy_int].TF=1

    ;do bound-checking for y indices,
    ;metric coeffs are periodic across y=2*PI
    iy_intP1=iy_int+1
    IF iy_intP1 EQ g.ny THEN iy_intP1=0

    mc_data[ix_int,iy_int].rxy = [ g.rxy[ix_int  ,iy_int  ] ,                  $
                                   g.rxy[ix_int+1,iy_int  ] ,                  $
                                   g.rxy[ix_int  ,iy_intP1] ,                  $
                                   g.rxy[ix_int+1,iy_intP1] ]

    mc_data[ix_int,iy_int].bxy = [ g.bxy[ix_int  ,iy_int  ] ,                  $
                                   g.bxy[ix_int+1,iy_int  ] ,                  $
                                   g.bxy[ix_int  ,iy_intP1] ,                  $
                                   g.bxy[ix_int+1,iy_intP1] ]

    mc_data[ix_int,iy_int].bpxy = [ g.bpxy[ix_int  ,iy_int  ] ,                $
                                    g.bpxy[ix_int+1,iy_int  ] ,                $
                                    g.bpxy[ix_int  ,iy_intP1] ,                $
                                    g.bpxy[ix_int+1,iy_intP1] ]

    mc_data[ix_int,iy_int].btxy = [ g.btxy[ix_int  ,iy_int  ] ,                $
                                    g.btxy[ix_int+1,iy_int  ] ,                $
                                    g.btxy[ix_int  ,iy_intP1] ,                $
                                    g.btxy[ix_int+1,iy_intP1] ]

    mc_data[ix_int,iy_int].sinty = [ g.sinty[ix_int  ,iy_int  ] ,              $
                                     g.sinty[ix_int+1,iy_int  ] ,              $
                                     g.sinty[ix_int  ,iy_intP1] ,              $
                                     g.sinty[ix_int+1,iy_intP1] ]

    mc_data[ix_int,iy_int].hthe = [ g.hthe[ix_int  ,iy_int  ] ,                $
                                    g.hthe[ix_int+1,iy_int  ] ,                $
                                    g.hthe[ix_int  ,iy_intP1] ,                $
                                    g.hthe[ix_int+1,iy_intP1] ]

    mc_data[ix_int,iy_int].bxcvx = [ g.bxcvx[ix_int  ,iy_int  ] ,              $
                                     g.bxcvx[ix_int+1,iy_int  ] ,              $
                                     g.bxcvx[ix_int  ,iy_intP1] ,              $
                                     g.bxcvx[ix_int+1,iy_intP1] ]

    mc_data[ix_int,iy_int].bxcvz = [ g.bxcvz[ix_int  ,iy_int  ] ,              $
                                     g.bxcvz[ix_int+1,iy_int  ] ,              $
                                     g.bxcvz[ix_int  ,iy_intP1] ,              $
                                     g.bxcvz[ix_int+1,iy_intP1] ]

  ENDIF

  rxy   = x_y_vec#mc_data[ix_int,iy_int].rxy
  bxy   = x_y_vec#mc_data[ix_int,iy_int].bxy
  bpxy  = x_y_vec#mc_data[ix_int,iy_int].bpxy
  btxy  = x_y_vec#mc_data[ix_int,iy_int].btxy
  sinty = x_y_vec#mc_data[ix_int,iy_int].sinty
  hthe  = x_y_vec#mc_data[ix_int,iy_int].hthe
  bxcvx = x_y_vec#mc_data[ix_int,iy_int].bxcvx
  bxcvz = x_y_vec#mc_data[ix_int,iy_int].bxcvz

  RETURN, { rxy:rxy, bxy:bxy, bpxy:bpxy, btxy:btxy, sinty:sinty, hthe:hthe,    $
            bxcvx:bxcvx, bxcvz:bxcvz                                       }
END
