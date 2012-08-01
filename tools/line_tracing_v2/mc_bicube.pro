; Bicubic interpolation for equilibrium field values on grid
; Input is x,y coordinates
; Output structure = {rxy, bxy, bpxy, btxy, sinty, hthe, bxcvx, bxcvy, bxcvz}
FUNCTION mc_bicube, x, y
  COMMON griddata, g, deltaZtor, Ntor
  COMMON metric_data, mc_data

  ;; Find lower left cell indices
  xind=INTERPOL(FINDGEN(g.nx),g.psixy[*,0],x)   ;float in [0,g.nx)
  ix_int=FIX(xind)
  IF ( (ix_int LT 1) OR (ix_int GE (g.nx-2)) ) THEN STOP,'Exit radial domain...'
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

  ;form x_y_vec [x^i*y^j]
  x_y_vec= [ [xd^0*yd^0],[xd^1*yd^0],[xd^2*yd^0],[xd^3*yd^0] ,                 $
             [xd^0*yd^1],[xd^1*yd^1],[xd^2*yd^1],[xd^3*yd^1] ,                 $
             [xd^0*yd^2],[xd^1*yd^2],[xd^2*yd^2],[xd^3*yd^2] ,                 $
             [xd^0*yd^3],[xd^1*yd^3],[xd^2*yd^3],[xd^3*yd^3] ]

  IF mc_data[ix_int,iy_int].TF EQ 0 THEN BEGIN
    mc_data[ix_int,iy_int].TF=1
    ;compute the normalized grid pts
    x0d=0.0                             ;normalized x_0
    x1d=1.0                             ;normalized x_1
    xmd=(g.psixy[ix_int-1,0]-x_0)/dx    ;normalized x_-1
    x2d=(g.psixy[ix_int+2,0]-x_0)/dx    ;normalized x_2
    xd_vec=[xmd,x0d,x1d,x2d]
    y0d=0.0                             ;normalized y_0
    y1d=1.0                             ;normalized y_1
    ymd=-1.0                            ;normalized y_-1
    y2d=2.0                             ;normalized y_2
    yd_vec=[ymd,y0d,y1d,y2d]

    ;do bound-checking for y indices,
    ;metric coeffs are periodic across y=2*PI
    iy_intM1=iy_int-1
    iy_intP1=iy_int+1
    iy_intP2=iy_int+2
    CASE iy_int OF
      0      : iy_intM1=g.ny-1
      g.ny-2 : iy_intP2=0
      g.ny-1 : BEGIN
                 iy_intP1=0
                 iy_intP2=1
               END
      else   :
    ENDCASE

    vals=[ [g.rxy[(ix_int-1):(ix_int+2),iy_intM1]],                            $
           [g.rxy[(ix_int-1):(ix_int+2),iy_int  ]],                            $
           [g.rxy[(ix_int-1):(ix_int+2),iy_intP1]],                            $
           [g.rxy[(ix_int-1):(ix_int+2),iy_intP2]] ]
    mc_data[ix_int,iy_int].rxy=interp_bicube(vals,xd_vec,yd_vec)

    vals=[ [g.bxy[(ix_int-1):(ix_int+2),iy_intM1]],                            $
           [g.bxy[(ix_int-1):(ix_int+2),iy_int  ]],                            $
           [g.bxy[(ix_int-1):(ix_int+2),iy_intP1]],                            $
           [g.bxy[(ix_int-1):(ix_int+2),iy_intP2]] ]
    mc_data[ix_int,iy_int].bxy=interp_bicube(vals,xd_vec,yd_vec)

    vals=[ [g.bpxy[(ix_int-1):(ix_int+2),iy_intM1]],                           $
           [g.bpxy[(ix_int-1):(ix_int+2),iy_int  ]],                           $
           [g.bpxy[(ix_int-1):(ix_int+2),iy_intP1]],                           $
           [g.bpxy[(ix_int-1):(ix_int+2),iy_intP2]] ]
    mc_data[ix_int,iy_int].bpxy=interp_bicube(vals,xd_vec,yd_vec)

    vals=[ [g.btxy[(ix_int-1):(ix_int+2),iy_intM1]],                           $
           [g.btxy[(ix_int-1):(ix_int+2),iy_int  ]],                           $
           [g.btxy[(ix_int-1):(ix_int+2),iy_intP1]],                           $
           [g.btxy[(ix_int-1):(ix_int+2),iy_intP2]] ]
    mc_data[ix_int,iy_int].btxy=interp_bicube(vals,xd_vec,yd_vec)

    vals=[ [g.sinty[(ix_int-1):(ix_int+2),iy_intM1]],                          $
           [g.sinty[(ix_int-1):(ix_int+2),iy_int  ]],                          $
           [g.sinty[(ix_int-1):(ix_int+2),iy_intP1]],                          $
           [g.sinty[(ix_int-1):(ix_int+2),iy_intP2]] ]
    mc_data[ix_int,iy_int].sinty=interp_bicube(vals,xd_vec,yd_vec)

    vals=[ [g.hthe[(ix_int-1):(ix_int+2),iy_intM1]],                           $
           [g.hthe[(ix_int-1):(ix_int+2),iy_int  ]],                           $
           [g.hthe[(ix_int-1):(ix_int+2),iy_intP1]],                           $
           [g.hthe[(ix_int-1):(ix_int+2),iy_intP2]] ]
    mc_data[ix_int,iy_int].hthe=interp_bicube(vals,xd_vec,yd_vec)

    vals=[ [g.bxcvx[(ix_int-1):(ix_int+2),iy_intM1]],                          $
           [g.bxcvx[(ix_int-1):(ix_int+2),iy_int  ]],                          $
           [g.bxcvx[(ix_int-1):(ix_int+2),iy_intP1]],                          $
           [g.bxcvx[(ix_int-1):(ix_int+2),iy_intP2]] ]
    mc_data[ix_int,iy_int].bxcvx=interp_bicube(vals,xd_vec,yd_vec)

    vals=[ [g.bxcvy[(ix_int-1):(ix_int+2),iy_intM1]],                          $
           [g.bxcvy[(ix_int-1):(ix_int+2),iy_int  ]],                          $
           [g.bxcvy[(ix_int-1):(ix_int+2),iy_intP1]],                          $
           [g.bxcvy[(ix_int-1):(ix_int+2),iy_intP2]] ]
    mc_data[ix_int,iy_int].bxcvy=interp_bicube(vals,xd_vec,yd_vec)

    vals=[ [g.bxcvz[(ix_int-1):(ix_int+2),iy_intM1]],                          $
           [g.bxcvz[(ix_int-1):(ix_int+2),iy_int  ]],                          $
           [g.bxcvz[(ix_int-1):(ix_int+2),iy_intP1]],                          $
           [g.bxcvz[(ix_int-1):(ix_int+2),iy_intP2]] ]
    mc_data[ix_int,iy_int].bxcvz=interp_bicube(vals,xd_vec,yd_vec)

    vals=[ [g.jpar0[(ix_int-1):(ix_int+2),iy_intM1]],                          $
           [g.jpar0[(ix_int-1):(ix_int+2),iy_int  ]],                          $
           [g.jpar0[(ix_int-1):(ix_int+2),iy_intP1]],                          $
           [g.jpar0[(ix_int-1):(ix_int+2),iy_intP2]] ]
    mc_data[ix_int,iy_int].jpar0=interp_bicube(vals,xd_vec,yd_vec)

  ENDIF

  rxy   = x_y_vec#mc_data[ix_int,iy_int].rxy
  bxy   = x_y_vec#mc_data[ix_int,iy_int].bxy
  bpxy  = x_y_vec#mc_data[ix_int,iy_int].bpxy
  btxy  = x_y_vec#mc_data[ix_int,iy_int].btxy
  sinty = x_y_vec#mc_data[ix_int,iy_int].sinty
  hthe  = x_y_vec#mc_data[ix_int,iy_int].hthe
  bxcvx = x_y_vec#mc_data[ix_int,iy_int].bxcvx
  bxcvy = x_y_vec#mc_data[ix_int,iy_int].bxcvy
  bxcvz = x_y_vec#mc_data[ix_int,iy_int].bxcvz
  jpar0 = x_y_vec#mc_data[ix_int,iy_int].jpar0

  RETURN, { rxy:rxy, bxy:bxy, bpxy:bpxy, btxy:btxy, sinty:sinty, hthe:hthe,    $
            bxcvx:bxcvx, bxcvy:bxcvy, bxcvz:bxcvz, jpar0:jpar0             }
END
