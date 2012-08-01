; Return the inverse of the 64x64 matrix A that satisfies Ax=b
;   x= coefficient vector (a_ijk)
;   b= function values vector 
FUNCTION make_tricube_mat, debug=debug

  A_01=make_tricube_vec(0,0,0,derivative='f')
  A_02=make_tricube_vec(1,0,0,derivative='f')
  A_03=make_tricube_vec(0,1,0,derivative='f')
  A_04=make_tricube_vec(1,1,0,derivative='f')
  A_05=make_tricube_vec(0,0,1,derivative='f')
  A_06=make_tricube_vec(1,0,1,derivative='f')
  A_07=make_tricube_vec(0,1,1,derivative='f')
  A_08=make_tricube_vec(1,1,1,derivative='f')

  A_09=make_tricube_vec(0,0,0,derivative='dx')
  A_10=make_tricube_vec(1,0,0,derivative='dx')
  A_11=make_tricube_vec(0,1,0,derivative='dx')
  A_12=make_tricube_vec(1,1,0,derivative='dx')
  A_13=make_tricube_vec(0,0,1,derivative='dx')
  A_14=make_tricube_vec(1,0,1,derivative='dx')
  A_15=make_tricube_vec(0,1,1,derivative='dx')
  A_16=make_tricube_vec(1,1,1,derivative='dx')

  A_17=make_tricube_vec(0,0,0,derivative='dy')
  A_18=make_tricube_vec(1,0,0,derivative='dy')
  A_19=make_tricube_vec(0,1,0,derivative='dy')
  A_20=make_tricube_vec(1,1,0,derivative='dy')
  A_21=make_tricube_vec(0,0,1,derivative='dy')
  A_22=make_tricube_vec(1,0,1,derivative='dy')
  A_23=make_tricube_vec(0,1,1,derivative='dy')
  A_24=make_tricube_vec(1,1,1,derivative='dy')

  A_25=make_tricube_vec(0,0,0,derivative='dz')
  A_26=make_tricube_vec(1,0,0,derivative='dz')
  A_27=make_tricube_vec(0,1,0,derivative='dz')
  A_28=make_tricube_vec(1,1,0,derivative='dz')
  A_29=make_tricube_vec(0,0,1,derivative='dz')
  A_30=make_tricube_vec(1,0,1,derivative='dz')
  A_31=make_tricube_vec(0,1,1,derivative='dz')
  A_32=make_tricube_vec(1,1,1,derivative='dz')

  A_33=make_tricube_vec(0,0,0,derivative='dxy')
  A_34=make_tricube_vec(1,0,0,derivative='dxy')
  A_35=make_tricube_vec(0,1,0,derivative='dxy')
  A_36=make_tricube_vec(1,1,0,derivative='dxy')
  A_37=make_tricube_vec(0,0,1,derivative='dxy')
  A_38=make_tricube_vec(1,0,1,derivative='dxy')
  A_39=make_tricube_vec(0,1,1,derivative='dxy')
  A_40=make_tricube_vec(1,1,1,derivative='dxy')

  A_41=make_tricube_vec(0,0,0,derivative='dxz')
  A_42=make_tricube_vec(1,0,0,derivative='dxz')
  A_43=make_tricube_vec(0,1,0,derivative='dxz')
  A_44=make_tricube_vec(1,1,0,derivative='dxz')
  A_45=make_tricube_vec(0,0,1,derivative='dxz')
  A_46=make_tricube_vec(1,0,1,derivative='dxz')
  A_47=make_tricube_vec(0,1,1,derivative='dxz')
  A_48=make_tricube_vec(1,1,1,derivative='dxz')

  A_49=make_tricube_vec(0,0,0,derivative='dyz')
  A_50=make_tricube_vec(1,0,0,derivative='dyz')
  A_51=make_tricube_vec(0,1,0,derivative='dyz')
  A_52=make_tricube_vec(1,1,0,derivative='dyz')
  A_53=make_tricube_vec(0,0,1,derivative='dyz')
  A_54=make_tricube_vec(1,0,1,derivative='dyz')
  A_55=make_tricube_vec(0,1,1,derivative='dyz')
  A_56=make_tricube_vec(1,1,1,derivative='dyz')

  A_57=make_tricube_vec(0,0,0,derivative='dxyz')
  A_58=make_tricube_vec(1,0,0,derivative='dxyz')
  A_59=make_tricube_vec(0,1,0,derivative='dxyz')
  A_60=make_tricube_vec(1,1,0,derivative='dxyz')
  A_61=make_tricube_vec(0,0,1,derivative='dxyz')
  A_62=make_tricube_vec(1,0,1,derivative='dxyz')
  A_63=make_tricube_vec(0,1,1,derivative='dxyz')
  A_64=make_tricube_vec(1,1,1,derivative='dxyz')

  A_=[ [A_01],[A_02],[A_03],[A_04],[A_05],[A_06],[A_07],[A_08],                $
       [A_09],[A_10],[A_11],[A_12],[A_13],[A_14],[A_15],[A_16],                $
       [A_17],[A_18],[A_19],[A_20],[A_21],[A_22],[A_23],[A_24],                $
       [A_25],[A_26],[A_27],[A_28],[A_29],[A_30],[A_31],[A_32],                $
       [A_33],[A_34],[A_35],[A_36],[A_37],[A_38],[A_39],[A_40],                $
       [A_41],[A_42],[A_43],[A_44],[A_45],[A_46],[A_47],[A_48],                $
       [A_49],[A_50],[A_51],[A_52],[A_53],[A_54],[A_55],[A_56],                $
       [A_57],[A_58],[A_59],[A_60],[A_61],[A_62],[A_63],[A_64]  ]

  A_inv=INVERT(A_)

  IF KEYWORD_SET(debug) THEN STOP
  RETURN, A_inv
END
