; Given a set of normalized coordinates (xd,yd,zd) (all between 0 and 1)
; Return the vector of values SUM(i=0,3)SUM(j=0,3)SUM(k=0,3)[x^i*y^j*z^k]
; or derivatives based on the keyword selected
FUNCTION make_tricube_vec, xd, yd, zd, derivative=derivative, debug=debug

IF NOT KEYWORD_SET(derivative) THEN derivative='f'

CASE derivative OF
   'f': BEGIN
  ;form x_y_z_vec [x^i*y^j*z^k] and its derivatives
     res= [[xd^0*yd^0*zd^0],[xd^1*yd^0*zd^0],[xd^2*yd^0*zd^0],[xd^3*yd^0*zd^0],$
           [xd^0*yd^1*zd^0],[xd^1*yd^1*zd^0],[xd^2*yd^1*zd^0],[xd^3*yd^1*zd^0],$
           [xd^0*yd^2*zd^0],[xd^1*yd^2*zd^0],[xd^2*yd^2*zd^0],[xd^3*yd^2*zd^0],$
           [xd^0*yd^3*zd^0],[xd^1*yd^3*zd^0],[xd^2*yd^3*zd^0],[xd^3*yd^3*zd^0],$
           [xd^0*yd^0*zd^1],[xd^1*yd^0*zd^1],[xd^2*yd^0*zd^1],[xd^3*yd^0*zd^1],$
           [xd^0*yd^1*zd^1],[xd^1*yd^1*zd^1],[xd^2*yd^1*zd^1],[xd^3*yd^1*zd^1],$
           [xd^0*yd^2*zd^1],[xd^1*yd^2*zd^1],[xd^2*yd^2*zd^1],[xd^3*yd^2*zd^1],$
           [xd^0*yd^3*zd^1],[xd^1*yd^3*zd^1],[xd^2*yd^3*zd^1],[xd^3*yd^3*zd^1],$
           [xd^0*yd^0*zd^2],[xd^1*yd^0*zd^2],[xd^2*yd^0*zd^2],[xd^3*yd^0*zd^2],$
           [xd^0*yd^1*zd^2],[xd^1*yd^1*zd^2],[xd^2*yd^1*zd^2],[xd^3*yd^1*zd^2],$
           [xd^0*yd^2*zd^2],[xd^1*yd^2*zd^2],[xd^2*yd^2*zd^2],[xd^3*yd^2*zd^2],$
           [xd^0*yd^3*zd^2],[xd^1*yd^3*zd^2],[xd^2*yd^3*zd^2],[xd^3*yd^3*zd^2],$
           [xd^0*yd^0*zd^3],[xd^1*yd^0*zd^3],[xd^2*yd^0*zd^3],[xd^3*yd^0*zd^3],$
           [xd^0*yd^1*zd^3],[xd^1*yd^1*zd^3],[xd^2*yd^1*zd^3],[xd^3*yd^1*zd^3],$
           [xd^0*yd^2*zd^3],[xd^1*yd^2*zd^3],[xd^2*yd^2*zd^3],[xd^3*yd^2*zd^3],$
           [xd^0*yd^3*zd^3],[xd^1*yd^3*zd^3],[xd^2*yd^3*zd^3],[xd^3*yd^3*zd^3]]
        END
  'dx': BEGIN
     res= [[0],[xd^0*yd^0*zd^0],[2*xd^1*yd^0*zd^0],[3*xd^2*yd^0*zd^0],         $
           [0],[xd^0*yd^1*zd^0],[2*xd^1*yd^1*zd^0],[3*xd^2*yd^1*zd^0],         $
           [0],[xd^0*yd^2*zd^0],[2*xd^1*yd^2*zd^0],[3*xd^2*yd^2*zd^0],         $
           [0],[xd^0*yd^3*zd^0],[2*xd^1*yd^3*zd^0],[3*xd^2*yd^3*zd^0],         $
           [0],[xd^0*yd^0*zd^1],[2*xd^1*yd^0*zd^1],[3*xd^2*yd^0*zd^1],         $
           [0],[xd^0*yd^1*zd^1],[2*xd^1*yd^1*zd^1],[3*xd^2*yd^1*zd^1],         $
           [0],[xd^0*yd^2*zd^1],[2*xd^1*yd^2*zd^1],[3*xd^2*yd^2*zd^1],         $
           [0],[xd^0*yd^3*zd^1],[2*xd^1*yd^3*zd^1],[3*xd^2*yd^3*zd^1],         $
           [0],[xd^0*yd^0*zd^2],[2*xd^1*yd^0*zd^2],[3*xd^2*yd^0*zd^2],         $
           [0],[xd^0*yd^1*zd^2],[2*xd^1*yd^1*zd^2],[3*xd^2*yd^1*zd^2],         $
           [0],[xd^0*yd^2*zd^2],[2*xd^1*yd^2*zd^2],[3*xd^2*yd^2*zd^2],         $
           [0],[xd^0*yd^3*zd^2],[2*xd^1*yd^3*zd^2],[3*xd^2*yd^3*zd^2],         $
           [0],[xd^0*yd^0*zd^3],[2*xd^1*yd^0*zd^3],[3*xd^2*yd^0*zd^3],         $
           [0],[xd^0*yd^1*zd^3],[2*xd^1*yd^1*zd^3],[3*xd^2*yd^1*zd^3],         $
           [0],[xd^0*yd^2*zd^3],[2*xd^1*yd^2*zd^3],[3*xd^2*yd^2*zd^3],         $
           [0],[xd^0*yd^3*zd^3],[2*xd^1*yd^3*zd^3],[3*xd^2*yd^3*zd^3]]
        END
  'dy': BEGIN
     res= [[0],[0],[0],[0],                                                    $
           [xd^0*yd^0*zd^0],[xd^1*yd^0*zd^0],[xd^2*yd^0*zd^0],[xd^3*yd^0*zd^0],$
   [xd^0*2*yd^1*zd^0],[xd^1*2*yd^1*zd^0],[xd^2*2*yd^1*zd^0],[xd^3*2*yd^1*zd^0],$
   [xd^0*3*yd^2*zd^0],[xd^1*3*yd^2*zd^0],[xd^2*3*yd^2*zd^0],[xd^3*3*yd^2*zd^0],$
           [0],[0],[0],[0],                                                    $
           [xd^0*yd^0*zd^1],[xd^1*yd^0*zd^1],[xd^2*yd^0*zd^1],[xd^3*yd^0*zd^1],$
   [xd^0*2*yd^1*zd^1],[xd^1*2*yd^1*zd^1],[xd^2*2*yd^1*zd^1],[xd^3*2*yd^1*zd^1],$
   [xd^0*3*yd^2*zd^1],[xd^1*3*yd^2*zd^1],[xd^2*3*yd^2*zd^1],[xd^3*3*yd^2*zd^1],$
           [0],[0],[0],[0],                                                    $
           [xd^0*yd^0*zd^2],[xd^1*yd^0*zd^2],[xd^2*yd^0*zd^2],[xd^3*yd^0*zd^2],$
   [xd^0*2*yd^1*zd^2],[xd^1*2*yd^1*zd^2],[xd^2*2*yd^1*zd^2],[xd^3*2*yd^1*zd^2],$
   [xd^0*3*yd^2*zd^2],[xd^1*3*yd^2*zd^2],[xd^2*3*yd^2*zd^2],[xd^3*3*yd^2*zd^2],$
           [0],[0],[0],[0],                                                    $
           [xd^0*yd^0*zd^3],[xd^1*yd^0*zd^3],[xd^2*yd^0*zd^3],[xd^3*yd^0*zd^3],$
   [xd^0*2*yd^1*zd^3],[xd^1*2*yd^1*zd^3],[xd^2*2*yd^1*zd^3],[xd^3*2*yd^1*zd^3],$
   [xd^0*3*yd^2*zd^3],[xd^1*3*yd^2*zd^3],[xd^2*3*yd^2*zd^3],[xd^3*3*yd^2*zd^3]]
        END
  'dz': BEGIN
     res= [[0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [xd^0*yd^0*zd^0],[xd^1*yd^0*zd^0],[xd^2*yd^0*zd^0],[xd^3*yd^0*zd^0],$
           [xd^0*yd^1*zd^0],[xd^1*yd^1*zd^0],[xd^2*yd^1*zd^0],[xd^3*yd^1*zd^0],$
           [xd^0*yd^2*zd^0],[xd^1*yd^2*zd^0],[xd^2*yd^2*zd^0],[xd^3*yd^2*zd^0],$
           [xd^0*yd^3*zd^0],[xd^1*yd^3*zd^0],[xd^2*yd^3*zd^0],[xd^3*yd^3*zd^0],$
   [xd^0*yd^0*2*zd^1],[xd^1*yd^0*2*zd^1],[xd^2*yd^0*2*zd^1],[xd^3*yd^0*2*zd^1],$
   [xd^0*yd^1*2*zd^1],[xd^1*yd^1*2*zd^1],[xd^2*yd^1*2*zd^1],[xd^3*yd^1*2*zd^1],$
   [xd^0*yd^2*2*zd^1],[xd^1*yd^2*2*zd^1],[xd^2*yd^2*2*zd^1],[xd^3*yd^2*2*zd^1],$
   [xd^0*yd^3*2*zd^1],[xd^1*yd^3*2*zd^1],[xd^2*yd^3*2*zd^1],[xd^3*yd^3*2*zd^1],$
   [xd^0*yd^0*3*zd^2],[xd^1*yd^0*3*zd^2],[xd^2*yd^0*3*zd^2],[xd^3*yd^0*3*zd^2],$
   [xd^0*yd^1*3*zd^2],[xd^1*yd^1*3*zd^2],[xd^2*yd^1*3*zd^2],[xd^3*yd^1*3*zd^2],$
   [xd^0*yd^2*3*zd^2],[xd^1*yd^2*3*zd^2],[xd^2*yd^2*3*zd^2],[xd^3*yd^2*3*zd^2],$
   [xd^0*yd^3*3*zd^2],[xd^1*yd^3*3*zd^2],[xd^2*yd^3*3*zd^2],[xd^3*yd^3*3*zd^2]]
        END
 'dxy': BEGIN
     res= [[0],[0],[0],[0],                                                    $
           [0],[xd^0*yd^0*zd^0],[2*xd^1*yd^0*zd^0],[3*xd^2*yd^0*zd^0],         $
           [0],[xd^0*2*yd^1*zd^0],[2*xd^1*2*yd^1*zd^0],[3*xd^2*2*yd^1*zd^0],   $
           [0],[xd^0*3*yd^2*zd^0],[2*xd^1*3*yd^2*zd^0],[3*xd^2*3*yd^2*zd^0],   $
           [0],[0],[0],[0],                                                    $
           [0],[xd^0*yd^0*zd^1],[2*xd^1*yd^0*zd^1],[3*xd^2*yd^0*zd^1],         $
           [0],[xd^0*2*yd^1*zd^1],[2*xd^1*2*yd^1*zd^1],[3*xd^2*2*yd^1*zd^1],   $
           [0],[xd^0*3*yd^2*zd^1],[2*xd^1*3*yd^2*zd^1],[3*xd^2*3*yd^2*zd^1],   $
           [0],[0],[0],[0],                                                    $
           [0],[xd^0*yd^0*zd^2],[2*xd^1*yd^0*zd^2],[3*xd^2*yd^0*zd^2],         $
           [0],[xd^0*2*yd^1*zd^2],[2*xd^1*2*yd^1*zd^2],[3*xd^2*2*yd^1*zd^2],   $
           [0],[xd^0*3*yd^2*zd^2],[2*xd^1*3*yd^2*zd^2],[3*xd^2*3*yd^2*zd^2],   $
           [0],[0],[0],[0],                                                    $
           [0],[xd^0*yd^0*zd^3],[2*xd^1*yd^0*zd^3],[3*xd^2*yd^0*zd^3],         $
           [0],[xd^0*2*yd^1*zd^3],[2*xd^1*2*yd^1*zd^3],[3*xd^2*2*yd^1*zd^3],   $
           [0],[xd^0*3*yd^2*zd^3],[2*xd^1*3*yd^2*zd^3],[3*xd^2*3*yd^2*zd^3]]
        END
 'dxz': BEGIN
     res= [[0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [0],[xd^0*yd^0*zd^0],[2*xd^1*yd^0*zd^0],[3*xd^2*yd^0*zd^0],         $
           [0],[xd^0*yd^1*zd^0],[2*xd^1*yd^1*zd^0],[3*xd^2*yd^1*zd^0],         $
           [0],[xd^0*yd^2*zd^0],[2*xd^1*yd^2*zd^0],[3*xd^2*yd^2*zd^0],         $
           [0],[xd^0*yd^3*zd^0],[2*xd^1*yd^3*zd^0],[3*xd^2*yd^3*zd^0],         $
           [0],[xd^0*yd^0*2*zd^1],[2*xd^1*yd^0*2*zd^1],[3*xd^2*yd^0*2*zd^1],   $
           [0],[xd^0*yd^1*2*zd^1],[2*xd^1*yd^1*2*zd^1],[3*xd^2*yd^1*2*zd^1],   $
           [0],[xd^0*yd^2*2*zd^1],[2*xd^1*yd^2*2*zd^1],[3*xd^2*yd^2*2*zd^1],   $
           [0],[xd^0*yd^3*2*zd^1],[2*xd^1*yd^3*2*zd^1],[3*xd^2*yd^3*2*zd^1],   $
           [0],[xd^0*yd^0*3*zd^2],[2*xd^1*yd^0*3*zd^2],[3*xd^2*yd^0*3*zd^2],   $
           [0],[xd^0*yd^1*3*zd^2],[2*xd^1*yd^1*3*zd^2],[3*xd^2*yd^1*3*zd^2],   $
           [0],[xd^0*yd^2*3*zd^2],[2*xd^1*yd^2*3*zd^2],[3*xd^2*yd^2*3*zd^2],   $
           [0],[xd^0*yd^3*3*zd^2],[2*xd^1*yd^3*3*zd^2],[3*xd^2*yd^3*3*zd^2]]
        END
 'dyz': BEGIN
     res= [[0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [xd^0*yd^0*zd^0],[xd^1*yd^0*zd^0],[xd^2*yd^0*zd^0],[xd^3*yd^0*zd^0],$
   [xd^0*2*yd^1*zd^0],[xd^1*2*yd^1*zd^0],[xd^2*2*yd^1*zd^0],[xd^3*2*yd^1*zd^0],$
   [xd^0*3*yd^2*zd^0],[xd^1*3*yd^2*zd^0],[xd^2*3*yd^2*zd^0],[xd^3*3*yd^2*zd^0],$
           [0],[0],[0],[0],                                                    $
   [xd^0*yd^0*2*zd^1],[xd^1*yd^0*2*zd^1],[xd^2*yd^0*2*zd^1],[xd^3*yd^0*2*zd^1],$
   [xd^0*4*yd^1*zd^1],[xd^1*4*yd^1*zd^1],[xd^2*4*yd^1*zd^1],[xd^3*4*yd^1*zd^1],$
   [xd^0*6*yd^2*zd^1],[xd^1*6*yd^2*zd^1],[xd^2*6*yd^2*zd^1],[xd^3*6*yd^2*zd^1],$
           [0],[0],[0],[0],                                                    $
   [xd^0*yd^0*3*zd^2],[xd^1*yd^0*3*zd^2],[xd^2*yd^0*3*zd^2],[xd^3*yd^0*3*zd^2],$
   [xd^0*6*yd^1*zd^2],[xd^1*6*yd^1*zd^2],[xd^2*6*yd^1*zd^2],[xd^3*6*yd^1*zd^2],$
   [xd^0*9*yd^2*zd^2],[xd^1*9*yd^2*zd^2],[xd^2*9*yd^2*zd^2],[xd^3*9*yd^2*zd^2]]
        END
'dxyz': BEGIN
     res= [[0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [0],[0],[0],[0],                                                    $
           [0],[xd^0*yd^0*zd^0],[2*xd^1*yd^0*zd^0],[3*xd^2*yd^0*zd^0],         $
           [0],[xd^0*2*yd^1*zd^0],[2*xd^1*2*yd^1*zd^0],[3*xd^2*2*yd^1*zd^0],   $
           [0],[xd^0*3*yd^2*zd^0],[2*xd^1*3*yd^2*zd^0],[3*xd^2*3*yd^2*zd^0],   $
           [0],[0],[0],[0],                                                    $
           [0],[xd^0*yd^0*2*zd^1],[2*xd^1*yd^0*2*zd^1],[3*xd^2*yd^0*2*zd^1],   $
     [0],[xd^0*2*yd^1*2*zd^1],[2*xd^1*2*yd^1*2*zd^1],[3*xd^2*2*yd^1*2*zd^1],   $
     [0],[xd^0*3*yd^2*2*zd^1],[2*xd^1*3*yd^2*2*zd^1],[3*xd^2*3*yd^2*2*zd^1],   $
           [0],[0],[0],[0],                                                    $
           [0],[xd^0*yd^0*3*zd^2],[2*xd^1*yd^0*3*zd^2],[3*xd^2*yd^0*3*zd^2],   $
     [0],[xd^0*2*yd^1*3*zd^2],[2*xd^1*2*yd^1*3*zd^2],[3*xd^2*2*yd^1*3*zd^2],   $
     [0],[xd^0*3*yd^2*3*zd^2],[2*xd^1*3*yd^2*3*zd^2],[3*xd^2*3*yd^2*3*zd^2]]
        END
  else: STOP,'Incorrect selection in call to make_tricube_vec!'
ENDCASE

IF KEYWORD_SET(debug) THEN STOP
RETURN, TRANSPOSE(res)
END
