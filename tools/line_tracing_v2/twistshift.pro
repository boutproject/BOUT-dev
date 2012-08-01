PRO TwistShift, x1,z1, x2,z2, zShift=zShift, debug=debug
; Inputs:  coordinates of entry location 
;  x1 [weber]
;  z1 [rad]
;
; Outputs: coordinates of exit location
;  ix2 [weber]
;  iz2 [rad]
;===============================================;
  
  COMMON griddata, g, deltaZtor, Ntor

  ;;-shift in the toroidal angle after full poloidal turn
  zShift=SafetyFactor(x1)*2*!DPI
  zNew=z1-zShift                ;; Add g.Shiftangle

  ;;-cast it into the range [0,deltaZtor]
  if keyword_set(DEBUG) then print, "toroidal shift [rad]/2PI at branch-cut:", zShift/(2*!DPI)
  z2=FMODULO(zNew,deltaZtor)

  ;;-no change in the radial coordinate
  x2=x1
;
END
