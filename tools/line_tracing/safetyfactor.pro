FUNCTION SafetyFactor, x, debug=debug
; Input:  x - [weber] radial coordinate
; Output: safety factor q interpolated to x
;---------------------------------------------;

  COMMON griddata, g, deltaZtor, Ntor

  ;;-calculate fractional index for this x
  xind=INTERPOL(FINDGEN(g.nx),g.psixy[*,0],x)

  res=INTERPOLATE(-g.SHIFTANGLE/(2*!PI),xind)

  IF KEYWORD_SET(DEBUG) THEN STOP
  RETURN, res
END
