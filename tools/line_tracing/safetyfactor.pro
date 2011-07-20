function SafetyFactor, x, debug=debug
;
; Input:  x - [weber] radial coordinate
; Output: safety factor q interpolated to x
;---------------------------------------------;

COMMON griddata, g

;;-calculate fractional index for this x
xmin=g.psixy[0,0]
xmax=g.psixy[g.nx-1,g.ny-1]
ix=(g.nx-1)*[x-xmin]/[xmax-xmin]

res=INTERPOLATE(-g.SHIFTANGLE/(2*!PI),ix)

;
;
;
if keyword_set(DEBUG) then STOP
return, res
end


