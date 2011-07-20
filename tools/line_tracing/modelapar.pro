function modelApar, x, y, z, deriv=deriv, debug=debug
;
; Input: fractional indices
; Output: {dApar/d(x), dApar/d(y), dApar/d(z)}
;====================================================;

xmid=-0.1 ;;-better get it from the structure g
ymid=!PI

alpx=1e-2
alpy=2e-2

;-if BOUT domain iz 2PI/5 then in the perturbation
;only nz=5,10,15,.. etc can be present
nz=15d0 

;;;Apar=SIN(nz*z)*EXP(-alpx*(x-xMid)^2)*EXP(-alpy*(y-yMid)^2)
dAdx=   SIN(nz*z)*EXP(-alpx*(x-xMid)^2)*EXP(-alpy*(y-yMid)^2)*(-2*(x-xMid)*alpx)
dAdy=   SIN(nz*z)*EXP(-alpx*(x-xMid)^2)*EXP(-alpy*(y-yMid)^2)*(-2*(y-yMid)*alpy)
dAdz=nz*COS(nz*z)*EXP(-alpx*(x-xMid)^2)*EXP(-alpy*(y-yMid)^2)


sFactor=1e-6
res={x:sFactor*dAdx, y:sFactor*dAdz, z:sFactor*dAdz}

;
;
;
if keyword_set(DEBUG) then STOP
return, res
end




