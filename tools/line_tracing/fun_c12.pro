function fun_C12, x,y,z
;
; Right-hand side of equation d(ix)/d(iy)=C1(ix,iy,iz), d(ix)/d(iy)=C2(ix,iy,iz)
;=====================================================;

shift = 0 ; Switch on/off shifted coordinates

COMMON griddata, g, deltaZtor, Ntor
yMin=0.0
yMax=2*!PI
xMin=g.psixy[0,0]
xMax=g.psixy[g.nx-1,g.ny-1]

;;-use closest grid point for metric coefs
ix_int=FIX((g.nx-1)*(x-xMin)/(xMax-xMin))
iy_int=FIX((g.ny-1)*(y-yMin)/(yMax-yMin))
;;iz_int=FIX(iz)

;;-make sure the index does not go outside of domain
if (ix_int ge g.nx-1) then begin
    ;print, "Index too large, correcting..."
    ix_int=g.nx-1
endif

if (ix_int lt 0) then begin
    ;print, "Index too small, correcting..."
    ix_int=0
endif


rxy=g.RXY[ix_int,iy_int]
bxy=g.BXY[ix_int,iy_int]
bpxy=g.BPXY[ix_int,iy_int]
btxy=g.BTXY[ix_int,iy_int]
sinty=g.SINTY[ix_int,iy_int]
hthe=g.HTHE[ix_int,iy_int]

IF shift GT 0 THEN sinty = 0.0 ; Remove sinty from metric


;;-metric tensor components vs. radial & poloidal indices
g_sup_11=(RXY*BPXY)^2
g_sup_12=0.0d0
g_sup_13=-SINTY*(RXY*BPXY)^2
g_sup_22=1./HTHE^2
g_sup_23=-(BTXY*HTHE/(BPXY*RXY))/(HTHE^2)
g_sup_33=(SINTY*RXY*BPXY)^2 + (BXY/(RXY*BPXY))^2

;;STOP

;;-derivatives of Apar in the index space
;dApar=modelApar2(x, y, z) ;;-analytic model for testing
;dApar=modelApar(x, y, z) ;;-analytic model for testing
dApar=boutApar(x, y, z, shift=shift)

dApar_dx=dApar.x
dApar_dy=dApar.y
dApar_dz=dApar.z

;;STOP

;;-subscripted coefs
f_1 = - (1/BXY)*(g_sup_13*dApar_dx + g_sup_23*dApar_dy + g_sup_33*dApar_dz)
f_3 =   (1/BXY)*(g_sup_11*dApar_dx + g_sup_13*dApar_dz)

;;-superscripted coefficients
f_sup_1 = g_sup_11*f_1 + g_sup_13*f_3
f_sup_3 = g_sup_13*f_1 + g_sup_33*f_3

;;-rhs of dynamic equations
dx_dy = f_sup_1*HTHE/BPXY
dz_dy = f_sup_3*HTHE/BPXY

;
;
return, [dx_dy, dz_dy]
end





