PRO FieldTrace, xIn=xIn, zIn=zIn, xout=xout, zout=zout, debug=debug
;
; Starting from point (xIn,zIn) at outer midplane make a 
; full poloidal circle and return to midplane at (xout,zout)
;
;
; Inputs: 
;    x0 [weber]
;    z0 [radian]
;
; Outputs:
;    x3 [weber]
;    z3 [radian]
;
;=======================================================;
format="('{x:',1e9.3,', y:',1e9.3,', z:',1e9.3,'}')"

;;-y-coordinate goes from 0 to 2PI, use Ny=64
Ny=64
y_imid= (0.0d0/Ny)*(2*!PI) ;-poloidal location of inner midplane
y_omid=(32.0d0/Ny)*(2*!PI) ;-poloidal location of outer midplane
y_bcut=(63.0d0/Ny)*(2*!PI) ;-poloidal location of last point below the cut

;Ny=128
;y_imid= (0.0d0/Ny)*(2*!PI) ;-poloidal location of inner midplane
;y_omid=(64.0d0/Ny)*(2*!PI) ;-poloidal location of outer midplane
;y_bcut=(127.0d0/Ny)*(2*!PI) ;-poloidal location of last point below the cut


;;-Step 1: integrate from outer midplane to last point below the cut
x0=xIn
z0=zIn
y0=y_omid
y1=y_bcut
IntegrationStep, x0, y0, z0, x1, y1, z1

if keyword_set(DEBUG) then begin
    str_In = STRING([x0,y0,z0], f=format)
    str_out = STRING([x1,y1,z1], f=format)
    print, "Step1 " + str_In + " -> " + str_out
endif


;;-Step 2, twist-shift jump between y_bcut and y_imid
;; from (x1,y1,z1) go to (x2,y2,z)
y2=y_imid
TwistShift, x1,z1, x2,z2, debug=debug

if keyword_set(DEBUG) then begin
    str_In = STRING([x1,y1,z1], f=format)
    str_out = STRING([x2,y2,z2], f=format) 
    print, "Step2 " + str_In + " -> " + str_out
endif


;;-Step 3, integrate from y2 to y3
y2=y_imid
y3=y_omid
IntegrationStep, x2, y2, z2, x3, y3, z3

if keyword_set(DEBUG) then begin
    str_In = STRING([x2,y2,z2], f=format)
    str_out = STRING([x3,y3,z3], f=format)
    print, "Step1 " + str_In + " -> " + str_out
endif

xout=x3
zout=z3
;
;
;
if keyword_set(DEBUG) then STOP
end

