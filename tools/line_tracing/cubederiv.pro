function CubeDeriv, cellData
;
; Given 8 values on vertices of unit cubic cell
; calculate partial derivatives inside
; assuming locally linear interpolation
;
; Inputs: array [2,2,2] of 8 values on vertices of
; unit cube
;
; Outputs: partial derivatives [d/dx, d/dy, d/dz]
;;================================================;


;;-assume locally linear interpolation
;f[0,0,0], f[1,0,0], f[0,1,0], f[1,1,0], f[0,0,1], f[1,0,1], f[0,1,1], f[1,1,1]

x=double([0,1,0,1,0,1,0,1])
y=double([0,0,1,1,0,0,1,1])
z=double([0,0,0,0,1,1,1,1])

res=pdiff_xyz(x, y, z, REFORM(cellData,8))

;
;
;
return, res
end



