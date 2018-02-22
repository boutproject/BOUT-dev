; Perpendicular Laplacian test
; 
; Tests the accuracy of the Delp2 operator on both uniform and
; non-uniform grids. The same coefficients are used in both 
; Delp2 and Laplacian inversion.
;
; Usage:
; =====
; 
; $ idl run_idl.pro
;

path = "data"

input = collect(var="input", path=path)
reference = collect(var="reference", path=path)
output = collect(var="result", path=path)

s = SIZE(input, /dim)
IF N_ELEMENTS(s) LT 3 THEN nz = 1 ELSE nz = s[2]

grid = file_import("test_delp2.grd.nc")

PRINT, "Mesh size: "+STR(grid.nx)+" x "+STR(grid.ny)+" x "+STR(nz)

sym=[1,2,4,6]

PRINT, ""
PRINT, "Meshes"
PRINT, "======"
PRINT, "Solid line: Analytic result"
PRINT, "Crosses: Uniform mesh"
PRINT, "Stars and dotted line: Linearly changing mesh"
PRINT, "Diamonds and dashed line: Quadratically changing mesh"
PRINT, "Squares and dot-dashed line: mesh with a jump in dx"
PRINT, ""

FOR t=0,3 DO BEGIN          &$
  inds = [0,4,8,12]+t       &$
  data = output[*,inds,*]   &$
  xpos = grid.xpos[*,inds]  &$
  ref = reference[*,inds,*] &$
  PLOT, xpos[*,0], ref[*,0], yr=[MIN(data), MAX(data)], $
     title="TEST "+STR(t) &$
  FOR i=0,3 DO BEGIN &$
    OPLOT, xpos[*,i], data[*,i], lines=i &$
    OPLOT, xpos[*,i], data[*,i], psym=sym[i] &$
  ENDFOR &$
  CURSOR, a,b,/down &$
ENDFOR
