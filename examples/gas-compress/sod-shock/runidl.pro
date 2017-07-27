
; path to the data
path = "data"

; Collect density and pressure

d = collect(var="density", path=path)
p = collect(var="pressure", path=path)
t = collect(var="t_array", path=path)

d = reform(d[2,*,0,*])
p = reform(p[2,*,0,*])

; Get the grid file
g = file_import("sod.grd.nc")

pos = findgen(g.ny)/FLOAT(g.ny-1)

showdata, d, ytitle="Density", xtitle="position"

nt = N_ELEMENTS(t)

set_plot, 'PS'
device, file="sod_result.ps"
plot, pos, d[*,nt-1], xtit="Position", ytit="Density (solid), Pressure (dashed)", title="Time = "+STR(t[nt-1])
oplot, pos, p[*,nt-1], lines=2
device, /close
set_plot, 'X'

