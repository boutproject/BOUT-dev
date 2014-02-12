path="data"
gridfile = path+"/cbm18_8_y064_x516_090309.pdb"

yind = 32
xind = 360

diamag = 1  ; Indicate if diamagnetic effects are included

; Fundamental quantities
mi = 2.*1.67e-27  ; Deuteron mass [kg]
e = 1.602e-19     ; electron charge [C]

gamma = 5./3.

; Import the grid file
g = file_import(gridfile)

; Get basic quantities from output files

zmin = collect(var="ZMIN", path=path)
zmax = collect(var="ZMAX", path=path)
zperiod = 1.0 / (zmax - zmin)
PRINT, 'ZPERIOD = ', zperiod

density = collect(var="density", path=path) ; Number density [m^-3]
Lbar = collect(var="Lbar", path=path) ; Length-scale [m]
Bbar = collect(var="Bbar", path=path)  ; Typical B field  [T]
Tbar = collect(var="Tbar", path=path) ; Time-scale [s]
Va   = collect(var="Va", path=path)   ; Alfven velocity [m/s]

wci = e*g.Bxy[*,yind] / mi                 ; ion cyclotron frequency [rad/s]

PRINT, "Ion cyclotron frequency: ", min(wci), " ->", max(wci), " rad/s"

cs = sqrt(gamma * g.pressure[*,yind] / (density*mi)) ; Sound speed [m/s]

PRINT, "Sound speed            : ", min(cs), " ->", max(cs), " m/s"

rhoi = cs / wci

PRINT, "Larmor radius          : ", min(rhoi), " ->", max(rhoi), " m"

qsafe = abs(g.shiftangle / (2.*!PI))

PRINT, "Safety factor          : ", min(qsafe), " ->", max(qsafe)

Raxis = (max(g.Rxy[0,*]) + min(g.Rxy[0,*]))/2.0

rminor = g.Rxy[*,yind] - Raxis

PRINT, "Minor radius           : ", min(rminor), " ->", max(rminor), " m"

; Diamagnetic drift velocity

vD = -DERIV(g.Rxy[*,yind], g.pressure[*,yind]) / (e*density*g.Bxy[*,yind]) > 0.01

PRINT, "Maximum diamagnetic velocity : ", max(vd), " m/s"

w_star = vD * (zperiod*qsafe / rminor)

PRINT, "Peak w_star (for fundamental): ", max(w_star), " rad/s"

; Get the pressure at the outboard midplane
psi = collect(var="Psi", y=yind, path=path)
t_array = collect(var="t_array", path=path)

IF diamag GT 0 THEN gamma_diamag, t_array, psi[xind,0,0,*], gamma=gamma, freq=freq 

STOP
