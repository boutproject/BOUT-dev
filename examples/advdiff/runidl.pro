loadct,39

v=collect(path='./data', var='V')
g=file_import('./slab.grd.nc')
rxy=g.rxy
zxy=g.zxy



nt=n_elements(v[0,0,0,*])

nlev=30
level=1+findgen(nlev)/(nlev-1)

for it=0,nt-1 do begin CONTOUR, reform(v[*,*,10,it]), rxy, zxy, lev=level, chars=3, /xst, /yst, /iso & wait,1.0
