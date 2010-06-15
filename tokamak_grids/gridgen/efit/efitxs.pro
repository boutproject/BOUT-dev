;
; Demonstration of the EFITXS package for accessing EFIT data in files aeqdsk and geqdsk
;
; MVU, 3-dec-09
;----------------------------------------------------------------------;


pro efitxs, psi, r, z, plot=plot, debug=debug

COMMON efitData, efd

;;-read data from files neqdsk and aeqdsk and load data to memory
status = CALL_EXTERNAL('efitxs.so', 'read_data_g') 
status = CALL_EXTERNAL('efitxs.so', 'read_data_a') 


;;-access some selected data through accessor funcitons
nxefit=0l
nyefit=0l
status = CALL_EXTERNAL('efitxs.so', 'get_dims_g', nxefit, nyefit) 
print, "In IDL nxefit,nyefit: ", nxefit, nyefit

;;-allocate memory according to efit grid size
psi=dblarr(nxefit,nyefit)

xdim=0d0 
zdim=0d0 
rcentr=0d0 
rgrid1=0d0 
zmid=0d0 
rmagx=0d0 
zmagx=0d0
simagx=0d0 
sibdry=0d0 
bcentr=0d0 

rseps1=0d0
zseps1=0d0
rseps2=0d0
zseps2=0d0



status = CALL_EXTERNAL('efitxs.so', 'get_data_g', $
                       nxefit, nyefit, psi, xdim, zdim,$
                       rcentr, rgrid1, zmid, rmagx, zmagx,$
                       simagx, sibdry, bcentr)


status = CALL_EXTERNAL('efitxs.so', 'get_data_a', rseps1, zseps1, rseps2, zseps2)


;-reproduce EFIT grid
r=dblarr(nxefit, nyefit)
z=dblarr(nxefit, nyefit)

for i=0,nxefit-1 do begin
    for j=0,nyefit-1 do begin
        r[i,j]=rgrid1 + xdim*i/(nxefit-1)
        z[i,j]=(zmid-0.5*zdim) + zdim*j/(nyefit-1)
    endfor
endfor


;;-put all efit data to global structure efd
efd={$
psi:psi,$
r:r,$
z:z,$
xdim:xdim,$
zdim:zdim,$
rcentr:rcentr,$
rgrid1:rgrid1,$
zmid:zmid,$
rmagx:rmagx,$
zmagx:zmagx,$
simagx:simagx,$
sibdry:sibdry,$
bcentr:bcentr,$
rseps1:rseps1,$
zseps1:zseps1,$
rseps2:rseps2,$
zseps2:zseps2}




if keyword_set(PLOT) then begin
    tek_color
    contour, psi, r, z, nlev=30,/xst,/yst, /iso, tit='Psi from eqd-files'
    contour, psi, r, z, lev=[sibdry],/over, col=2
    contour, psi, r, z, lev=[simagx],/over, col=2 
    plots, rseps1*1e-2, zseps1*1e-2, psym=4, col=3
    plots, rmagx, zmagx, psym=4, col=3
endif


;
;
;
if keyword_set(DEBUG) then STOP
end
