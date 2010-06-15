;
; Read EFIT data
;----------------------------------------------------------------------;


pro read_efit, psi, r, z, load=load, debug=debug
;
;
;
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


;;-grid range, [m]
rmin=rgrid1
rmax=rmin+xdim
zmin=(zmid-0.5*zdim)
zmax=zmin+zdim

if (1) then begin
    ;;-make the bottom at zero, like UEDGE does
    zmax=zmax-zmin
    zmin=0d0
endif


;-reproduce EFIT grid
r=dblarr(nxefit, nyefit)
z=dblarr(nxefit, nyefit)

for i=0,nxefit-1 do begin
    for j=0,nyefit-1 do begin
        r[i,j]=rmin + (rmax-rmin)*i/(nxefit-1)
        z[i,j]=zmin + (zmax-zmin)*j/(nyefit-1)
    endfor
endfor



if keyword_set(LOAD) then begin
    print, 'Loading psi data...'
    Gint_Load, psi, nxefit, rmin, rmax, zmin, zmax, rcentr, bcentr
endif


;
;
;
if keyword_set(DEBUG) then STOP
end
