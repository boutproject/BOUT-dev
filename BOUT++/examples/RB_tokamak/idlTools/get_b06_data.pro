pro get_b06_data, arms=arms, du=du, t=t, wci=wci, ni=ni, phi=phi, rho=rho, jpar=jpar, debug=debug, file=file
;
;
;

if not keyword_set(FILE) then file='BOUT_FLUC_xyzt.pdb'
print, "Restoring data from ", file
d=pd_import(file)

gridFile="uedge.grd.pdb"
print, "Importing BOUT-06 grid from ", gridfile
du=pd_import(gridFile)

wci=d.wci
t=d.t_array/wci ;;-[s]

ni  =d.ni_xyzt[*,*,0:d.ngz-2,*]
phi =d.phi_xyzt[*,*,0:d.ngz-2,*]
rho =d.rho_xyzt[*,*,0:d.ngz-2,*]
jpar=d.jpar_xyzt[*,*,0:d.ngz-2,*]


;if keyword_set(ARMS) then begin
    MOMENT_XYZT, ni,   rms=rms_ni
    MOMENT_XYZT, phi,  rms=rms_phi
    MOMENT_XYZT, rho,  rms=rms_rho
    MOMENT_XYZT, jpar, rms=rms_jpar


    arms={$
           ni_xyzt:rms_ni,$
           phi_xyzt:rms_phi,$
           rho_xyzt:rms_rho,$
           jpar_xyzt:rms_jpar}
;endif

;
;
;
if keyword_set(DEBUG) then STOP
end
