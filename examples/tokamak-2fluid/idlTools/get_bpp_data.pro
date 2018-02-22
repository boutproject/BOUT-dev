pro get_bpp_data, arms=arms, du=du, t=t, wci=wci, ni=ni, phi=phi, rho=rho, jpar=jpar, debug=debug, file=file
;
;
;

du=pd_import("data/uedge.grd.pdb")

d = collect(path="data", var="Ni")
moment_xyzt, d, rms=rms_ni

d = collect(path="data", var="phi")
moment_xyzt, d, rms=rms_phi

d = collect(path="data", var="rho")
moment_xyzt, d, rms=rms_rho

d = collect(path="data", var="jpar")
moment_xyzt, d, rms=rms_jpar


arms={$
       ni_xyzt:rms_ni[*,*,1:*],$
       phi_xyzt:rms_phi[*,*,1:*],$
       rho_xyzt:rms_rho[*,*,1:*],$
       jpar_xyzt:rms_jpar[*,*,1:*]}


tt = collect(path="data", var="t_array")
wci = collect(path="data", var="wci")
t = tt[1:*] / wci

;
;
;
if keyword_set(DEBUG) then STOP
end
