path="."


fdu = path+'/uedge.grd.pdb'

path=path+"/data"

du=file_import(fdu)

ni = collect(var="Ni", path=path)
rho = collect(var="rho", path=path)
phi = collect(var="phi", path=path)
jpar = collect(var="jpar", path=path)

t_array = collect(var="t_array", path=path)
wci = collect(var="wci", path=path)
rho_s = collect(var="rho_s", path=path)

; Remove the first time-point (as BOUT-06)
ni = ni[*,*,*,1:*]
rho = rho[*,*,*,1:*]
phi = phi[*,*,*,1:*]
jpar=jpar[*,*,*,1:*]
t_array = t_array[1:*]

trange = N_ELEMENTS(t_array)

; Put into a structure to mimic BOUT_FLUC_xyzt
d = {ni_xyzt:ni, rho_xyzt:rho, phi_xyzt:phi, jpar_xyzt:jpar, $
     t_array:t_array, wci:wci, rho_s:rho_s, trange:trange}

fit_time, d,du

