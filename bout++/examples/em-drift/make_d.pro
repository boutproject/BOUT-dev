function make_d, path=path
 
if not keyword_set(PATH) then path='' 

ni=collect(path=path, var="Ni") 
phi=collect(path=path, var="phi")
ni=ni[2,*,*,*] 
phi=phi[2,*,*,*] 


du=pd_import("uedge.grd.pdb")
new_zmax=pd_read(path+"/BOUT.dmp.0.pdb", "ZMAX") ;-fraction of 2PI
rho_s = pd_read(path+"/BOUT.dmp.0.pdb", "rho_s")
wci  = pd_read(path+"/BOUT.dmp.0.pdb", "wci")
t_array = pd_read(path+"/BOUT.dmp.0.pdb", "t_array")

zeff = pd_read(path+"/BOUT.dmp.0.pdb", "Zeff")
AA = pd_read(path+"/BOUT.dmp.0.pdb", "AA")

old_zmax=new_zmax/(rho_s/du.hthe0)
print, 'new_zmax=', new_zmax, "; old_zmax=", old_zmax

d={ni_xyzt:ni, phi_xyzt:phi, rho_s:rho_s, zmax:old_zmax, Zeff:zeff, AA:AA, t_array:t_array, wci:wci}

return,d
end
