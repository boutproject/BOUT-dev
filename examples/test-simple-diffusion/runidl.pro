;Make sure BOUT_TOP is set, e.g.,
;export BOUT_TOP=$HOME/bout/trunk/BOUT++/
;
;If using PACT make sure pdb2idl.so is set, e.g.,
;rm -f pdb2idl.so; ln -s $BOUT_TOP/../PDB2IDL/pdb2idl.so
;========================================================;

.run $BOUT_TOP/../idllib/pdb2idl.pro
SPAWN, "printenv BOUT_TOP", BOUT_TOP
!path=!path+":"+BOUT_TOP+"/../idllib"


path="./"
fdu = path+'data/advect.grd.cdl'
du=file_import(fdu)

path=path+"/data"
rho = collect(var="density", path=path)


sz=size(rho)

nx=sz[1]
ny=sz[2]
nz=sz[3]
nt=sz[4]


tit='1D Advection Test'
for it=0,nt-1 do begin plot, rho[nx/2,*,nz/2-1,it], yr=[0,2],/xst,/yst, tit=tit & wait,0.1
