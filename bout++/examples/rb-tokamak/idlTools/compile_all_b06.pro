$rm *.so
$ln -s $HOME/BoutDataAnalysis/PDB2IDL/pdb2idl.so .

pwd=getenv('PWD')
cd,'$HOME/BoutDataAnalysis/'

.run PDB2IDL/pdb2idl.pro
.run PDB2IDL/reformat.pro
.run Plots2D/moment_xyzt.pro
.run Plots2D/allrms.pro

.run UEgrid/curvature.pro 
.run GRD2PDB/grd2pdb
.run UEgrid/read_uedgegrd

.run UEgrid/read_uedata.pro 
.run UEgrid/read_uedata3.pro 

cd,pwd
