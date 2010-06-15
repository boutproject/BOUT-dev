$rm -f *.so; ln -s ../../../PDB2IDL/pdb2idl.so .
.run ../../../idllib/pdb2idl.pro
.run ../../../idllib/moment_xyzt.pro
.run idlTools/get_bpp_data
.run idlTools/show_nphi
