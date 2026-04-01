#!/usr/bin/env bash
        
MPIEXEC=mpirun
NP=8

#-compile/build local executable
make

#-run the case       
echo Running with NP = $NP       
$MPIEXEC -np $NP ./mhd


#-check the result
#idl runidl.pro
