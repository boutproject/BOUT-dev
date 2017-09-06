#!/bin/bash

QSUB=mpich2sub # Common values are qsub, mpisub  
NP=4


#-compile/build local executable
make

#- Some queueing systems don't like name starting with number
mv 2fluid twofluid

#-run the case       
echo Running with NP = $NP       

$QSUB $NP ./twofluid -d data_1
$QSUB $NP ./twofluid -d data_10

#-check the result
#idl runidl.pro
