#!/bin/bash

MPIEXEC=mpirun
NP=16

# compile
make

echo Running with NP = $NP

# run case
$MPIEXEC -np $NP ./ue_bmark

