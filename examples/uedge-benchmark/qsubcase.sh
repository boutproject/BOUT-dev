#!/usr/bin/env bash

QSUB=mpich2sub # Common values are qsub, mpisub  
NP=16

# compile
make

echo Running with NP = $NP

# run case
$QSUB $NP ./ue_bmark

