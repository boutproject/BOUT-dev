#!/bin/bash

QSUB=mpich2sub # Common values are qsub, mpisub  
NP=16

# compile
make
mv 2fluid twofluid

echo Running with NP = $NP

# run case
$QSUB $NP ./twofluid 

