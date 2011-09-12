#!/bin/bash

QSUB=mpich2sub # Common values are qsub, mpisub  
NP=32

#-compile/build local executable
make

$QSUB $NP ./elm_pb
