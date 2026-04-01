#!/usr/bin/env bash
NP=1
FLAGS="-q -q -q -q"
EXE=tuning_regionblocksize

make

mpirun -np ${NP} ./${EXE} ${FLAGS}
