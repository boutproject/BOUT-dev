#!/usr/bin/env bash

MPIEXEC=mpirun
NP=2

rm -rf data*

label=( "0.15" "0.25" "0.50" "1.00" )
zmax=( 3.0246e-4 5.0410e-4 1.0082e-3 2.0164e-3 )

for ((i=0; i <= 3 ; i++))
do
  mkdir data_${label[$i]}
  ln -s data_${label[$i]} data
  sed "s/ZMAX = 3.0246e-4/ZMAX = ${zmax[$i]}/g" BOUT.inp > data/BOUT.inp
  $MPIEXEC -np $NP ./2fluid
  rm -f data
done
