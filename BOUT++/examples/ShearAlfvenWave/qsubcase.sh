#!/bin/bash

QSUB=mpich2sub # Common values are qsub, mpisub  
NP=4

#-compile/build local executable
make

#- Some queueing systems don't like name starting with number
mv 2fluid twofluid

rm -rf data*

label=( "0.15" "0.25" "0.50" "1.00" )
zmax=( 3.0246e-4 5.0410e-4 1.0082e-3 2.0164e-3 )

for ((i=0; i <= 3 ; i++))
do
  d=data_${label[$i]}
  mkdir ${d}
  sed "s/ZMAX = 3.0246e-4/ZMAX = ${zmax[$i]}/g" BOUT.inp > ${d}/BOUT.inp
  $QSUB $NP ./twofluid -d ${d}
done
