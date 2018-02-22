#!/bin/bash
#
# Create a set of runs and submit to a queueing system e.g. Sun Grid Engine
#
        
QSUB=mpich2sub # Common values are qsub, mpisub
NP=4

#-compile/build local executable
make

#- Some queueing systems don't like name starting with number
mv 2fluid twofluid

#-run the case
echo Running with NP = $NP 

rm -rf data*

zlist=( 1 2 4 8 16 32 64 128 256 )


for zval in ${zlist[@]}
do
  d=data_${zval}
  mkdir ${d}
  sed "s/Zeff = 128.0/Zeff = ${zval}.0/g" BOUT.inp > ${d}/BOUT.inp
  if [ $zval -lt 128 ]
      then
      # reduce time-step. At large times these cases produce noise
      sed "s/TIMESTEP = 5e3/TIMESTEP = 1e3/g" ${d}/BOUT.inp > ${d}/tmp
      mv -f ${d}/tmp ${d}/BOUT.inp
  fi
  
  $QSUB $NP ./twofluid -d ${d}
done

#-check the result
#idl runidl.pro
