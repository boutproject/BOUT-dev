#!/bin/bash

#A simple example of using the iterator test to explore scaling of the different approaches

GRID_SIZES=(4 8 16 32 64 128)
EXE=bracket
NP=1
FLAGS="-q -q -q -q performance:profileMode=true"

#Make first run
currGrid=${GRID_SIZES[0]}
mpirun -np ${NP} ./${EXE} ${FLAGS} mesh:nx=${currGrid} mesh:ny=${currGrid} mesh:nz=${currGrid}
#Do other values
for currGrid in ${GRID_SIZES[@]:1}
do
    mpirun -np ${NP} ./${EXE} ${FLAGS} mesh:nx=${currGrid} mesh:ny=${currGrid} mesh:nz=${currGrid} performance:includeHeader=false
done


#Make first run
mpirun -np ${NP} ./${EXE} ${FLAGS} 
#Do other values
for NP in 2 4
do
    mpirun -np ${NP} ./${EXE} ${FLAGS} performance:includeHeader=false
done

