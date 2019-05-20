#!/bin/bash

NO_ARGS=0 
OPTERROR=65

NP=4
test ".$MPIRUN" = . && MPIRUN="mpirun -np"

# Usage: scriptname -options
# Note: dash (-) necessary


while getopts ":n:np" Option
do
  case $Option in
    n ) NP=$OPTARG;;
    * ) ;;   # DEFAULT
  esac
done

#-compile/build local executable
make

#-run the case       
echo Running with NP = $NP       
$MPIRUN $NP ./simple_diff


