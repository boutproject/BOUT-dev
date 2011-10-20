#!/bin/bash

NO_ARGS=0 
OPTERROR=65

if [ $# -eq "$NO_ARGS" ]  # Script invoked with no command-line args?
then
    MPIEXEC="mpirun -np "
    NP=4
fi  
# Usage: scriptname -options
# Note: dash (-) necessary


while getopts ":n:np" Option
do
  case $Option in
    n ) MPIEXEC="mpiexec -np";NP=$OPTARG;;
    * ) ;;   # DEFAULT
  esac
done

#-compile/build local executable
make

#-run the case       
echo Running with NP = $NP       
$MPIEXEC $NP ./simple_diff


