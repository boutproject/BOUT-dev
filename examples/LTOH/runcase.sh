#!/bin/bash

NO_ARGS=0 
OPTERROR=65

if [ $# -eq "$NO_ARGS" ]  # Script invoked with no command-line args?
then
	MPIEXEC=
  NP=
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

echo Compiling

# compile
make

echo Running with NP = $NP

# run case
$MPIEXEC $NP ./elm_pb

