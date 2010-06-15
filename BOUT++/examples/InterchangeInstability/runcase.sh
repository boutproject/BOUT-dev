#!/bin/sh

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
    n ) MPIEXEC="mpirun -np";NP=$OPTARG;;
    * ) ;;   # DEFAULT
  esac
done


#-compile/build local executable
make

#-run the case       
echo Running with NP = $NP       

ln -s data_1 data
$MPIEXEC $NP ./2fluid
rm -f data

ln -s data_10 data
$MPIEXEC $NP ./2fluid
rm -f data


#-check the result
#idl runidl.pro
