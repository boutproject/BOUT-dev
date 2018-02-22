#!/bin/sh
#
# Usage example: ./runcase.sh -m "aprun -n" -n 2
#================================================#


NO_ARGS=0 
OPTERROR=65


##-setting defaults
MPIEXEC="mpirun -np "
NP=4


while getopts m:n: option; do
    case $option in
n) 
   NP=$OPTARG;
   ;;
m) 
   MPIEXEC=$OPTARG;
   ;;
  esac
done


echo "Running command: " $MPIEXEC $NP


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
idl runidl.pro

