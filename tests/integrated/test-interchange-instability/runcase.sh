#!/bin/sh
#
# Usage example: ./runcase.sh -m "aprun -n" -n 2
#================================================#


NO_ARGS=0 
OPTERROR=65


##-setting defaults
test ".$MPIRUN" = . && MPIRUN="mpirun -np"
NP=4


while getopts m:n: option; do
    case $option in
n) 
   NP=$OPTARG;
   ;;
m) 
   MPIRUN=$OPTARG;
   ;;
  esac
done


echo "Running command: " $MPIRUN $NP


#-compile/build local executable
make

#-run the cases
echo Running with NP = $NP       

$MPIRUN $NP ./2fluid -d data_1
rm -f data

$MPIRUN $NP ./2fluid -d data_10
rm -f data


#-check the result
idl runidl.pro
