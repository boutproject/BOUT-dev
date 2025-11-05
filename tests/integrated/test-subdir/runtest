#!/bin/sh

make || exit

test ".$MPIRUN" = . && MPIRUN="mpirun -np"

$MPIRUN 1 ./subdirs
