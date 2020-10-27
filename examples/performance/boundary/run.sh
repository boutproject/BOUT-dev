#!/bin/bash

for ord in o2 o3 o4
do
    for bndry in dirichlet_ neumann_ free_ dirichlet_nu_ neumann_nu_ free_nu_
    do
	out=$(./boundary -q -q -q f:bndry_all=${bndry}$ord 2>&1 ntests=10000) && echo "$out" || echo ${bndry}$ord failed
    done
done
