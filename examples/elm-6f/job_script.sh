#!/bin/bash
#SBATCH -N 4
#SBATCH -A mp2
#SBATCH -C cpu
#SBATCH -J elm-6f 
#SBATCH -q debug
#SBATCH -o data/output
#SBATCH -e data/error
#SBATCH -t 00:30:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=true


#run the application:
srun -n 512 ./elm_6f_v0 -d data
#srun -n 512 ./elm_6f_v0 -d data restart append
