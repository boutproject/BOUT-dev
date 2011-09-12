#PBS -q regular
#PBS -l mppwidth=128
#PBS -l mppnppn=4
#PBS -l walltime=5:55:00
#PBS -j eo
#PBS -V

cd $PBS_O_WORKDIR
aprun -n 128 -N 4 ./elm_pb
