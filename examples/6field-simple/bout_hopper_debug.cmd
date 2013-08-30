#PBS -q debug
#PBS -l mppwidth=128
#PBS -l mppnppn=24
#PBS -l walltime=0:10:00
#PBS -N my_job
#PBS -e my_job.$PBS_JOBID.err
#PBS -e my_job.$PBS_JOBID.out
#PBS -V

cd $PBS_O_WORKDIR
aprun -n 128 -N 24 ./elm_6f

