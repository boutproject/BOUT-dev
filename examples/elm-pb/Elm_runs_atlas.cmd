 #!/bin/csh
# Sample MOAB script to be submitted with msub
#MSUB -V 
#MSUB -q pbatch                # pool to use
#MSUB -A gyroturb
#MSUB -N den8S5R7 # sets job name (limit of 7 characters)
#MSUB -e xboutn15.err # sets output log name
#MSUB -o xboutn15.log # sets output log name
#MSUB -m abe
#MSUB -l partition=atlas       # machine to run on
#MSUB -l walltime=15:00:00  # sets maximum total CPU time
#MSUB -l nodes=16:ppn=8

##MSUB -l depend=34876

cd /p/lscratchb/xu/bout_Hall_hyperResist/elm_runs/phi0_diamag_Vt0_0source_0viscos_cbm18_dens8/V_VE0_VD0_hyperR/Nonlinear_n15_lund1e5

srun -n128 -p pbatch /p/lscratcha/xu/bout/trunk/BOUT++/examples/highbeta_peeling_ballooning/elm_pb

date
echo "ALL DONE"
