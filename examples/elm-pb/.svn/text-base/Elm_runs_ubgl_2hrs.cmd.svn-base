#!/bin/csh
# Sample MOAB script to be submitted with msub
#MSUB -V
#MSUB -q pdebug                # pool to use
#MSUB -A gyroturb
#MSUB -N xdivS8R4 # sets job name (limit of 7 characters)
#MSUB -e xboutn15.err # sets output log name
#MSUB -o xboutn15.log # sets output log name
#MSUB -m abe
#MSUB -l partition=ubgl       # machine to run on
#MSUB -l walltime=00:30:00  # sets maximum total CPU time
#MSUB -l nodes=128

##MSUB -l depend=34876

cd /p/gscratchc/xu/iaea2010/bout/trunk/BOUT++/examples/highbeta_peeling_ballooning

mpirun -verbose 1 -mode CO -exe /p/gscratchc/xu/iaea2010/bout/trunk/BOUT++/examples/highbeta_peeling_ballooning/elm_pb  -cwd /p/gscratchc/xu/iaea2010/bout/trunk/BOUT++/examples/highbeta_peeling_ballooning

 date
 echo "ALL DONE"
