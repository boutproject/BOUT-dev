#!/usr/bin/env bash
# Use this script on a Cori GPU node
# Note you'll also need module load cgpu to allocate a cgpu node while on a Normal Cori node
modulepathadd() {
   if [[ "$MODULEPATH" =~ (^|:)"$1"(:|$) ]]
    then
        return 0
    fi
    export MODULEPATH=$1:$MODULEPATH
    }

module purge
module load cgpu/1.0
module load gcc/8.3.0
module load cuda/10.1.243
module load openmpi/4.0.3
module load cmake/3.14.4 
module load intel/19.1.3.304

export SPACK_ROOT=/global/project/projectdirs/bout/BOUT-GPU/spack
# the following for spack admin - add/del packages/modules
#. $SPACK_ROOT/share/spack/setup-env.sh
#spack -C /global/project/projectdirs/bout/BOUT-GPU/env module tcl refresh -y
pp=$SPACK_ROOT/share/spack/modules/cray-cnl7-skylake_avx512
modulepathadd "$pp"
echo $pp
module refresh

export cc=/usr/common/software/sles15_cgpu/gcc/8.3.0/bin/gcc
export cxx=/usr/common/software/sles15_cgpu/gcc/8.3.0/bin/g++
export CC=/usr/common/software/sles15_cgpu/gcc/8.3.0/bin/g++
export f77=/usr/common/software/sles15_cgpu/gcc/8.3.0/bin/gfortran
export fc=/usr/common/software/sles15_cgpu/gcc/8.3.0/bin/gfortran
export ftn=/usr/common/software/sles15_cgpu/gcc/8.3.0/bin/gfortran

#petsc check_cray_modules looks for the following to be set; we subvert require ProgEnv
unset CRAY_SITE_LIST_DIR

# The following for petsc without hypre
module load petsc-3.13.0-gcc-8.3.0_cgpu-ym7gwgx

module load netcdf-cxx4-4.3.1-gcc-8.3.0_cgpu-ptxvbr5
module load netcdf-c-4.7.4-gcc-8.3.0_cgpu-hpuuuxa
module load hdf5-1.10.6-gcc-8.3.0_cgpu-twbl2eg
module load fftw-3.3.8-gcc-8.3.0_cgpu-3nkroqh
module load sundials-5.1.0-gcc-8.3.0_cgpu-nxrldjq
module load superlu-dist-6.4.0-gcc-8.3.0_cgpu-gnzckur
module load openblas-0.3.12-gcc-8.3.0_cgpu-j7oii36
module load parmetis-4.0.3-gcc-8.3.0_cgpu-6gstsni
module load metis-5.1.0-gcc-8.3.0_cgpu-h3ldxba
module load autoconf-2.69-gcc-8.3.0_cgpu-y7sih6k
module load automake-1.16.2-gcc-8.3.0_cgpu-3hi6mou
module load python-3.7-gcc-8.3.0_cgpu-xkofizg
module load py-setuptools-50.3.2-gcc-8.3.0_cgpu-yhqv5q2
module load py-numpy-1.19.4-gcc-8.3.0_cgpu-r4h5h5y
module load py-netcdf4-1.4.2-gcc-8.3.0_cgpu-xtj5s5w
module load py-cython-0.29.21-gcc-8.3.0_cgpu-snty3m7
module load py-cftime-1.0.3.4-gcc-8.3.0_cgpu-e3iigwj

# allocate 10/20 cores 1/2 gpus -- note use make -j 10 or -j 20 as appropriate
alias sxd1='salloc -C gpu -N 1 -t 60 -c 10 -G 1 -A mp2'
alias sxd2='salloc -C gpu -N 1 -t 60 -c 20 -G 2 -A mp2'
export SLURM_CPUS_PER_TASK=10
module list


