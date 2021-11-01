#!/usr/bin/env bash
modulepathadd() {
   if [[ "$MODULEPATH" =~ (^|:)"$1"(:|$) ]]
    then
        return 0
    fi
    export MODULEPATH=$1:$MODULEPATH
    }

pathadd() {
   if [[ "$PATH" =~ (^|:)"$1"(:|$) ]]
    then
        return 0
    fi
    export PATH=$1:$PATH
    }

libraryadd() {
   if [[ "$LD_LIBRARY_PATH" =~ (^|:)"$1"(:|$) ]]
    then
        return 0
    fi
    export PATH=$1:$LD_LIBRARY_PATH
    }

module purge
source /opt/cray/pe/cpe/21.09/restore_lmod_system_defaults.sh
module load craype-x86-rome
module load libfabric/1.11.0.4.79
module load craype-network-ofi
module load perftools-base/21.09.0
module load xpmem/2.2.40-7.0.1.0_3.1__g1d7a24d.shasta
module load cray-dsmml/0.2.1
module load cray-libsci/21.08.1.2
module load PrgEnv-gnu
module load cpe-cuda
module load cuda/11.1.1
module load cray-pmi/6.0.13
module load cray-pmi-lib/6.0.13
module load cray-mpich/8.1.9
module load xalt/2.10.2
module load darshan/3.2.1  
module load python/3.8-anaconda-2020.11
module load py-cffi/1.14.3 
module load cray-fftw/3.3.8.11
module load cray-hdf5-parallel/1.12.0.7
module load cray-parallel-netcdf/1.12.1.7

export SPACK_ROOT=/global/cfs/cdirs/bout/BOUT-GPU/spack_perlmutter/spack
# the following for spack admin - add/del packages/modules
. $SPACK_ROOT/share/spack/setup-env.sh
#spack env activate bout-gnu
pp=$SPACK_ROOT/share/spack/modules/cray-sles15-zen2
modulepathadd "$pp"
#echo $pp
module refresh

pathadd /opt/cray/pe/mpich/8.1.9/ofi/gnu/9.1/bin/

module load cmake-3.20.5-gcc-9.3.0-6oxlyfa
module load metis-5.1.0-gcc-9.3.0-ozqdhv3
module load openblas-0.3.18-gcc-9.3.0-go4x7jk
module load parmetis-4.0.3-gcc-9.3.0-m53sma7
module load raja-0.14.0-gcc-9.3.0-es5wwa5
module load umpire-6.0.0-gcc-9.3.0-wsaedkc
module load sundials-5.8.0-gcc-9.3.0-jqainc7
module load netcdf-c-4.8.1-gcc-9.3.0-yx3q7wo
module load netcdf-cxx4-4.3.1-gcc-9.3.0-howgwk7

alias sxd='salloc -C gpu -N 1 -t 60 -A mp2_g'

module list


