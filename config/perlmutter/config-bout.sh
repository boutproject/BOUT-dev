#!/usr/bin/env bash
# this script is designed to be run in the BOUT source dir
arch=$(uname -m)
compiler=gcc
scratch_dir=`pwd`
build_prefix=${scratch_dir}/build/${arch}-${compiler}/
install_prefix=${scratch_dir}/install/${arch}-${compiler}/
source_prefix=${scratch_dir}/
tpl_prefix=/global/cfs/cdirs/bout/BOUT-GPU/tpl_11p/
tpl_install_prefix=${tpl_prefix}/install/${arch}-${compiler}/

pkg=BOUT-dev
echo 'continue with script ' $pkg
if [[ "$1" == "-c" ]]; then
    pkg=$2
    dir=${build_prefix}/$pkg
    echo Removing $dir ...
    rm -rf $dir
fi

if [[ "$compiler" == "xl" ]]; then
    cc=xlc
    cpp=xlC
    fc=xlf
    extra_cflags=""
elif [[ "$compiler" == "clang" ]]; then
    cc=clang
    cpp=clang++
#    fc=xlf
    fc=mpifort
    extra_cflags=""
elif [[ "$compiler" == "gcc" ]]; then
    cc=mpicc
    cpp=mpicxx
    fc=mpifort
    extra_cflags=""
elif [[ "$compiler" == "CC" ]]; then
    cc=cc
    cpp=CC
    fc=ftn
    extra_cflags=""
else
    exit 1
fi

build_dir=${build_prefix}/${pkg}
install_dir=${install_prefix}/${pkg}
source_dir=${source_prefix} 
echo 'Build in ' ${build_dir}
echo 'Install in ' ${install_dir}
echo 'Source in ' ${source_dir}

if [ "$pkg" == "BOUT-dev" ]; then
    echo 'enter BOUT-dev script'
    source_dir=${source_prefix} 
    mkdir -p $build_dir && cd $build_dir
    cmake  \
          -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_CXX_STANDARD=17 \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_CXX_FLAGS="-w " \
          -DCMAKE_CUDA_FLAGS="-w " \
          -DBOUT_USE_NETCDF=On \
          -DBOUT_USE_FFTW=On \
          -DBOUT_USE_LAPACK=On \
          -DBOUT_USE_NLS=On \
          -DBOUT_USE_PETSC=Off \
          -DBOUT_USE_PVODE=On \
          -DBOUT_USE_SUNDIALS=On \
          -DENABLE_GTEST_DEATH_TESTS=On \
          -DBOUT_ENABLE_RAJA=On \
          -DBOUT_ENABLE_UMPIRE=On \
          -DBOUT_ENABLE_MPI=On \
          -DBOUT_ENABLE_OPENMP=Off \
          -DBOUT_ENABLE_CUDA=On \
          -DCMAKE_CUDA_STANDARD=17 \
          -DCUDA_ARCH="compute_80,code=sm_80" \
          -DBOUT_USE_HYPRE=On \
          -DHYPRE_DIR="${tpl_prefix}/hypre_dir/hypre_autoconf/install" \
          -DHYPRE_CUDA=On \
          -DBUILD_SHARED_LIBS=Off \
          -DBUILD_TESTING=On \
          -DCHECK=1 \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=On \
          -DCMAKE_VERBOSE_MAKEFILE=Off \
          -Dgtest_disable_pthreads=ON \
          $source_dir
fi
