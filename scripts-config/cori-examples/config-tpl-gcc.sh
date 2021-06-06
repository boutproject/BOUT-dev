#!/usr/bin/env bash
arch=$(uname -m)
compiler=gcc
base=/global/project/projectdirs/bout/BOUT-GPU/tpl
build_prefix=${base}/build/${arch}-${compiler}/
install_prefix=${base}/install/${arch}-${compiler}/
source_prefix=${base}

#for automake version of hypre be sure to export CC and CXX environment variables as hypre will pick up XL for the default
# export CC=mpicc
# export CXX=mpiCC
# ./configure --prefix=~/workspace/hypre_automake/install --with-cuda --enable-unified-memory --enable-debug --enable-cublas HYPRE_CUDA_SM=70

pkg=$1
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
    cc=mpiclang
    cpp=mpiclang++
    fc=mpigfortran
    #fc=mpifort
    extra_cflags=""
elif [[ "$compiler" == "gcc" ]]; then
    cc=mpicc
    cpp=mpiCC
    fc=mpigfortran
    extra_cflags=""
else
    exit 1
fi

build_dir=${build_prefix}/${pkg}
install_dir=${install_prefix}/${pkg}
if [ "$pkg" == "raja" ]; then
    source_dir=${source_prefix}/${pkg}
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DRAJA_CXX_STANDARD_FLAG="default" \
          -DCMAKE_BUILD_TYPE=Release\
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DENABLE_OPENMP=Off \
          -DENABLE_CUDA=On \
          -DCUDA_ARCH=sm_70 \
          -DENABLE_TESTS=ON \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=On \
          -DCMAKE_VERBOSE_MAKEFILE=On \
          $source_dir
elif [ "$pkg" == "umpire" ]; then
    source_dir=${source_prefix}/Umpire
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_CXX_STANDARD=14 \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=Release \
          -DENABLE_CUDA=On \
          -DCMAKE_CUDA_FLAGS="-arch sm_70 --expt-extended-lambda --expt-relaxed-constexpr" \
          -DENABLE_TESTS=On \
          -DENABLE_TOOLS=On \
          -DENABLE_NUMA=On \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=On \
          -DCMAKE_VERBOSE_MAKEFILE=On \
          $source_dir
fi
