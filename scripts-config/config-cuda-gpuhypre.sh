#!/usr/bin/env bash
arch=$(uname -m)
compiler=gcc
scratch_dir=/p/gpfs1/fisher47/bout/BOUT_build_cuda/
build_prefix=${scratch_dir}/build/${arch}-${compiler}/
install_prefix=${scratch_dir}/install/${arch}-${compiler}/
source_prefix=/p/gpfs1/fisher47/bout/


# The following modules are provided by the system
module --ignore-cache load spectrum-mpi/rolling-release
module --ignore-cache load cmake/3.14.5
module --ignore-cache load cuda/10.1.243
module --ignore-cache load gcc/8.3.1
module --ignore-cache load sundials/4.1.0
module --ignore-cache load lapack/3.8.0-gcc-4.9.3

# The following modules are installed by you using spack
module --ignore-cache load fftw-3.3.8-gcc-8.3.1-vlusxnt
module --ignore-cache load hdf5-1.10.1-gcc-8.3.1-rd54eh6
module --ignore-cache load python-3.7.7-gcc-8.3.1-4vwg3n5
module --ignore-cache load netcdf-c-4.7.3-gcc-8.3.1-fpxr24s 
module --ignore-cache load netcdf-cxx4-4.3.1-gcc-8.3.1-ybne6oo
module --ignore-cache load petsc-3.12.3-gcc-8.3.1-triqq62
module --ignore-cache load py-cftime-1.0.3.4-gcc-8.3.1-bbyzoow
module --ignore-cache load py-cython-0.29.16-gcc-8.3.1-u4srpzl
module --ignore-cache load py-pybind11-2.5.0-gcc-8.3.1-tyma4yw
module --ignore-cache load py-numpy-1.18.3-gcc-8.3.1-uh6npig
module --ignore-cache load py-scipy-1.4.1-gcc-8.3.1-4qmxxkv
module --ignore-cache load py-netcdf4-1.4.2-gcc-8.3.1-ujcrt7k


# PETSc spack install signature
#spack install petsc@3.12.3%gcc@8.3.1 +fftw +metis +superlu-dist ~hypre +mpi ^hdf5@1.10.1+cxx+hl+mpi+pic+shared ^spectrum-mpi@rolling-release%gcc@8.3.1 ^metis%gcc@8.3.1+real64 ^fftw%gcc@8.3.1

#netcdf-cxx4 signature
#spack install netcdf-cxx4 ^spectrum-mpi@rolling-release%gcc@8.3.1 ^hdf5@1.10.1+cxx+hl+mpi+pic+shared 


# Note spectrum-mpi is picked up by an edit you provide to .spack/linux/packages.yaml which for gcc8.3.1 is
#packages:
#  all:
#    providers:
#       mpi: [spectrum-mpi, mvapich2, openmpi]
#  spectrum-mpi:
#    buildable: False
#    paths:
#      spectrum-mpi@rolling-release%gcc@8.3.1 arch=linux-rhel7-ppc64le : /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1


# use module display package-name to show what package provides including path info

# be sure to provide  ^spectrum-mpi@rolling-release%gcc@8.3.1 for packages needing mpi
# else you may receive python error '<' not supported between instances of 'str' and 'NoneType'

#for PETSC be sure to module load in your .profile to avoid system ambiguities
# or before typing make -j

#for automake version of hypre be sure to export CC and CXX environment variables as hypre will pick up XL for the default
# export CC=mpicc
# export CXX=mpiCC
# ./configure --prefix=/g/g0/holger/workspace/hypre_automake/install --with-cuda --enable-unified-memory --enable-debug --enable-cublas HYPRE_CUDA_SM=70

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
    #fc=mpifort
    extra_cflags=""
else
    exit 1
fi

echo 'continue with script'

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
    cc=clang
    cpp=clang++
#    fc=xlf
    fc=mpifort
    extra_cflags=""
elif [[ "$compiler" == "gcc" ]]; then
    cc=mpicc
    cpp=mpiCC
    fc=mpifort
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
          -DCMAKE_BUILD_TYPE=RelWithDebInfo \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DENABLE_OPENMP=On \
          -DENABLE_CUDA=On \
          -DENABLE_TESTS=ON \
          $source_dir
elif [ "$pkg" == "umpire" ]; then
    source_dir=${source_prefix}/Umpire
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=RelWithDebInfo \
          -DENABLE_CUDA=On \
          -DENABLE_TESTS=On \
          -DENABLE_TOOLS=On \
          -DENABLE_NUMA=On \
          $source_dir
elif [ "$pkg" == "hypre" ]; then
    source_dir=${source_prefix}/${pkg}/src
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_BUILD_TYPE=RelWithDebInfo \
          -DCMAKE_PREFIX_PATH="${install_prefix}/raja/share/raja/cmake;${install_prefix}/umpire/share/umpire/cmake" \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DHYPRE_BUILD_TESTS=On \
          -DENABLE_GTEST_DEATH_TESTS=On \
          -DENABLE_RAJA=Off \
          -DENABLE_MPI=On \
          -DENABLE_OPENMP=Off \
          -DENABLE_CUDA=On \
          -DCUDA_ARCH=sm_70 \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=On \
          -DCMAKE_VERBOSE_MAKEFILE=On \
          -Dgtest_disable_pthreads=ON \
          $source_dir
elif [ "$pkg" == "BOUT-dev" ]; then
    source_dir=${source_prefix}/${pkg}
    cd $source_dir 
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=RelWithDebInfo \
          -DCMAKE_PREFIX_PATH="${install_prefix}/raja/share/raja/cmake;${install_prefix}/HDF5/share/cmake;${install_prefix}/umpire/share/umpire/cmake" \
          -DUSE_NETCDF=On \
          -DNCXX4_CONFIG:FILEPATH=/var/tmp/fisher47/spack-stage/spack-stage-netcdf-cxx4-4.3.1-nh4qne3jjk77a4tmok2xbi57gjhfiohq/bin/ncxx4-config \
          -DNC_CONFIG:FILEPATH=/var/tmp/fisher47/spack-stage/spack-stage-netcdf-c-4.7.3-y7mfxm7lkqk76k46micqosfnv5oyup4l/bin/nc-config \
          -DUSE_FFTW=On \
          -DUSE_LAPACK=On \
          -DUSE_NLS=On \
          -DENABLE_PETSC=On \
          -DPETSC_DIR=/var/tmp/fisher47/spack-stage/spack-stage-petsc-3.13.0-4273anzwyd5cqhkx6fceninoycg4wo66 \
          -DUSE_PVODE=On \
          -DUSE_SUNDIALS=On \
          -DENABLE_GTEST_DEATH_TESTS=On \
          -DENABLE_RAJA=Off \
          -DENABLE_MPI=On \
          -DENABLE_OPENMP=Off \
          -DENABLE_CUDA=On \
          -DENABLE_HYPRE=On \
          -DHYPRE_DIR="/p/gpfs1/fisher47/bout/hypre_gpu/install" \
          -DBUILD_SHARED_LIBS=Off \
          -DENABLE_GTEST=On \
          -DENABLE_GMOCK=On \
          -DBUILD_GTEST=Off \
          -DBUILD_TESTING=On \
          -DENABLE_TESTS=On \
          -DENABLE_EXAMPLES=OFF \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=On \
          -DCMAKE_VERBOSE_MAKEFILE=On \
          -Dgtest_disable_pthreads=ON \
          -DCMAKE_INSTALL_RPATH="/var/tmp/fisher47/spack-stage/spack-stage-petsc-3.13.0-4273anzwyd5cqhkx6fceninoycg4wo66/lib;/var/tmp/fisher47/spack-stage/spack-stage-hdf5-1.10.6-6qtujf5zc5iwwyow6ub52jfxyzej5gse/lib" \
          -DCMAKE_BUILD_RPATH="/var/tmp/fisher47/spack-stage/spack-stage-petsc-3.13.0-4273anzwyd5cqhkx6fceninoycg4wo66/lib;/var/tmp/fisher47/spack-stage/spack-stage-hdf5-1.10.6-6qtujf5zc5iwwyow6ub52jfxyzej5gse/lib" \
          $source_dir
fi
