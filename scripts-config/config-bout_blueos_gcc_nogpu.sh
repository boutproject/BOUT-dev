#!/usr/bin/env bash


arch=$(uname -m)
compiler=gcc
scratch_dir=${HOME}/workspace/BOUT_build_serial/
build_prefix=${scratch_dir}/build/${arch}-${compiler}/
install_prefix=${scratch_dir}/install/${arch}-${compiler}/
source_prefix=${HOME}/workspace/

module --ignore-cache load spectrum-mpi/rolling-release
module --ignore-cache load cmake/3.14.5
module --ignore-cache load cuda/10.1.243
module --ignore-cache load gcc/8.3.1
module --ignore-cache load python/3.7.2
module --ignore-cache load netcdf/4.7.0
module --ignore-cache load sundials/4.1.0
module --ignore-cache load fftw-3.3.8-gcc-8.3.1-vlusxnt
module --ignore-cache load lapack/3.8.0-gcc-4.9.3
module --ignore-cache load python-3.7.6-gcc-8.3.1-usivcqa
module --ignore-cache load petsc-3.12.3-gcc-8.3.1-jnsim2o
module --ignore-cache load py-setuptools-41.4.0-gcc-8.3.1-d4wih3g
module --ignore-cache load py-cftime-1.0.3.4-gcc-8.3.1-q6ofwn4
module --ignore-cache load py-cython-0.29.14-gcc-8.3.1-5sfsoak
module --ignore-cache load py-pybind11-2.5.0-gcc-8.3.1-4hcy5vc
module --ignore-cache load py-numpy-1.18.2-gcc-8.3.1-6wn32qx
module --ignore-cache load py-scipy-1.4.1-gcc-8.3.1-cck6efe


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

#          -DOpenMP_Fortran_FLAGS="-qsmp" \
#          --debug-trycompile \
#          -DRAJA_CUDA_ARCH=sm_60 \
#          -DPETSC_DIR=${HOME}/spack/opt/spack/linux-linuxmint19-haswell/gcc-7.3.0/petsc-3.12.3-s77fwdjblxkslbp2ymv5uehkmimdbyka \
#          -DSUNDIALS_ROOT=${HOME}/spack/opt/spack/linux-linuxmint19-haswell/gcc-7.3.0/sundials-5.1.0-jawjhwlqup3qlusyygjoj4r4xxw7qvlf \

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
          -DENABLE_CUDA=Off \
          -DENABLE_TESTS=ON \
          $source_dir
elif [ "$pkg" == "HDF5" ]; then
    source_dir=${source_prefix}/hdf5-1.10.1
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=Release \
          -DHDF5_BUILD_EXAMPLES=OFF \
          $source_dir
elif [ "$pkg" == "umpire" ]; then
    source_dir=${source_prefix}/Umpire
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=RelWithDebInfo \
          -DENABLE_CUDA=Off \
          -DENABLE_TESTS=On \
          -DENABLE_TOOLS=On \
          -DENABLE_NUMA=On \
          $source_dir
elif [ "$pkg" == "hypre" ]; then
    source_dir=${source_prefix}/${pkg}/src
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_BUILD_TYPE=RelWithDebInfo \
          -DCMAKE_PREFIX_PATH="${install_prefix}/raja/share/raja/cmake;${install_prefix}/umpire/share/umpire/cmake" \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DHYPRE_BUILD_TESTS=On \
          -DENABLE_GTEST_DEATH_TESTS=On \
          -DENABLE_RAJA=Off \
          -DENABLE_MPI=On \
          -DENABLE_OPENMP=Off \
          -DENABLE_CUDA=Off \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=On \
          -DCMAKE_VERBOSE_MAKEFILE=On \
          -Dgtest_disable_pthreads=ON \
          $source_dir
elif [ "$pkg" == "BOUT-dev" ]; then
    source_dir=${source_prefix}/${pkg}
    cd $source_dir && {
        branch=$(git rev-parse --abbrev-ref HEAD)
        [[ "$branch" != "hypre-laplacexy-blt" ]] && exit 1
    }
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=Debug \
          -DCMAKE_PREFIX_PATH="${install_prefix}/raja/share/raja/cmake;${install_prefix}/HDF5/share/cmake;${install_prefix}/umpire/share/umpire/cmake" \
          -DUSE_NETCDF=On \
          -DNCXX4_CONFIG:FILEPATH=/usr/tcetmp/packages/petsc/netcdf-4.7.0/bin/ncxx4-config \
          -DNC_CONFIG:FILEPATH=/usr/tcetmp/packages/petsc/netcdf-4.7.0/bin/nc-config \
          -DUSE_FFTW=On \
          -DUSE_LAPACK=On \
          -DUSE_NLS=On \
          -DENABLE_PETSC=On \
          -DPETSC_DIR=${HOME}/workspace/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/petsc-3.12.3-jnsim2ouogo4zbyohdtgw3apmfwcwf5s \
          -DUSE_PVODE=On \
          -DUSE_SUNDIALS=On \
          -DENABLE_GTEST_DEATH_TESTS=On \
          -DENABLE_RAJA=On \
          -DENABLE_MPI=On \
          -DENABLE_OPENMP=On \
          -DENABLE_CUDA=Off \
          -DENABLE_HYPRE=On \
          -DHYPRE_DIR="${install_prefix}/hypre" \
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
          $source_dir
fi
