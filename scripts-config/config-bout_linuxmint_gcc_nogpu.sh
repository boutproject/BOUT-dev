#!/bin/bash -i

arch=$(uname -m)
compiler=gcc
scratch_dir=${HOME}/workspace/BOUT_build/
build_prefix=${scratch_dir}/build/${arch}-${compiler}/
install_prefix=${scratch_dir}/install/${arch}-${compiler}/
source_prefix=${HOME}/workspace/

module load gmake-4.2.1-gcc-7.3.0-l7i7nma
module load gcc-7.3.0-gcc-7.3.0-zhmzu4d
module load cmake-3.16.2-gcc-7.3.0-u54r6tn
module load openmpi-3.1.5-gcc-7.3.0-354646n
module load netcdf-cxx4-4.3.1-gcc-7.3.0-ggnkvuj
module load fftw-3.3.8-gcc-7.3.0-yxqd4gn
module load openblas-0.3.8-gcc-7.3.0-xmnfitp
module load hdf5-1.10.1-gcc-7.3.0-miduohb
module load python-3.7.6-gcc-7.3.0-cfhy5mb
module load py-cftime-1.0.3.4-gcc-7.3.0-5ow6zpj
module load py-cython-0.29.14-gcc-7.3.0-4wmiozn
module load py-netcdf4-1.4.2-gcc-7.3.0-v55waqr
module load py-numpy-1.18.1-gcc-7.3.0-ic45rrw
module load py-scipy-1.4.1-gcc-7.3.0-773vnfv
module load sundials-5.1.0-gcc-7.3.0-jawjhwl
module load hypre-2.18.2-gcc-7.3.0-33inqou
module load superlu-dist-6.1.1-gcc-7.3.0-jok2jvn
module load petsc-3.12.3-gcc-7.3.0-s77fwdj

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
        [[ "$branch" != "feature/uvm-raja-next" ]] && exit 1
    }
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=Debug \
          -DCMAKE_PREFIX_PATH="${install_prefix}/raja/share/raja/cmake;${install_prefix}/HDF5/share/cmake;${install_prefix}/umpire/share/umpire/cmake" \
          -DUSE_NETCDF=On \
          -DNCXX4_CONFIG:FILEPATH=${HOME}/spack/opt/spack/linux-linuxmint19-haswell/gcc-7.3.0/netcdf-cxx4-4.3.1-ggnkvuj77tr2yoq4vaa5xtnfattqavq4/bin/ncxx4-config \
          -DNC_CONFIG:FILEPATH=/usr/bin/nc-config \
          -DUSE_FFTW=On \
          -DUSE_LAPACK=On \
          -DUSE_NLS=On \
          -DUSE_PETSC=On \
          -DUSE_PVODE=On \
          -DPETSC_DIR=${HOME}/spack/opt/spack/linux-linuxmint19-haswell/gcc-7.3.0/petsc-3.12.3-s77fwdjblxkslbp2ymv5uehkmimdbyka \
          -DUSE_SUNDIALS=On \
          -DSUNDIALS_ROOT=${HOME}/spack/opt/spack/linux-linuxmint19-haswell/gcc-7.3.0/sundials-5.1.0-jawjhwlqup3qlusyygjoj4r4xxw7qvlf \
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
