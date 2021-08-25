#!/usr/bin/env bash
# this script is designed to be run in the BOUT source dir
arch=$(uname -m)
compiler=gcc
scratch_dir=`pwd`
build_prefix=${scratch_dir}/build/${arch}-${compiler}/
install_prefix=${scratch_dir}/install/${arch}-${compiler}/
source_prefix=${scratch_dir}/
tpl_prefix=/global/project/projectdirs/bout/BOUT-GPU/tpl/
tpl_install_prefix=${tpl_prefix}/install/${arch}-${compiler}/
module_prefix=/global/project/projectdirs/bout/BOUT-GPU/spack/opt/spack/cray-cnl7-skylake_avx512/gcc-8.3.0_cgpu/

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
    cpp=mpiCC
    fc=mpifort
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
          -DCMAKE_CXX_STANDARD=14 \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_CXX_FLAGS="-w " \
          -DCMAKE_CUDA_FLAGS="-w " \
          -DCMAKE_PREFIX_PATH="${tpl_install_prefix}/raja/share/raja/cmake;${tpl_install_prefix}/umpire/share/umpire/cmake" \
          -DNCXX4_CONFIG:FILEPATH=${module_prefix}/netcdf-cxx4-4.3.1-ptxvbr5iimq3lcapnzs5tw7heniv7mha/bin/ncxx4-config \
          -DNC_CONFIG:FILEPATH=${module_prefix}/netcdf-c-4.7.4-hpuuuxa5vze5qwvqhdzxlpkrigjghgtu/bin/nc-config \
          -DNetCDF_ROOT=${module_prefix}/netcdf-c-4.7.4-hpuuuxa5vze5qwvqhdzxlpkrigjghgtu/ \
          -DUSE_NETCDF=On \
          -DUSE_FFTW=On \
          -DFFTW_ROOT=${module_prefix}/fftw-3.3.8-3nkroqhdwtudny5aifsjujxmzdvdz3jw/ \
          -DUSE_LAPACK=On \
          -DUSE_NLS=On \
          -DENABLE_PETSC=On \
          -DPETSC_DIR=${module_prefix}/petsc-3.13.0-ihyubgkjxq7aldwlwpomlq2wbhkj5lql \
          -DUSE_PVODE=On \
          -DUSE_SUNDIALS=On \
          -DENABLE_GTEST_DEATH_TESTS=On \
          -DENABLE_RAJA=Off \
          -DENABLE_UMPIRE=Off \
          -DENABLE_MPI=On \
          -DENABLE_OPENMP=On \
          -DENABLE_CUDA=Off \
          -DCUDA_ARCH=sm_70 \
          -DCMAKE_CUDA_STANDARD=14 \
          -DBUILD_SHARED_LIBS=Off \
          -DENABLE_GTEST=On \
          -DENABLE_GMOCK=On \
          -DBUILD_GTEST=Off \
          -DBUILD_TESTING=On \
          -DENABLE_TESTS=On \
          -DPACKAGE_TESTS=On \
          -DCHECK=1 \
          -DENABLE_EXAMPLES=OFF \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=On \
          -DCMAKE_VERBOSE_MAKEFILE=Off \
          -Dgtest_disable_pthreads=ON \
          -DCMAKE_INSTALL_RPATH="${module_prefix}/petsc-3.13.0-ihyubgkjxq7aldwlwpomlq2wbhkj5lql/lib;${module_prefix}/hdf5-1.10.6-twbl2egvtk5bmvx4bmlvpnugciapg46s/lib;${module_prefix}/netcdf-cxx4-4.3.1-ptxvbr5iimq3lcapnzs5tw7heniv7mha/lib;${module_prefix}/netcdf-c-4.7.4-hpuuuxa5vze5qwvqhdzxlpkrigjghgtu/lib;${module_prefix}/fftw-3.3.8-3nkroqhdwtudny5aifsjujxmzdvdz3jw/lib;${module_prefix}/sundials-5.1.0-nxrldjqsuekh3nmm4soadzjlwy3ggwz4/lib64" \
          -DCMAKE_BUILD_RPATH="${module_prefix}/petsc-3.13.0-ihyubgkjxq7aldwlwpomlq2wbhkj5lql/lib;${module_prefix}/hdf5-1.10.6-twbl2egvtk5bmvx4bmlvpnugciapg46s/lib;${module_prefix}/netcdf-cxx4-4.3.1-ptxvbr5iimq3lcapnzs5tw7heniv7mha/lib;${module_prefix}/netcdf-c-4.7.4-hpuuuxa5vze5qwvqhdzxlpkrigjghgtu/lib;${module_prefix}/fftw-3.3.8-3nkroqhdwtudny5aifsjujxmzdvdz3jw/lib;${module_prefix}/sundials-5.1.0-nxrldjqsuekh3nmm4soadzjlwy3ggwz4/lib64" \
          $source_dir
fi
