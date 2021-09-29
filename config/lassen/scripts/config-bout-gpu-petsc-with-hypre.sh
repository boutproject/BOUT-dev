#!/usr/bin/env bash
arch=$(uname -m)
compiler=gcc
env_prefix=/usr/WS2/BOUT-GPU/lassen/env/
build_base_dir=${env_prefix}/BOUT_build_cuda_ph
build_prefix=${build_base_dir}/build/${arch}-${compiler}/
install_prefix=${build_base_dir}/install/${arch}-${compiler}/
source_prefix=${env_prefix}/BOUT-dev
tpl_prefix=${env_prefix}/tpl/
tpl_install_prefix=${tpl_prefix}/install/${arch}-${compiler}/
module_prefix=/usr/WS2/BOUT-GPU/lassen/env/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/

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
echo 'Build in ' ${build_dir}
echo 'Install in ' ${install_dir}

if [ "$pkg" == "BOUT-dev" ]; then
    echo 'enter BOUT-dev script'
    source_dir=${source_prefix} 
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_CXX_STANDARD=14 \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_PREFIX_PATH="${tpl_install_prefix}/raja/share/raja/cmake;${tpl_install_prefix}/umpire/share/umpire/cmake" \
          -DNCXX4_CONFIG:FILEPATH=${module_prefix}/netcdf-cxx4-4.3.1-sj65c6a4kyaaw3d6dopj6txamfkli3dc/bin/ncxx4-config \
          -DNC_CONFIG:FILEPATH=${module_prefix}/netcdf-c-4.7.4-7x3quj36clrfnejvxqyb5ky6gbco73zr/bin/nc-config \
          -DUSE_NETCDF=On \
          -DUSE_FFTW=On \
          -DUSE_LAPACK=On \
          -DUSE_NLS=On \
          -DENABLE_PETSC=On \
          -DPETSC_DIR=${module_prefix}/petsc-3.13.0-w3s32s2mhufbunrx7qu3i6pjm5h2f7g5/ \
          -DUSE_PVODE=On \
          -DUSE_SUNDIALS=On \
          -DENABLE_RAJA=On \
          -DENABLE_UMPIRE=Off \
          -DENABLE_MPI=On \
          -DENABLE_OPENMP=Off \
          -DENABLE_CUDA=On \
          -DCUDA_ARCH=sm_70 \
          -DCMAKE_CUDA_STANDARD=14 \
          -DBUILD_SHARED_LIBS=Off \
          -DENABLE_GTEST_DEATH_TESTS=On \
          -DENABLE_GTEST=On \
          -DENABLE_GMOCK=On \
          -DBUILD_GTEST=Off \
          -DBUILD_TESTING=On \
          -DENABLE_TESTS=On \
          -DPACKAGE_TESTS=ON \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=On \
          -DCMAKE_VERBOSE_MAKEFILE=On \
          -Dgtest_disable_pthreads=ON \
          -DCMAKE_INSTALL_RPATH="${module_prefix}/petsc-3.13.0-w3s32s2mhufbunrx7qu3i6pjm5h2f7g5/lib;${module_prefix}/hdf5-1.10.6-zkwocrtngfhf5nw6xk5c4cbekrdycznu/lib" \
          -DCMAKE_BUILD_RPATH="${module_prefix}/petsc-3.13.0-w3s32s2mhufbunrx7qu3i6pjm5h2f7g5/lib;${module_prefix}/hdf5-1.10.6-zkwocrtngfhf5nw6xk5c4cbekrdycznu/lib" \
          $source_dir
fi
