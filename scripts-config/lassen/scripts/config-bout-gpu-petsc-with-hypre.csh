#!/usr/bin/env tcsh

set arch=`uname -m`
set compiler=gcc

set myenv_prefix=/usr/workspace/BOUT-GPU/lassen/Xu/env/
set hjenv_prefix=/usr/WS2/BOUT-GPU/lassen/env/

source ${myenv_prefix}/scripts/setup-lassen-gpu-petsc_with_hypre.csh
set build_base_dir=${myenv_prefix}/BOUT_build_cuda_ph
set build_prefix=${build_base_dir}/build/${arch}-${compiler}/
set install_prefix=${build_base_dir}/install/${arch}-${compiler}/
set source_prefix=${myenv_prefix}/BOUT-dev

set tpl_prefix=${hjenv_prefix}/tpl/
set tpl_install_prefix=${tpl_prefix}/install/${arch}-${compiler}/

set hjtpl_prefix=${hjenv_prefix}/tpl/
set module_prefix=${hjenv_prefix}/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/

set pkg=BOUT-dev
echo 'continue with script ' $pkg
if ( $#argv > 0 ) then
   if ( X$1 =~ X"-c" ) then
      set dir=${build_prefix}/$pkg
      rm -rf $dir
      echo "cleanup of build dir" ${dir}
   endif
endif 
echo "arch=${arch}" 

if ( $compiler == "xl" ) then
    set cc=xlc
    set cpp=xlC
    set fc=xlf
    set extra_cflags=""
else if ( $compiler == "clang" ) then 
    set cc=clang
    set cpp=clang++
    set fc=mpifort
    set extra_cflags=""
else if ( $compiler == "gcc" ) then
    set cc=mpicc
    set cpp=mpiCC
    set fc=mpifort
    set extra_cflags=""
else
    exit 1
endif


set build_dir=${build_prefix}/${pkg}
set install_dir=${install_prefix}/${pkg}
echo 'Build in ' ${build_dir}
echo 'Install in ' ${install_dir}

if ( "$pkg" == "BOUT-dev" ) then
    echo 'enter BOUT-dev script'
    set source_dir=${source_prefix} 
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
endif
