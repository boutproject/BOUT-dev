#!/usr/bin/env tcsh
set arch=`uname -m`
set compiler=gcc
set env_prefix=/usr/WS2/BOUT-GPU/lassen/env/
source ${env_prefix}/scripts/setup-lassen-gpu.csh
set base=${env_prefix}/tpl
set build_prefix=${base}/build/${arch}-${compiler}/
set install_prefix=${base}/install/${arch}-${compiler}/
set source_prefix=${base}

if ( $#argv > 0 ) then
   if ( X$1 =~ X"-c" ) then
      set pkg=$2
      set dir=${build_prefix}/$pkg
      #rm -rf $dir
      echo "cleanup of build dir" ${dir}
   else
      set pkg=$1
   endif
else
   echo "This script needs to specify a package argument : e.g. raja or umpire"
   exit 1
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
if ( "$pkg" == "raja" ) then
    echo "enter raja script"
    set source_dir=${source_prefix}/${pkg}
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
else if ( "$pkg" == "umpire" ) then
    echo "enter umpire script"
    set source_dir=${source_prefix}/Umpire
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_CXX_STANDARD=14 \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=Release \
          -DENABLE_CUDA=On \
          -DCUDA_ARCH=sm_70 \
          -DENABLE_TESTS=On \
          -DENABLE_TOOLS=On \
          -DENABLE_NUMA=On \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=On \
          -DCMAKE_VERBOSE_MAKEFILE=On \
          $source_dir
else
   echo "Package must be one of raja or umpire"
   exit 1
endif
