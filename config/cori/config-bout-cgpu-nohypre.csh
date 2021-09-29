#!/usr/bin/env tcsh
# this script is designed to be run in the BOUT source dir
set arch = `uname -m`
set compiler =  gcc
set scratch_dir = `pwd`
set build_prefix = ${scratch_dir}/build/${arch}-${compiler}/
set install_prefix = ${scratch_dir}/install/${arch}-${compiler}/
set source_prefix = ${scratch_dir}/
set tpl_prefix = /global/project/projectdirs/bout/BOUT-GPU/tpl/
set tpl_install_prefix = ${tpl_prefix}/install/${arch}-${compiler}/
set module_prefix = /global/project/projectdirs/bout/BOUT-GPU/spack/opt/spack/cray-cnl7-skylake_avx512/gcc-8.3.0_cgpu/

set pkg = BOUT-dev
echo 'continue with script ' $pkg

if ( "$1" == "-c" ) then
    set pkg = $2
    set dir = ${build_prefix}/$pkg
    echo Removing $dir ...
    rm -rf $dir
endif

if ( "$compiler" == "xl" ) then
    set cc = xlc
    set cpp = xlC
    set fc = xlf
    set extra_cflags = ""
else if ( "$compiler" == "clang" ) then
    set cc = clang
    set cpp = clang++
#    set fc = xlf
    set fc = mpifort
    set = extra_cflags ""
else if ( "$compiler" == "gcc" ) then
    set cc = mpicc
    set cpp = mpiCC
    set fc = mpifort
    set extra_cflags = ""
else
    exit 1
endif

set build_dir = ${build_prefix}/${pkg}
set install_dir = ${install_prefix}/${pkg}
set source_dir = ${source_prefix} 
echo 'Build in ' ${build_dir}
echo 'Install in ' ${install_dir}
echo 'Source in ' ${source_dir}

if ( "$pkg" == "BOUT-dev" ) then
    echo 'enter BOUT-dev script'
    set source_dir = ${source_prefix} 
    mkdir -p $build_dir && cd $build_dir
    cmake  \
          -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_CXX_STANDARD=14 \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_CXX_FLAGS="-w " \
          -DCMAKE_CUDA_FLAGS="-w " \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_PREFIX_PATH="${tpl_install_prefix}/raja/share/raja/cmake;${tpl_install_prefix}/umpire/share/umpire/cmake" \
          -DNCXX4_CONFIG:FILEPATH=${module_prefix}/netcdf-cxx4-4.3.1-ptxvbr5iimq3lcapnzs5tw7heniv7mha/bin/ncxx4-config \
          -DNC_CONFIG:FILEPATH=${module_prefix}/netcdf-c-4.7.4-hpuuuxa5vze5qwvqhdzxlpkrigjghgtu/bin/nc-config \
          -DNetCDF_ROOT=${module_prefix}/netcdf-c-4.7.4-hpuuuxa5vze5qwvqhdzxlpkrigjghgtu/ \
          -DBOUT_USE_NETCDF=On \
          -DBOUT_USE_FFTW=On \
          -DFFTW_ROOT=${module_prefix}/fftw-3.3.8-3nkroqhdwtudny5aifsjujxmzdvdz3jw/ \
          -DBOUT_USE_LAPACK=On \
          -DBOUT_USE_NLS=On \
          -DBOUT_USE_PETSC=On \
          -DPETSC_DIR=${module_prefix}/petsc-3.13.0-ym7gwgxfutg4m7ap3bz5tbvsitvfq2w2 \
          -DBOUT_USE_PVODE=On \
          -DBOUT_USE_SUNDIALS=On \
          -DBOUT_ENABLE_RAJA=On \
          -DBOUT_ENABLE_UMPIRE=On \
          -DBOUT_ENABLE_MPI=On \
          -DBOUT_ENABLE_OPENMP=Off \
	  -DBOUT_ENABLE_WARNINGS=Off \
          -DBOUT_ENABLE_CUDA=On \
          -DCUDA_ARCH="compute_70,code=sm_70" \
          -DCMAKE_CUDA_STANDARD=14 \
          -DBOUT_USE_HYPRE=Off \
          -DHYPRE_DIR="${tpl_prefix}/hypre_dir/hypre_autoconf/install" \
          -DHYPRE_CUDA=Off \
          -DBUILD_SHARED_LIBS=Off \
          -DBOUT_TESTS=On \
          -DCHECK=1 \
          -DBOUT_BUILD_EXAMPLES=OFF \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=On \
          -DCMAKE_VERBOSE_MAKEFILE=On \
          -Dgtest_disable_pthreads=On \
          -DCMAKE_INSTALL_RPATH="${module_prefix}/petsc-3.13.0-ym7gwgxfutg4m7ap3bz5tbvsitvfq2w2/lib;${module_prefix}/hdf5-1.10.6-twbl2egvtk5bmvx4bmlvpnugciapg46s/lib;${module_prefix}/netcdf-cxx4-4.3.1-ptxvbr5iimq3lcapnzs5tw7heniv7mha/lib;${module_prefix}/netcdf-c-4.7.4-hpuuuxa5vze5qwvqhdzxlpkrigjghgtu/lib;${module_prefix}/fftw-3.3.8-3nkroqhdwtudny5aifsjujxmzdvdz3jw/lib;${module_prefix}/sundials-5.1.0-nxrldjqsuekh3nmm4soadzjlwy3ggwz4/lib64" \
          -DCMAKE_BUILD_RPATH="${module_prefix}/petsc-3.13.0-ym7gwgxfutg4m7ap3bz5tbvsitvfq2w2/lib;${module_prefix}/hdf5-1.10.6-twbl2egvtk5bmvx4bmlvpnugciapg46s/lib;${module_prefix}/netcdf-cxx4-4.3.1-ptxvbr5iimq3lcapnzs5tw7heniv7mha/lib;${module_prefix}/netcdf-c-4.7.4-hpuuuxa5vze5qwvqhdzxlpkrigjghgtu/lib;${module_prefix}/fftw-3.3.8-3nkroqhdwtudny5aifsjujxmzdvdz3jw/lib;${module_prefix}/sundials-5.1.0-nxrldjqsuekh3nmm4soadzjlwy3ggwz4/lib64" \
          $source_dir
endif
