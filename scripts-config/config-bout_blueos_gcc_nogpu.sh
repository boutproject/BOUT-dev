#!/usr/bin/env bash
arch=$(uname -m)
compiler=gcc
scratch_dir=${HOME}/workspace/BOUT_build_serial/
build_prefix=${scratch_dir}/build/${arch}-${compiler}/
install_prefix=${scratch_dir}/install/${arch}-${compiler}/
source_prefix=${HOME}/workspace/

# The following modules are provided by the system
module --ignore-cache load spectrum-mpi/rolling-release
module --ignore-cache load cmake/3.14.5
module --ignore-cache load cuda/10.1.243
module --ignore-cache load gcc/8.3.1
module --ignore-cache load sundials/4.1.0
module --ignore-cache load lapack/3.8.0-gcc-4.9.3

# The following modules are installed by you using spack
module --ignore-cache load fftw-3.3.8-gcc-8.3.1-vlusxnt
module --ignore-cache load hdf5-1.10.1-gcc-8.3.1-xkc527f
module --ignore-cache load python-3.7.6-gcc-8.3.1-usivcqa
module --ignore-cache load netcdf-c-4.7.3-gcc-8.3.1-usnrhsd # auto installed as part of netcdf-cxx4 install
module --ignore-cache load netcdf-cxx4-4.3.1-gcc-8.3.1-uj77ss3
module --ignore-cache load petsc-3.12.3-gcc-8.3.1-ut4eyhs
module --ignore-cache load py-setuptools-41.4.0-gcc-8.3.1-d4wih3g
module --ignore-cache load py-cftime-1.0.3.4-gcc-8.3.1-q6ofwn4
module --ignore-cache load py-cython-0.29.14-gcc-8.3.1-5sfsoak
module --ignore-cache load py-pybind11-2.5.0-gcc-8.3.1-4hcy5vc
module --ignore-cache load py-numpy-1.18.2-gcc-8.3.1-6wn32qx
module --ignore-cache load py-scipy-1.4.1-gcc-8.3.1-cck6efe
module --ignore-cache load py-netcdf4-1.4.2-gcc-8.3.1-t6cuidv

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
          -DENABLE_CUDA=Off \
          -DENABLE_TESTS=ON \
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
          -DENABLE_OPENMP=On \
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
          -DCMAKE_BUILD_TYPE=RelWithDebInfo \
          -DCMAKE_PREFIX_PATH="${install_prefix}/raja/share/raja/cmake;${install_prefix}/HDF5/share/cmake;${install_prefix}/umpire/share/umpire/cmake" \
          -DUSE_NETCDF=On \
          -DNCXX4_CONFIG:FILEPATH=${HOME}/workspace/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/netcdf-cxx4-4.3.1-uj77ss3iyalvvd3mesfvdv5pihxafxap/bin/ncxx4-config \
          -DNC_CONFIG:FILEPATH=${HOME}/workspace/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/netcdf-c-4.7.3-usnrhsddn4n6vko5lvvg63vpwdi25pfg/bin/nc-config \
          -DUSE_FFTW=On \
          -DUSE_LAPACK=On \
          -DUSE_NLS=On \
          -DENABLE_PETSC=Off \
          -DPETSC_DIR=${HOME}/workspace/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/petsc-3.12.3-ut4eyhszto3njzigfesmrfqnbhegp7iu \
          -DUSE_PVODE=On \
          -DUSE_SUNDIALS=On \
          -DENABLE_GTEST_DEATH_TESTS=On \
          -DENABLE_RAJA=Off \
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
          -DCMAKE_INSTALL_RPATH="${HOME}/workspace/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/petsc-3.12.3-ut4eyhszto3njzigfesmrfqnbhegp7iu/lib;${HOME}/workspace/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/hdf5-1.10.1-xkc527ftnhfp2zz3j4v7vo5l7ypse6m3/lib" \
          -DCMAKE_BUILD_RPATH="${HOME}/workspace/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/petsc-3.12.3-ut4eyhszto3njzigfesmrfqnbhegp7iu/lib;${HOME}/workspace/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/hdf5-1.10.1-xkc527ftnhfp2zz3j4v7vo5l7ypse6m3/lib" \
          $source_dir
fi
