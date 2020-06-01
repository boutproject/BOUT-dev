#!/usr/bin/env bash
arch=$(uname -m)
compiler=gcc
scratch_dir=${HOME}/workspace/BOUT_build_cuda/
build_prefix=${scratch_dir}/build/${arch}-${compiler}/
install_prefix=${scratch_dir}/install/${arch}-${compiler}/
source_prefix=${HOME}/workspace/

# The following modules are provided by the system
module --ignore-cache load spectrum-mpi/rolling-release
module --ignore-cache load cmake/3.14.5
module --ignore-cache load cuda/10.1.243
module --ignore-cache load gcc/8.3.1
#module --ignore-cache load sundials/4.1.0
module --ignore-cache load sundials-5.1.0-gcc-8.3.1-g7h45nh
module --ignore-cache load lapack/3.8.0-gcc-4.9.3

# The following modules are installed by you using spack
module --ignore-cache load openblas-0.3.9-gcc-8.3.1-hx4gart
module --ignore-cache load fftw-3.3.8-gcc-8.3.1-vlusxnt
module --ignore-cache load hdf5-1.10.1-gcc-8.3.1-xkc527f
module --ignore-cache load metis-5.1.0-gcc-8.3.1-wbsy3pr
module --ignore-cache load parmetis-4.0.3-gcc-8.3.1-7lwyokt
module --ignore-cache load superlu-dist-6.1.1-gcc-8.3.1-43evi4b
module --ignore-cache load python-3.7.6-gcc-8.3.1-usivcqa
module --ignore-cache load netcdf-c-4.7.3-gcc-8.3.1-usnrhsd # auto installed as part of netcdf-cxx4 install
module --ignore-cache load zlib-1.2.11-gcc-8.3.1-drkbhfs
module --ignore-cache load netcdf-cxx4-4.3.1-gcc-8.3.1-uj77ss3
#module --ignore-cache load petsc-3.12.3-gcc-8.3.1-ut4eyhs 
#module --ignore-cache load petsc-3.12.3-gcc-8.3.1-jnsim2o
module --ignore-cache load py-setuptools-41.4.0-gcc-8.3.1-d4wih3g
module --ignore-cache load py-cftime-1.0.3.4-gcc-8.3.1-q6ofwn4
module --ignore-cache load py-cython-0.29.14-gcc-8.3.1-5sfsoak
module --ignore-cache load py-pybind11-2.5.0-gcc-8.3.1-4hcy5vc
module --ignore-cache load py-numpy-1.18.2-gcc-8.3.1-6wn32qx
module --ignore-cache load py-scipy-1.4.1-gcc-8.3.1-cck6efe
module --ignore-cache load py-netcdf4-1.4.2-gcc-8.3.1-t6cuidv

# Until we get PETSc linked against cuda-enabled Hypre via autoconf build of PETSc source we rely on two different spack installs of PETSc
# One with Hypre enabled and the other without.  The version of PETSc without Hypre would allow you to link BOUT++ with cuda-enabled Hypre
# Use this script to point to the relevant spack module via PETSC_DIR in the BOUT config


# PETSc spack install signature without Hypre
#spack install petsc@3.12.3%gcc@8.3.1 +fftw +metis +superlu-dist ~hypre +mpi ^hdf5@1.10.1+cxx+hl+mpi+pic+shared ^spectrum-mpi@rolling-release%gcc@8.3.1 ^metis%gcc@8.3.1+real64 ^fftw%gcc@8.3.1

# PETSc install with Hypre
#spack install petsc@3.12.3%gcc@8.3.1 +fftw +metis +superlu-dist +hypre +mpi ^hdf5@1.10.1+cxx+hl+mpi+pic+shared ^spectrum-mpi@rolling-release%gcc@8.3.1 ^metis%gcc@8.3.1+real64 ^fftw%gcc@8.3.1


#netcdf-cxx4 signature
#spack install netcdf-cxx4 ^spectrum-mpi@rolling-release%gcc@8.3.1 ^hdf5@1.10.1+cxx+hl+mpi+pic+shared 

#spack install py-netcdf4  ^hdf5@1.10.1+cxx+hl+mpi+pic+shared ^spectrum-mpi@rolling-release%gcc@8.3.1 ^netcdf-c@4.7.3%gcc@8.3.1

#patchelf --replace-needed libnetcdf.so.7 libnetcdf.so.15 /usr/WS2/holger/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/py-netcdf4-1.4.2-t6cuidvypso7wo73jvka25ea3ltv6rs4/lib/python3.7/site-packages/netCDF4/_netCDF4.cpython-37m-powerpc64le-linux-gnu.so

#patchelf --replace-needed libhdf5.so.8 libhdf5.so.101 /usr/WS2/holger/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/py-netcdf4-1.4.2-t6cuidvypso7wo73jvka25ea3ltv6rs4/lib/python3.7/site-packages/netCDF4/_netCDF4.cpython-37m-powerpc64le-linux-gnu.so 
#patchelf --replace-needed libhdf5_hl.so.8 libhdf5_hl.so.100 /usr/WS2/holger/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/py-netcdf4-1.4.2-t6cuidvypso7wo73jvka25ea3ltv6rs4/lib/python3.7/site-packages/netCDF4/_netCDF4.cpython-37m-powerpc64le-linux-gnu.so 

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

###          --with-netcdf-dir="/usr/WS2/holger/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/netcdf-c-4.7.3-usnrhsddn4n6vko5lvvg63vpwdi25pfg/" \
###          --with-hypre-include="${HOME}/workspace/hypre_automake/install/include" \
###          --with-hypre-lib=[${HOME}/workspace/hypre_automake/install/lib/libHYPRE.a,libHYPRE.a] \
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
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_BUILD_TYPE=RelWithDebInfo \
          -DCMAKE_PREFIX_PATH="${install_prefix}/raja/share/raja/cmake" \
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
elif [ "$pkg" == "petsc" ]; then
    source_dir=${source_prefix}/${pkg}
    cd $source_dir
    mkdir -p $install_dir
    ./configure --prefix=${install_dir} \
          --CXX=$cpp \
          --CC=$cc \
          --FC=$fc \
          --with-cuda=0 \
          --with-blaslapack-dir="${LAPACK_DIR}/.." \
          --with-openblas-dir="/usr/WS2/holger/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/openblas-0.3.9-hx4garttzvdqltwiclgkfguqfpmkjnyi/" \
          --with-hdf5-dir="/usr/WS2/holger/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/hdf5-1.10.1-xkc527ftnhfp2zz3j4v7vo5l7ypse6m3/" \
          --download-netcdf \
          --download-mumps \
          --download-hypre \
          --download-scalapack \
          --with-metis-dir="/usr/WS2/holger/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/metis-5.1.0-wbsy3praw43mfeyat3x7jz2dym3bifsl/" \
          --with-parmetis-dir="/usr/WS2/holger/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/parmetis-4.0.3-7lwyokt7zd2cbvxl2mglmoho5xmnze4c/" \
          --with-superlu_dist-dir="/usr/WS2/holger/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/superlu-dist-6.1.1-43evi4baw2e7hkrgecwvtspvomo4wltw/" \
          --with-fftw-dir="/usr/WS2/holger/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/fftw-3.3.8-vlusxntezsqbongx3pyjwlfg6ddjbbxw/"  \
          --with-zlib-dir="/usr/WS2/holger/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/zlib-1.2.11-drkbhfszh6qdkeofmyrnbil3kgfnfz4s/" 
elif [ "$pkg" == "BOUT-dev" ]; then
    source_dir=${source_prefix}/${pkg}
    cd $source_dir && {
        branch=$(git rev-parse --abbrev-ref HEAD)
        [[ "$branch" != "next-outerloop-blt" ]] && exit 1
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
          -DENABLE_PETSC=On \
          -DPETSC_DIR=${HOME}/workspace/spack/opt/spack/linux-rhel7-power9le/gcc-8.3.1/petsc-3.12.3-jnsim2ouogo4zbyohdtgw3apmfwcwf5s \
          -DUSE_PVODE=On \
          -DUSE_SUNDIALS=On \
          -DENABLE_GTEST_DEATH_TESTS=On \
          -DENABLE_RAJA=On \
          -DENABLE_UMPIRE=On \
          -DENABLE_MPI=On \
          -DENABLE_OPENMP=Off \
          -DENABLE_CUDA=On \
          -DCUDA_ARCH=sm_70 \
          -DENABLE_HYPRE=Off \
          -DHYPRE_DIR="${HOME}/workspace/hypre_automake/install" \
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
